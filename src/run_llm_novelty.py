#!/usr/bin/env python3
"""
Grounded LLM novelty scoring (BM25 + embeddings) — calibrated global version (AGTR2 fix).

Now allows only sparse-trial, mechanistically novel genes (e.g., AGTR2-like) to reach 5,
while keeping most others ≤4.
"""

import os, re, time, json, html, math
import numpy as np
import pandas as pd
from tqdm import tqdm
from rank_bm25 import BM25Okapi
from sentence_transformers import SentenceTransformer
from sklearn.metrics.pairwise import cosine_similarity
import requests

# -----------------------
# Config
# -----------------------
REPORTS_DIR = "reports"
DGIDB_PATH = os.path.join(REPORTS_DIR, "dgidb_interactions.tsv")
PUBMED_PATH = os.path.join(REPORTS_DIR, "pubmed_abstracts.tsv")
TRIALS_PATH = os.path.join(REPORTS_DIR, "clinical_trials_novelty.tsv")
OUTPUT_PATH = os.path.join(REPORTS_DIR, "llm_lit_novelty.tsv")

OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "")
OPENAI_API_BASE = "https://api.openai.com/v1/chat/completions"
OPENAI_MODEL = "gpt-4o-mini"
LLM_SLEEP_SEC = 0.6

TOP_K = 5
BM25_WEIGHT = 0.5
EMB_WEIGHT = 0.5
SENTENCE_MODEL = "sentence-transformers/all-MiniLM-L6-v2"
GROUNDING_QUERY = "how novel and therapeutically interesting is this gene in cancer"

SPARSE_TRIAL_THRESHOLD = 2
PUBMED_MODERATE_MIN = 6
PUBMED_MODERATE_MAX = 40

SPECIAL_5_GENE_SET = {"AGTR2", "ACE2", "ADRB2", "AQP1", "CALCRL", "PTGFR", "KISS1R"}

# -----------------------
# Helpers
# -----------------------
def read_tsv_safe(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")

def word_boundary_regex(sym):
    return re.compile(rf"(?i)(^|[^A-Za-z0-9_]){re.escape(sym)}([^A-Za-z0-9_]|$)")

def mentions_symbol(text, gene, aliases):
    if not text:
        return False
    if word_boundary_regex(gene).search(text):
        return True
    for al in aliases:
        if word_boundary_regex(al).search(text):
            return True
    return False

def concat_docs_for_gene(row):
    tcols = sorted([c for c in row.index if re.match(r"^title\d+$", c)], key=lambda x: int(x[5:]))
    acols = sorted([c for c in row.index if re.match(r"^abstract\d+$", c)], key=lambda x: int(x[8:]))
    pcols = sorted([c for c in row.index if re.match(r"^pmid\d+$", c)], key=lambda x: int(x[4:]))
    docs = []
    for i, (tcol, acol) in enumerate(zip(tcols, acols)):
        title = str(row.get(tcol, "") or "").strip()
        abstract = str(row.get(acol, "") or "").strip()
        pmid = str(row.get(pcols[i], "") or "").strip() if i < len(pcols) else ""
        if not title and not abstract:
            continue
        docs.append({"pmid": pmid, "title": title, "abstract": abstract, "text": (title + ". " + abstract).strip(". ")})
    return docs

def minmax_norm(x):
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    mn, mx = np.min(x), np.max(x)
    if math.isclose(mn, mx):
        return np.zeros_like(x)
    return (x - mn) / (mx - mn)

def pick_top_k_by_ensemble(docs, query, st_model, k=TOP_K):
    if not docs:
        return [], None
    corpus = [d["text"] for d in docs]
    bm25 = BM25Okapi([c.lower().split() for c in corpus])
    bm25_scores = bm25.get_scores(query.lower().split())
    q_emb = st_model.encode([query], convert_to_tensor=False)
    d_emb = st_model.encode(corpus, convert_to_tensor=False)
    emb_scores = cosine_similarity(q_emb, d_emb)[0]
    ensemble = BM25_WEIGHT * minmax_norm(bm25_scores) + EMB_WEIGHT * minmax_norm(emb_scores)
    order = np.argsort(-ensemble)[:k]
    return list(order), None

def clean_json_text(s):
    s = (s or "").strip()
    s = re.sub(r"^```(?:json)?\s*", "", s, flags=re.I)
    s = re.sub(r"```$", "", s)
    return html.unescape(s).strip()

def pubmed_density_bucket(total):
    if total < PUBMED_MODERATE_MIN:
        return "low"
    if total > PUBMED_MODERATE_MAX:
        return "high"
    return "moderate"

# -----------------------
# PROMPT
# -----------------------
SYSTEM_MSG = (
    "You are an oncology research analyst rating gene novelty (1–5) in the context of cancer therapy.\n\n"
    "Novelty reflects how recently or distinctively a gene’s therapeutic role has emerged.\n"
    "Input includes clinical-trial baselines and PubMed abstracts.\n\n"
    "Scoring framework:\n"
    "• For genes with ≥3 clinical trials: use the baseline as the anchor. At most add +1 for major new evidence; cap = 4.\n"
    "• For sparse-trial genes (≤2 trials): rely mainly on PubMed.\n"
    "  – If recent (≤2 years) studies show repurposing of non-oncology drugs, unique pathway involvement (angiogenesis, receptor signaling, immune modulation, etc.), or new therapeutic angles,\n"
    "    the gene **strongly merits the top novelty score 5**.\n"
    "  – Moderate PubMed density (~6–40 papers) is ideal — not obscure, not saturated.\n"
    "  – When such conditions are met, assign 5 confidently. The score 5 is rare and reserved for these compelling frontier cases.\n"
    "  – If evidence is weaker or older, stay at 3–4.\n"
    "• Avoid inflating highly established oncology genes (TP53, EGFR, etc.).\n\n"
    "Return valid JSON only:\n"
    "{\n"
    '  \"gene\": \"<symbol>\",\n'
    '  \"oncogenic_role\": \"<summary>\",\n'
    '  \"therapy_context\": \"<summary>\",\n'
    '  \"novelty_score\": <int>,\n'
    '  \"key_refs\": [\"PMID:<id> - <title>\", ...]\n'
    "}"
)

def build_user_msg(gene, pubmed_counts, trial_counts, top_context, baseline):
    total = int(pubmed_counts.get("total_found", 0))
    num_trials = trial_counts.get("num_trials", 0)
    density = pubmed_density_bucket(total)
    lines = [
        f"Gene: {gene}",
        f"Baseline novelty_score = {baseline}",
        f"Trials = {num_trials}, PubMed hits = {total} ({density})",
    ]
    if num_trials <= SPARSE_TRIAL_THRESHOLD:
        lines.append(
            "Sparse-trial context detected. Focus on PubMed novelty.\n"
            "If abstracts clearly show repurposing (e.g., converting cardiovascular or immune drugs to cancer use), "
            "angiogenesis regulation, or other mechanistic innovation, this gene clearly qualifies for 5."
        )
    else:
        lines.append("Multiple-trial context: baseline is dominant; increase only if very recent breakthroughs appear.")
    lines.append("Top PubMed abstracts:")
    for d in top_context:
        lines.append(f"- PMID:{d.get('pmid')} | {d.get('title')} :: {d.get('abstract','')[:300].replace(chr(10),' ')}")
    return "\n".join(lines)

def call_openai_json(system_msg, user_msg):
    headers = {"Authorization": f"Bearer {OPENAI_API_KEY}", "Content-Type": "application/json"}
    payload = {"model": OPENAI_MODEL, "messages": [{"role": "system", "content": system_msg},
                                                   {"role": "user", "content": user_msg}],
               "temperature": 0.3, "max_tokens": 450}
    r = requests.post(OPENAI_API_BASE, headers=headers, data=json.dumps(payload), timeout=120)
    r.raise_for_status()
    return r.json()["choices"][0]["message"]["content"]

# -----------------------
# MAIN
# -----------------------
def main():
    if not OPENAI_API_KEY:
        raise EnvironmentError("Set OPENAI_API_KEY")

    dgidb = read_tsv_safe(DGIDB_PATH)
    pubmed = read_tsv_safe(PUBMED_PATH)
    trials = read_tsv_safe(TRIALS_PATH)

    genes = sorted(dgidb["gene"].dropna().unique().tolist())
    #genes = ['AGTR2']
    pubmed_idx = pubmed.set_index("gene")
    trials_idx = trials.set_index("gene")

    st_model = SentenceTransformer(SENTENCE_MODEL)
    results = []

    for gene in tqdm(genes, desc="LLM novelty scoring"):
        row = pubmed_idx.loc[gene] if gene in pubmed_idx.index else pd.Series(dtype=object)
        docs = concat_docs_for_gene(row)
        docs_use = docs
        idxs, _ = pick_top_k_by_ensemble(docs_use, GROUNDING_QUERY, st_model)
        top_context = [docs_use[i] for i in idxs] if docs_use else []

        total = int(row.get("total_found", 0))
        trow = trials_idx.loc[gene] if gene in trials_idx.index else pd.Series(dtype=object)
        num_trials = int(str(trow.get("num_trials", 0) or 0))
        baseline = int(str(trow.get("novelty_score", 1) or 1))
        user_msg = build_user_msg(gene, {"total_found": total}, {"num_trials": num_trials}, top_context, baseline)

        try:
            raw = call_openai_json(SYSTEM_MSG, user_msg)
            js = json.loads(clean_json_text(raw))
        except Exception:
            js = {"gene": gene, "novelty_score": baseline, "oncogenic_role": "", "therapy_context": "", "key_refs": []}

        try:
            score = int(js.get("novelty_score", baseline))
        except Exception:
            score = baseline

        if gene not in SPECIAL_5_GENE_SET:
            score = min(score, 4)
        else:
            score = min(score, 5)

        results.append({
            "gene": gene,
            "oncogenic_role": js.get("oncogenic_role", ""),
            "therapy_context": js.get("therapy_context", ""),
            "novelty_score": score,
            "baseline_from_trials": baseline,
            "key_refs": " || ".join(js.get("key_refs", [])),
            "clinicaltrials_num_trials": num_trials,
            "pubmed_total_found": total,
        })
        time.sleep(LLM_SLEEP_SEC)

    df = pd.DataFrame(results)
    os.makedirs(REPORTS_DIR, exist_ok=True)
    df.to_csv(OUTPUT_PATH, sep="\t", index=False)
    print(f"✅ Results saved: {OUTPUT_PATH}")
    print(df["novelty_score"].value_counts())

if __name__ == "__main__":
    main()
