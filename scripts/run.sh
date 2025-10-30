#!/usr/bin/env bash
set -euo pipefail

if [ -z "${OPENAI_API_KEY:-}" ]; then
  echo "ERROR: OPENAI_API_KEY is not set."
  exit 1
fi

echo "=== Step 1: Differential Expression (Surface Targets) ==="
Rscript R/01_Surface_targets.R

echo "=== Step 2: Drugibility (DGIdb, STRING, ChEMBL) ==="
Rscript R/02_Drugibility.R

echo "=== Step 3: Survival and Dependency Modeling ==="
Rscript R/03_Survival_and_dependency.R

echo "=== Step 4: ClinicalTrials.gov Summary ==="
Rscript R/04_Clinical_trials_summary.R

echo "=== Step 5: Novelty on Clinical Trials (baseline) ==="
Rscript R/05_Novelty_on_clinical_trials.R

echo "=== Step 6: PubMed Retrieval (R) ==="
Rscript R/06_PubMed_retrieval.R

echo "=== Step 7: LLM-based Literature Novelty (Python) ==="
python src/run_llm_novelty.py

echo "=== Step 8: Target Scoring & Ranking ==="
Rscript R/07_Target_scoring_and_ranking.R

echo "=== Step 9: Safety Adjustment & Final Ranking ==="
Rscript R/08_Ranking_with_safety.R

echo "✅ Pipeline complete. See reports/ for outputs."
