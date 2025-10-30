# 🧬 Gene-NLP Pipeline

A hybrid **R + Python** pipeline for candidate-target discovery and literature-based novelty scoring.  
The workflow integrates multi-omics data, druggability information, survival and dependency modeling, clinical-trial mining, PubMed retrieval, and LLM-based text analysis.

---

## 📦 Features

- R-based statistical and data-integration modules  
- Python-based LLM novelty analysis using the OpenAI API  
- Conda environment for reproducibility  
- Docker image for one-command execution  
- **Git LFS** management of large datasets

---

## ✅ Quick Start Summary
```
git lfs install
git clone https://github.com/1boChen/gene-nlp.git
cd gene-nlp
cp .env.example .env   # add your OpenAI key
docker build -t gene-nlp .
docker run --rm \
  -e OPENAI_API_KEY=$(grep OPENAI_API_KEY .env | cut -d '=' -f2) \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/reports:/app/reports" \
  gene-nlp
```
  
---

## ⚙️ 1 · Prerequisites

Choose **one** setup method:

### 🐋 Docker (recommended)

- Install [Docker Desktop](https://docs.docker.com/get-docker/).

### 🐍 Conda (for local dev)

- Install [Miniconda / Anaconda](https://docs.conda.io/en/latest/miniconda.html).  
  Then:

  ```bash
  conda env create -f env/environment.yml
  conda activate gene_nlp
  ```
---

## 🧰 2 · Clone the Repository (with Git LFS)

This repo uses Git Large File Storage
 for big data files.
Install and initialize LFS before cloning or pulling:

```
git lfs install
git clone https://github.com/1boChen/gene-nlp.git
cd gene-nlp
git lfs pull
```

If you skip git lfs install, large files will appear as tiny pointer text files instead of real data.

---

## 🔑 3 · Configure the OpenAI API Key

Copy the example file and add your key:
```
cp .env.example .env
```

Edit .env:
```
OPENAI_API_KEY=sk-your-key-here
```

---

## 🚀 4 · Run with Docker

Build the image
```
docker build -t gene-nlp:latest .
```

Execute the pipeline
```
docker run --rm \
  -e OPENAI_API_KEY=$(grep OPENAI_API_KEY .env | cut -d '=' -f2) \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/external:/app/external" \
  -v "$(pwd)/reports:/app/reports" \
  gene-nlp:latest
```

All results will appear in reports/.

---

## 🧩 5 · Run Manually (Conda Environment)

```
# activate environment
conda activate gene_nlp
# export key
export $(grep -v '^#' .env | xargs)
# run
bash scripts/run.sh
```

---

## 📄 6 · Input Data

Place the following in data/:
- expression.csv:	Gene × sample matrix (gene in first column)
- phenotype.csv	:Sample metadata (sample_id, response)

---

🧠 7 · Outputs

All intermediate and final tables are written to reports/, including:

Differential-expression results

Druggability & STRING summaries

Survival / dependency models

Clinical-trial summaries

PubMed abstracts + LLM novelty scores

Final target ranking

