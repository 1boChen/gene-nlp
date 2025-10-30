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

## 🔑 3 · Configure the OpenAI API Key

Copy the example file and add your key:
```
cp .env.example .env
```

Edit .env:
```
OPENAI_API_KEY=sk-your-key-here
```


