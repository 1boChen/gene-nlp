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

## 🗂️ Repository Layout

gene-nlp
├── data # input data (CSV)
├── external # large reference datasets (LFS)
├── reports # analysis outputs
├── R # numbered R scripts (01–08)
├── src # Python scripts (LLM module)
├── env/environment.yml # Conda environment definition
├── scripts/run.sh # unified pipeline launcher
├── Dockerfile # container recipe
├── .env.example # OpenAI key template
├── .gitattributes # Git LFS tracking rules
└── README.md # this file
