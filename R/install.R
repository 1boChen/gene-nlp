# Ensure CRAN and Bioconductor repos
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  BioC_mirror = "https://bioconductor.org"
)

# Explicitly define target library (inside micromamba env)
target_lib <- "/opt/conda/envs/gene_nlp/lib/R/library"
dir.create(target_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(target_lib, .libPaths()))

# Install BiocManager if missing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = target_lib)

# Initialize BiocManager and set Bioconductor version
BiocManager::install(version = "3.22", ask = FALSE, update = TRUE, lib = target_lib)

# CRAN + Bioconductor packages
cran_pkgs <- c(
  "data.table", "tidyverse", "httr", "jsonlite", "xml2", "ghql",
  "glue", "readr", "dplyr", "stringr", "lubridate", "purrr",
  "survival", "glmnet", "timeROC", "rsample", "yardstick", "scales"
)

bio_pkgs <- c("AnnotationDbi", "org.Hs.eg.db")

# Install CRAN packages
for (p in cran_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org", lib = target_lib)

# Install Bioconductor packages
for (p in bio_pkgs)
  if (!requireNamespace(p, quietly = TRUE))
    BiocManager::install(p, ask = FALSE, update = TRUE, lib = target_lib)
