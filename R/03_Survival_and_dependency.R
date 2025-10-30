renv::init()  # optional
pkgs <- c("survival","glmnet","survAUC","timeROC","rsample","yardstick","data.table","tidyverse")
library(survival); library(glmnet); library(timeROC); library(rsample); library(data.table); library(tidyverse)

druggable_candidates <- fread("reports/dgidb_interactions.tsv")     # gene, coef
genes <- unique(druggable_candidates$gene) %>% sort()

expr_file <- "external/expression_surv.csv"     # genes x samples; counts or CPM/TPM
pheno_file <- "external/phenotype_surv.csv" 
if (file.exists(expr_file)) {
  E <- fread(expr_file) %>% na.omit() %>% column_to_rownames('gene') %>% na.omit()
  E <- E[genes,] %>% na.omit()
  pheno <- fread(pheno_file) %>%
    arrange(match(sample_id, colnames(E)))
  
  x <- t(E)  # samples x genes
  y <- Surv(time=pheno$os_time, event=pheno$os_event)
  
  set.seed(7)
  #cvfit <- cv.glmnet(x, y, family="cox", alpha=0.5, nfolds=10, standardize=TRUE)
  #fit <- glmnet(x, y, family="cox", alpha=0.5, lambda=cvfit$lambda.min)
  path <- glmnet(x, y, family="cox", alpha=0.5)
  hits <- as.matrix(coef(path)) %>% as.data.frame()
  #coefs <- as.matrix(coef(fit)); hits <- coefs[coefs!=0,,drop=FALSE]
  cox_hits <- data.table(gene=rownames(hits), coef=as.numeric(hits$s99)) %>%
    filter(!coef==0)
  
  fwrite(cox_hits, "reports/cox_gene_hits.tsv", sep="\t")
}


# --- DepMap gene effect (optional local file) ---
# Expect a CSV with rows=genes, cols=cell lines; compute per-gene mean
dep_file <- "external/CRISPRGeneEffect.csv"
if (file.exists(dep_file)) {
  dep <- fread(dep_file) %>%
    column_to_rownames('V1') %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('gene') %>%
    mutate(gene=gsub(' (.*)', '', gene)) %>%
    as.data.table() 
  # assume first column are gene symbols; adjust if your header differs
  dep_long <- melt(dep, id.vars = "gene", variable.name = "cell_line", value.name = "gene_effect")
  dep_sum <- dep_long[gene %in% genes, .(mean_gene_effect = mean(gene_effect, na.rm=TRUE)), by = gene]
  fwrite(dep_sum, "reports/depmap_gene_effect_summary.tsv", sep="\t")
}
