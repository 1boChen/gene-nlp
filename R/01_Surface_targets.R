.libPaths(c("/opt/conda/envs/gene_nlp/lib/R/library", .libPaths()))
library(data.table); library(tidyverse); library(org.Hs.eg.db)

dir.create("reports", showWarnings = FALSE, recursive = TRUE)
expr_file <- "data/expression.csv"     # genes x samples; counts or CPM/TPM
pheno_file <- "data/phenotype.csv" 
E <- fread(expr_file) %>% na.omit() %>% column_to_rownames('gene') %>% na.omit()
pheno <- fread(pheno_file) %>%
  arrange(match(sample_id, colnames(E)))
x <- t(E); y <- pheno$response # 0/1

de_results <- data.frame(
  gene = rownames(E),
  p_value = sapply(rownames(E), function(gene) {
    t.test(E[gene, y == 0], E[gene, y == 1])$p.value
  }),
  fold_change = sapply(rownames(E), function(gene) {
    mean(E[gene, y == 1]) / mean(E[gene, y == 0])
  })
) %>%
  mutate(
    adj_p_value = p.adjust(p_value, method = "BH"),  # Benjamini-Hochberg correction
    significant = adj_p_value < 0.05
  ) %>%
  arrange(adj_p_value)

surface_genes <- unique(c(
  # Cell surface genes
  AnnotationDbi::select(org.Hs.eg.db,
                        keys = "GO:0009986",
                        columns = "SYMBOL",
                        keytype = "GO")$SYMBOL,
  # Plasma membrane genes  
  AnnotationDbi::select(org.Hs.eg.db,
                        keys = "GO:0005886",
                        columns = "SYMBOL",
                        keytype = "GO")$SYMBOL
))

# Filter for significant surface proteins
surface_candidates <- de_results %>%
  filter(p_value < 0.05 & gene %in% surface_genes) %>%
  arrange(adj_p_value)

fwrite(surface_candidates, "reports/surface_candidates.tsv", sep="\t")
detach("package:org.Hs.eg.db", unload = TRUE)

