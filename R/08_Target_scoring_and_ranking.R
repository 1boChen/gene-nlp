library(data.table); library(jsonlite)

cox  <- fread("reports/cox_gene_hits.tsv")
dg   <- if (file.exists("reports/dgidb_interactions.tsv")) fread("reports/dgidb_interactions.tsv") else data.table()
dep  <- if (file.exists("reports/depmap_gene_effect_summary.tsv")) fread("reports/depmap_gene_effect_summary.tsv") else data.table()
lit_clean  <- if (file.exists("reports/llm_lit_novelty.tsv")) fread("reports/llm_lit_novelty.tsv") else data.table()

dt <- merge(cox[, .(gene, coef)], lit_clean[, .(gene, novelty_score)], by="gene", all=TRUE)
if (nrow(dep)) dt <- merge(dt, dep, by="gene", all.x=TRUE)

# normalize components
dt[, cox_score  := tanh(replace(coef, is.na(coef), 0))]
dt[, dep_score  := pmax(0, -replace(mean_gene_effect, is.na(mean_gene_effect), 0))]  # more essential -> higher
# composite priority (tune weights as desired)
dt[, priority := (0.1*cox_score + 0.1*dep_score + 0.5*novelty_score)]
setorder(dt, -priority)
fwrite(dt, "reports/target_ranked_list.tsv", sep="\t")
