pkgs <- c("httr2","jsonlite","data.table","xml2","ghql","tidyverse","glue")
library(httr2); library(jsonlite); library(data.table); library(xml2); library(ghql); library(tidyverse);library(glue)

# read your hits (from steps 2–3)
surface_candidates <- fread("reports/surface_candidates.tsv")     # gene, coef
genes <- unique(surface_candidates$gene) %>% sort()

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
# 1) GraphQL query (parameterized)
gql <- '
query($genes: [String!]!) {
  genes(names: $genes) {
    nodes {
      name
      interactions {
        drug { name conceptId }
        interactionScore
        interactionTypes { type directionality }
        interactionAttributes { name value }
        publications { pmid }
        sources { sourceDbName }
      }
    }
  }
}'

dg_list <- list()
for (i in seq_along(genes)){
  tryCatch(
    {
      req <- request("https://dgidb.org/api/graphql") |>
      req_headers(`Content-Type` = "application/json") |>
      req_body_json(list(query = gql, variables = list(genes = genes[i]))) |>
      req_timeout(60)
    resp <- req_perform(req)
    out  <- resp_body_json(resp, simplifyVector = TRUE)
    dg_list[[i]] <- out$data$genes$nodes %>%
      unnest(cols = where(is.list))  %>%
      unnest(cols = sources) %>%
      unnest(cols = interactionTypes) %>%
      unnest_wider(drug, names_sep = "_drug_") %>%
      dplyr::rename('gene'='name',
             'source'='sourceDbName',
             'drug' = 'drug_drug_name')  %>%
      dplyr::select(gene,drug,interactionScore,type,source)
    }, 
    error = function(e){print(paste0('skipping ', genes[i]))}
  )
}
dg_rows <- rbindlist(dg_list)
fwrite(dg_rows, "reports/dgidb_interactions.tsv", sep="\t")

# --- ChEMBL target search (no key) ---
chembl_list <- list()
for (i in seq_along(genes)){
  tryCatch(
    {
      gene <- genes[i]
      url <- glue("https://www.ebi.ac.uk/chembl/api/data/target/search.json?q={gene}")
      resp <- request(url) |> req_timeout(30) |> req_perform()
      items <- resp_body_json(resp, simplifyVector = TRUE)$targets
      chembl_list[[i]] <- as.data.table(items[, c("target_chembl_id","pref_name")])[
        , `:=`(gene = gene)][, .(gene, chembl_target = target_chembl_id, pref_name)]
    }, 
    error = function(e){print(paste0('skipping ', genes[i]))}
  )
}
chembl_rows <- rbindlist(chembl_list)
fwrite(chembl_rows, "reports/chembl_targets.tsv", sep="\t")

# --- STRING neighbors (9606 human) ---
string_neighbors_one <- function(gene, species=9606) {
  
  url <- glue("https://string-db.org/api/json/network?identifiers={gene}&species=9606")
  resp <- request(url) |> req_timeout(30) |> req_perform()
  
  resp <- request("https://string-db.org/api/json/network") |>
    req_url_query(identifiers = gene, species = species) |>
    req_timeout(30) |>
    req_perform()
  if (resp_status(resp) != 200) return(NULL)
  js <- resp_body_json(resp, simplifyVector = TRUE)
  if (length(js)==0) return(NULL)
  as.data.table(js)[, .(gene = preferredName_A, partner = preferredName_B, score)]
}

str_list <- list()
for (i in seq_along(genes)){
  tryCatch(
    {
      gene <- genes[i]
      species <- 9606
      resp <- request("https://string-db.org/api/json/interaction_partners") |>
        req_url_query(identifiers = gene, species = species) |>
        req_timeout(30) |>
        req_perform()
      js <- resp_body_json(resp, simplifyVector = TRUE, check_type = FALSE)
      str_list[[i]] <- as.data.table(js)[, .(gene = preferredName_A, partner = preferredName_B, score)]
      
      
    }, 
    error = function(e){print(paste0('skipping ', genes[i]))}
  )
}
str_rows <- rbindlist(str_list)
fwrite(str_rows, "reports/string_neighbors.tsv", sep="\t")
