###########################################################
# PubMed abstracts by gene (up to 50 per gene), with PMIDs
# - Records total_found and num_retrieved
# - Adds pmid1..pmid50 alongside title1..title50, abstract1..abstract50
# - Writes to reports/pubmed_abstracts.tsv
###########################################################

# 0) Install & load packages
required_packages <- c("httr", "jsonlite", "xml2", "readr", "dplyr", "stringr")
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(dplyr)
library(stringr)

# 1) Config
input_path <- "reports/dgidb_interactions.tsv"
output_dir <- dirname(input_path)
output_path <- file.path(output_dir, "pubmed_abstracts.tsv")

NCBI_TOOL  <- "gene_pubmed_fetcher"
NCBI_EMAIL <- ""  # optional

MAX_PER_GENE <- 50
REQUEST_DELAY_SEC <- 0.4

# 2) Helpers ---------------------------------------------------------------

clean_text <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  x <- paste(x, collapse = " ")
  x <- str_replace_all(x, "[\\r\\n\\t]+", " ")
  str_squish(x)
}

build_query <- function(gene) {
  gene_escaped <- gsub('"', '\\"', gene)
  sprintf(
    '((cancer[Title/Abstract] OR leukemia[Title/Abstract] OR lymphoma[Title/Abstract] OR tumor[Title/Abstract] OR carcinoma[Title/Abstract]) AND (therapy[Title/Abstract] OR treatment[Title/Abstract] OR drug[Title/Abstract] OR inhibitor[Title/Abstract] OR "clinical trial"[Publication Type] OR "ClinicalTrials.gov"[Title/Abstract])) AND ("%s"[Title/Abstract])',
    gene_escaped
  )
}

pubmed_esearch <- function(query, retmax = MAX_PER_GENE) {
  esearch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  params <- list(
    db = "pubmed",
    term = query,
    retmax = retmax,
    retmode = "json",
    tool = NCBI_TOOL
  )
  if (nzchar(NCBI_EMAIL)) params$email <- NCBI_EMAIL

  resp <- GET(esearch_url, query = params)
  if (http_error(resp)) {
    warning("esearch failed: HTTP ", status_code(resp))
    return(list(ids = character(0), total_found = 0))
  }
  Sys.sleep(REQUEST_DELAY_SEC)

  js <- fromJSON(content(resp, as = "text", encoding = "UTF-8"))
  ids <- js$esearchresult$idlist
  total <- as.integer(js$esearchresult$count)
  if (is.null(ids)) ids <- character(0)
  list(ids = ids, total_found = total)
}

# Return data.frame(pmid, title, abstract) for the given PMIDs
pubmed_efetch <- function(pmids) {
  if (length(pmids) == 0) return(data.frame(pmid = character(0), title = character(0), abstract = character(0)))

  efetch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  params <- list(
    db = "pubmed",
    id = paste(pmids, collapse = ","),
    retmode = "xml",
    tool = NCBI_TOOL
  )
  if (nzchar(NCBI_EMAIL)) params$email <- NCBI_EMAIL

  resp <- GET(efetch_url, query = params)
  if (http_error(resp)) {
    warning("efetch failed: HTTP ", status_code(resp))
    return(data.frame(pmid = character(0), title = character(0), abstract = character(0)))
  }
  Sys.sleep(REQUEST_DELAY_SEC)

  doc <- read_xml(content(resp, as = "text", encoding = "UTF-8"))
  arts <- xml_find_all(doc, ".//PubmedArticle")
  if (length(arts) == 0) arts <- xml_find_all(doc, ".//PubmedBookArticle")

  pmids_out <- character(length(arts))
  titles <- character(length(arts))
  abstracts <- character(length(arts))

  for (i in seq_along(arts)) {
    art <- arts[[i]]

    # PMID
    pm_node <- xml_find_first(art, ".//MedlineCitation/PMID")
    pmid <- if (!is.na(pm_node) && !is.null(pm_node)) xml_text(pm_node) else ""

    # Title
    title_node <- xml_find_first(art, ".//Article/ArticleTitle")
    if (is.na(title_node) || is.null(title_node)) {
      title_node <- xml_find_first(art, ".//Book/BookTitle")
    }
    title <- if (!is.na(title_node) && !is.null(title_node)) xml_text(title_node) else ""

    # Abstract (concat multiple AbstractText sections, with labels)
    abs_nodes <- xml_find_all(art, ".//Abstract/AbstractText")
    if (length(abs_nodes) > 0) {
      abs_parts <- vapply(abs_nodes, function(n) {
        label <- xml_attr(n, "Label")
        txt <- xml_text(n)
        if (!is.null(label) && nzchar(label)) sprintf("%s: %s", label, txt) else txt
      }, FUN.VALUE = character(1))
      abstract <- paste(abs_parts, collapse = " ")
    } else {
      abstract <- ""
    }

    pmids_out[i] <- clean_text(pmid)
    titles[i] <- clean_text(title)
    abstracts[i] <- clean_text(abstract)
  }

  data.frame(pmid = pmids_out, title = titles, abstract = abstracts, stringsAsFactors = FALSE)
}

# 3) Load genes ------------------------------------------------------------
dgidb <- read_tsv(input_path, show_col_types = FALSE)
if (!"gene" %in% names(dgidb)) stop("The input TSV must contain a 'gene' column.")
genes <- unique(dgidb$gene)
cat(sprintf("Found %d unique genes.\n", length(genes)))

title_cols <- paste0("title", 1:MAX_PER_GENE)
abstract_cols <- paste0("abstract", 1:MAX_PER_GENE)
pmid_cols <- paste0("pmid", 1:MAX_PER_GENE)

# 4) Iterate ---------------------------------------------------------------
rows <- vector("list", length(genes))

for (i in seq_along(genes)) {
  gene <- genes[i]
  cat(sprintf("[%d/%d] %s ... ", i, length(genes), gene))

  query <- build_query(gene)
  es <- pubmed_esearch(query, retmax = MAX_PER_GENE)
  pmids <- es$ids
  total_found <- es$total_found

  df <- pubmed_efetch(pmids)
  num_retrieved <- nrow(df)

  # Cap & pad
  if (num_retrieved > MAX_PER_GENE) df <- df[seq_len(MAX_PER_GENE), , drop = FALSE]

  t_vec <- rep("", MAX_PER_GENE)
  a_vec <- rep("", MAX_PER_GENE)
  p_vec <- rep("", MAX_PER_GENE)
  if (nrow(df) > 0) {
    k <- nrow(df)
    p_vec[1:k] <- df$pmid
    t_vec[1:k] <- df$title
    a_vec[1:k] <- df$abstract
  }

  row_df <- data.frame(
    gene = gene,
    total_found = total_found,
    num_retrieved = min(num_retrieved, MAX_PER_GENE),
    as.list(setNames(as.list(p_vec), pmid_cols)),
    as.list(setNames(as.list(t_vec), title_cols)),
    as.list(setNames(as.list(a_vec), abstract_cols)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  rows[[i]] <- row_df
  cat(sprintf("total: %d, retrieved: %d\n", total_found, min(num_retrieved, MAX_PER_GENE)))
}

final_df <- bind_rows(rows)

# 5) Save TSV --------------------------------------------------------------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write_tsv(final_df, output_path)
cat(sprintf("\n✅ Saved PubMed results to: %s\n", normalizePath(output_path)))
