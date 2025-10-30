library(httr)
library(jsonlite)
library(readr)
library(dplyr)
library(stringr)

# ---- Helper: safely get nested fields ----
get_in <- function(x, path) {
  tryCatch({
    for (p in path) x <- x[[p]]
    if (is.null(x)) NA_character_ else as.character(x)
  }, error = function(e) NA_character_)
}

# ---- Core function: query clinical trials for one gene ----
get_clinical_trials_for_gene <- function(gene_name, max_titles = 50, sleep_sec = 1) {
  base_url <- "https://clinicaltrials.gov/api/v2/studies"
  params <- list(
    "query.term" = gene_name,
    "pageSize"   = max_titles,
    "countTotal" = "true"
  )

  resp <- GET(base_url, query = params)
  if (http_error(resp)) {
    warning(sprintf("Trials fetch failed for %s (HTTP %s)", gene_name, status_code(resp)))
    return(list(num_trials = NA_integer_, titles = "", start_dates = "", last_update_dates = "", phases = ""))
  }

  dat <- fromJSON(content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  total <- if (!is.null(dat$totalCount)) as.integer(dat$totalCount) else 0
  studies <- dat$studies
  if (is.null(studies) || length(studies) == 0) {
    return(list(num_trials = total, titles = "", start_dates = "", last_update_dates = "", phases = ""))
  }
  if (!is.list(studies[[1]])) studies <- list(studies)

  titles <- character()
  start_dates <- character()
  last_updates <- character()
  phases <- character()

  for (st in studies) {
    title <- get_in(st, c("protocolSection", "identificationModule", "briefTitle"))
    start_date <- get_in(st, c("protocolSection", "statusModule", "startDateStruct", "date"))
    last_update <- get_in(st, c("protocolSection", "statusModule", "lastUpdatePostDateStruct", "date"))

    # ✅ Extract phase info correctly
    phase_val <- tryCatch({
      phase_data <- st$protocolSection$designModule$phases
      if (is.null(phase_data)) {
        NA_character_
      } else {
        phase_vec <- unlist(phase_data, use.names = FALSE)
        paste(unique(phase_vec), collapse = ", ")
      }
    }, error = function(e) NA_character_)

    titles <- c(titles, ifelse(is.na(title), "", title))
    start_dates <- c(start_dates, ifelse(is.na(start_date), "", start_date))
    last_updates <- c(last_updates, ifelse(is.na(last_update), "", last_update))
    phases <- c(phases, ifelse(is.na(phase_val), "", phase_val))
  }

  Sys.sleep(sleep_sec)  # be gentle to API
  list(
    num_trials = total,
    titles = paste(titles, collapse = " || "),
    start_dates = paste(start_dates, collapse = " || "),
    last_update_dates = paste(last_updates, collapse = " || "),
    phases = paste(phases, collapse = " || ")
  )
}

# ---- Main: read input genes and write results ----
input_path  <- "reports/dgidb_interactions.tsv"
output_path <- "reports/clinical_trials.tsv"

dg <- read_tsv(input_path, show_col_types = FALSE)
stopifnot("gene" %in% names(dg))
genes <- unique(dg$gene)
cat(sprintf("Found %d unique genes.\n", length(genes)))

rows <- vector("list", length(genes))
for (i in seq_along(genes)) {
  g <- genes[i]
  cat(sprintf("[%d/%d] %s\n", i, length(genes), g))
  info <- get_clinical_trials_for_gene(g, max_titles = 50, sleep_sec = 1)
  rows[[i]] <- data.frame(
    gene = g,
    num_trials = info$num_trials,
    titles = info$titles,
    start_dates = info$start_dates,
    last_update_dates = info$last_update_dates,
    phases = info$phases,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

# ---- Combine and clean output ----
out <- bind_rows(rows)

# Normalize phase labels (PHASE1 → Phase 1)
out <- out %>%
  mutate(phases = str_replace_all(phases, "PHASE([0-9])", "Phase \\1"))

# ---- Save results ----
write_tsv(out, output_path)
cat(sprintf("✅ Saved: %s\n", normalizePath(output_path)))
