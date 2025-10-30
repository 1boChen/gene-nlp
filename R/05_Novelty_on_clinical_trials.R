###############################################
# Novelty scoring (year-based, phase-aware)
###############################################

required_packages <- c("dplyr", "readr", "stringr", "lubridate", "purrr")
library(dplyr)
library(readr)
library(stringr)
library(lubridate)
library(purrr)

# ---- Input ----
input_path <- "reports/clinical_trials.tsv"
output_path <- "reports/clinical_trials_novelty.tsv"

trials <- read_tsv(input_path, show_col_types = FALSE)

# ---- Helper: split multi-values ----
split_multi <- function(x) {
  if (is.na(x) || str_trim(x) == "") return(character(0))
  parts <- unlist(strsplit(x, " \\|\\| ", fixed = FALSE))
  parts <- str_trim(parts)
  if (length(parts) == 0) return("")
  if (all(parts == "")) return(rep("", max(1, length(parts))))
  parts
}

# ---- Helper: robust date parser ----
parse_dates <- function(x) {
  as_date(
    suppressWarnings(
      parse_date_time(x, orders = c("Ymd", "Ym", "Y"))
    )
  )
}

# ---- Helper: phase ranking ----
phase_rank <- function(phase_str) {
  if (is.na(phase_str) || str_trim(phase_str) == "") return(0)
  s <- toupper(phase_str)
  if (str_detect(s, "PHASE\\s*IV|PHASE\\s*4|PHASE4")) return(4)
  if (str_detect(s, "PHASE\\s*III|PHASE\\s*3|PHASE3")) return(3)
  if (str_detect(s, "PHASE\\s*II|PHASE\\s*2|PHASE2")) return(2)
  if (str_detect(s, "EARLY\\s*PHASE\\s*1|EARLY\\s*PHASE\\s*I|PHASE\\s*I\\b|PHASE\\s*1\\b|PHASE1")) return(1)
  return(0)
}

# ---- Main scoring ----
get_novelty_score <- function(num_trials, start_dates, last_update_dates, phases) {
  if (is.na(num_trials) || num_trials == 0) return(0)

  start_v  <- split_multi(start_dates)
  update_v <- split_multi(last_update_dates)
  phase_v  <- split_multi(phases)

  n <- max(length(start_v), length(update_v), length(phase_v), 1)
  length(start_v)  <- n
  length(update_v) <- n
  length(phase_v)  <- n
  start_v[is.na(start_v)]  <- ""
  update_v[is.na(update_v)] <- ""
  phase_v[is.na(phase_v)]  <- ""

  start_dt  <- parse_dates(start_v)
  update_dt <- parse_dates(update_v)
  recency_dt <- coalesce(update_dt, start_dt)
  recency_year <- year(recency_dt)
  valid_idx <- which(!is.na(recency_year))
  if (length(valid_idx) == 0) return(1)

  # Recency thresholds
  any_after_2024 <- any(recency_year >= 2024)
  any_after_2023 <- any(recency_year >= 2023)

  if (!any_after_2023) return(1)
  if (!any_after_2024) return(2)

  # Trials in past two years
  recent_idx <- which(recency_year >= 2024)
  recent_phases <- phase_v[recent_idx]

  # ---- Phase-based logic ----
  if (length(recent_phases) == 0) {
    # Default behavior: no phase info but recent trial -> moderate novelty
    return(3)
  }

  phase_vals <- map_dbl(recent_phases, phase_rank)
  max_phase_recent <- ifelse(length(phase_vals) == 0, 0, max(phase_vals, na.rm = TRUE))

  if (max_phase_recent >= 4) return(5)
  if (max_phase_recent == 3) return(4)
  if (max_phase_recent %in% c(1, 2)) return(3)

  return(2)
}

# ---- Apply ----
if (!"start_dates" %in% names(trials)) {
  trials <- trials %>%
    mutate(
      novelty_score = pmap_dbl(
        list(num_trials, "", last_update_dates, phases),
        get_novelty_score
      )
    )
} else {
  trials <- trials %>%
    mutate(
      novelty_score = pmap_dbl(
        list(num_trials, start_dates, last_update_dates, phases),
        get_novelty_score
      )
    )
}

# ---- Save ----
write_tsv(trials, output_path)
cat(sprintf("✅ Novelty scoring complete. Saved: %s\n", normalizePath(output_path)))
