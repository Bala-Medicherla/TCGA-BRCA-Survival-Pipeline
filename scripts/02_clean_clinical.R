# ------------------------------------------------------------
# Script: 02_clean_clinical.R
# Purpose:
#   Prepare an analysis-ready clinical dataset for survival analysis.
#   This includes:
#     1) Defining Overall Survival (OS) time and censoring
#     2) Cleaning a simple stage grouping (I/II/III/IV)
#     3) Standardizing age (years)
#     4) Saving:
#        - analysis dataset
#        - missingness summary
#        - exclusions summary
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_file <- file("results/logs/cleaning_log.txt", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("Starting clinical cleaning + survival endpoint derivation...\n\n")

# ------------------------------------------------------------
# 1) Load raw clinical data downloaded in Script 01
# ------------------------------------------------------------
clinical_raw <- readRDS("data_processed/clinical_raw.rds")
cat("Raw clinical rows:", nrow(clinical_raw), "\n")
cat("Raw clinical cols:", ncol(clinical_raw), "\n\n")

# Helper: find first column name that exists in the dataset
pick_first_existing <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  return(NA_character_)
}

# ------------------------------------------------------------
# 2) Identify survival-related columns (TCGA naming can vary)
# ------------------------------------------------------------
days_to_death_col <- pick_first_existing(clinical_raw, c("days_to_death"))
days_to_lfu_col   <- pick_first_existing(clinical_raw, c("days_to_last_follow_up", "days_to_last_followup"))

if (is.na(days_to_death_col) || is.na(days_to_lfu_col)) {
  stop("Expected survival columns not found. Check names(clinical_raw).")
}

cat("Using survival columns:\n")
cat(" - days_to_death:", days_to_death_col, "\n")
cat(" - days_to_last_follow_up:", days_to_lfu_col, "\n\n")

# Stage columns may vary
stage_col <- pick_first_existing(clinical_raw, c("ajcc_pathologic_stage", "tumor_stage"))
if (!is.na(stage_col)) cat("Using stage column:", stage_col, "\n\n")

# Age column may vary
age_col <- pick_first_existing(clinical_raw, c("age_at_diagnosis", "age"))
if (!is.na(age_col)) cat("Using age column:", age_col, "\n\n")

# ------------------------------------------------------------
# 3) Derive Overall Survival (OS)
# ------------------------------------------------------------
dat <- clinical_raw %>%
  mutate(
    vital_status = if ("vital_status" %in% names(clinical_raw)) as.character(vital_status) else NA_character_,
    os_event = case_when(
      !is.na(vital_status) & str_to_lower(vital_status) == "dead"  ~ 1L,
      !is.na(vital_status) & str_to_lower(vital_status) == "alive" ~ 0L,
      TRUE ~ NA_integer_
    ),
    days_to_death = suppressWarnings(as.numeric(.data[[days_to_death_col]])),
    days_to_last_follow_up = suppressWarnings(as.numeric(.data[[days_to_lfu_col]])),
    os_time_days = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up)
  )

# ------------------------------------------------------------
# 4) Clean stage into a simple group (I/II/III/IV)
# ------------------------------------------------------------
if (!is.na(stage_col)) {
  dat <- dat %>%
    mutate(
      stage_raw = str_squish(as.character(.data[[stage_col]])),
      stage_group = case_when(
        str_detect(str_to_lower(stage_raw), "stage i\\b")   ~ "I",
        str_detect(str_to_lower(stage_raw), "stage ii\\b")  ~ "II",
        str_detect(str_to_lower(stage_raw), "stage iii\\b") ~ "III",
        str_detect(str_to_lower(stage_raw), "stage iv\\b")  ~ "IV",
        TRUE ~ NA_character_
      )
    )
} else {
  dat <- dat %>% mutate(stage_group = NA_character_)
}

# ------------------------------------------------------------
# 5) Standardize age (years)
#   - TCGA sometimes stores age in days, sometimes in years.
#   - We use a simple rule: if > 200, treat as days and convert.
# ------------------------------------------------------------
if (!is.na(age_col)) {
  dat <- dat %>%
    mutate(
      age_raw = suppressWarnings(as.numeric(.data[[age_col]])),
      age_years = ifelse(!is.na(age_raw) & age_raw > 200, age_raw / 365.25, age_raw)
    )
} else {
  dat <- dat %>% mutate(age_years = NA_real_)
}

# ------------------------------------------------------------
# 6) Summarize exclusions + missingness (transparent reporting)
# ------------------------------------------------------------
exclusions <- tibble(
  rule = c("Missing/invalid OS time (<=0 or NA)", "Missing event status (vital_status unknown)"),
  n_excluded = c(sum(is.na(dat$os_time_days) | dat$os_time_days <= 0),
                 sum(is.na(dat$os_event)))
)

write_csv(exclusions, "results/tables/exclusions.csv")

missingness <- tibble(
  field = c("os_time_days", "os_event", "stage_group", "age_years"),
  n_missing = c(sum(is.na(dat$os_time_days)),
                sum(is.na(dat$os_event)),
                sum(is.na(dat$stage_group)),
                sum(is.na(dat$age_years)))
) %>%
  mutate(pct_missing = round(100 * n_missing / nrow(dat), 2))

write_csv(missingness, "results/tables/missingness.csv")

cat("Saved exclusions summary to results/tables/exclusions.csv\n")
cat("Saved missingness summary to results/tables/missingness.csv\n\n")

# ------------------------------------------------------------
# 7) Create final analysis dataset
#   Minimum required fields: OS time + event
# ------------------------------------------------------------
dat_clean <- dat %>%
  filter(!is.na(os_time_days), os_time_days > 0, !is.na(os_event)) %>%
  mutate(
    os_event = as.integer(os_event),
    stage_group = factor(stage_group, levels = c("I", "II", "III", "IV"))
  )

cat("Final analysis dataset rows:", nrow(dat_clean), "\n\n")

# Save outputs
saveRDS(dat_clean, "data_processed/clinical_surv.rds")
write_csv(dat_clean, "data_processed/clinical_surv.csv")

cat("Saved analysis dataset:\n")
cat(" - data_processed/clinical_surv.rds\n")
cat(" - data_processed/clinical_surv.csv\n")

sink(type = "output")
sink(type = "message")
close(log_file)