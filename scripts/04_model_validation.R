# ------------------------------------------------------------
# Script: 04_model_validation.R
# Purpose:
#   Perform a simple internal validation of the Cox model using
#   a train/test split and report Harrell's C-index on the test set.
#
# Why this matters:
#   In observational datasets, it is easy to overfit.
#   A basic internal validation step provides a more realistic
#   estimate of model performance.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(readr)
})

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_file <- file("results/logs/validation_log.txt", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("Starting internal validation (train/test)...\n\n")

set.seed(2026)  # reproducible split

# Load analysis dataset
dat <- readRDS("data_processed/clinical_surv.rds") %>%
  filter(!is.na(age_years), !is.na(stage_group))

cat("Rows available for validation:", nrow(dat), "\n\n")

# 70/30 split
n <- nrow(dat)
train_idx <- sample(seq_len(n), size = floor(0.7 * n))

train <- dat[train_idx, ]
test  <- dat[-train_idx, ]

cat("Train size:", nrow(train), "\n")
cat("Test size :", nrow(test), "\n\n")

# Fit Cox model on training set
cox_fit <- coxph(
  Surv(os_time_days, os_event) ~ age_years + stage_group,
  data = train
)

# Linear predictor (risk score) for test set
lp_test <- predict(cox_fit, newdata = test, type = "lp")

# Harrell's C-index on test set using concordance()
c_obj <- concordance(Surv(os_time_days, os_event) ~ lp_test, data = test)
c_index <- unname(c_obj$concordance)

cat("Test-set Harrell's C-index:", round(c_index, 4), "\n\n")

# Save results to a simple text file
out <- c(
  paste0("Train N: ", nrow(train)),
  paste0("Test N: ", nrow(test)),
  paste0("Test-set Harrell's C-index: ", round(c_index, 4))
)

writeLines(out, "results/tables/c_index.txt")
cat("Saved C-index summary to results/tables/c_index.txt\n\n")

cat("Validation completed successfully.\n")

sink(type = "output")
sink(type = "message")
close(log_file)