# ------------------------------------------------------------
# Script: 04_model_validation.R
# Purpose:
#   Perform a simple internal validation of the Cox model using
#   1) a train/test split
#   2) bootstrap validation (out-of-bag C-index)
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
on.exit({
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}, add = TRUE)

cat("Starting internal validation (train/test)...\n\n")

set.seed(2026)  # reproducible split and bootstrap

# Load analysis dataset
dat <- readRDS("data_processed/clinical_surv.rds") %>%
  filter(!is.na(age_years), !is.na(stage_group))

cat("Rows available for validation:", nrow(dat), "\n\n")

if (nrow(dat) == 0 || sum(dat$os_event, na.rm = TRUE) == 0) {
  stop("Insufficient data/events for validation.")
}

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

# ------------------------------------------------------------
# 2) Bootstrap validation (out-of-bag C-index)
# ------------------------------------------------------------
cat("Starting bootstrap validation (out-of-bag C-index)...\n\n")

bootstrap_iters <- 200
oob_cindex <- numeric(bootstrap_iters)

for (i in seq_len(bootstrap_iters)) {
  boot_idx <- sample(seq_len(n), replace = TRUE)
  oob_idx <- setdiff(seq_len(n), unique(boot_idx))
  
  if (length(oob_idx) < 5) {
    oob_cindex[i] <- NA_real_
    next
  }
  
  boot_train <- dat[boot_idx, ]
  boot_oob <- dat[oob_idx, ]
  
  boot_fit <- coxph(
    Surv(os_time_days, os_event) ~ age_years + stage_group,
    data = boot_train
  )
  
  lp_oob <- predict(boot_fit, newdata = boot_oob, type = "lp")
  c_obj_oob <- concordance(Surv(os_time_days, os_event) ~ lp_oob, data = boot_oob)
  oob_cindex[i] <- unname(c_obj_oob$concordance)
}

oob_cindex_clean <- oob_cindex[!is.na(oob_cindex)]

boot_summary <- c(
  paste0("Bootstrap iterations: ", bootstrap_iters),
  paste0("OOB C-index mean: ", round(mean(oob_cindex_clean), 4)),
  paste0("OOB C-index SD: ", round(sd(oob_cindex_clean), 4)),
  paste0("OOB C-index N: ", length(oob_cindex_clean))
)

writeLines(boot_summary, "results/tables/bootstrap_c_index.txt")
cat("Saved bootstrap C-index summary to results/tables/bootstrap_c_index.txt\n\n")

cat("Bootstrap validation completed successfully.\n")

