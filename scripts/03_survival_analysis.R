# ------------------------------------------------------------
# Script: 03_survival_analysis.R
# Purpose:
#   Generate core survival analysis outputs using the cleaned TCGA-BRCA
#   clinical dataset:
#     1) Kaplan–Meier curves by stage + log-rank test
#     2) Cox proportional hazards model (age + stage)
#     3) Proportional hazards assumption check
#   Outputs are saved under results/figures and results/tables.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
  library(survminer)
})

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs",    showWarnings = FALSE, recursive = TRUE)

log_file <- file("results/logs/analysis_log.txt", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}, add = TRUE)

cat("Starting survival analysis (KM + Cox + diagnostics)...\n\n")

# ------------------------------------------------------------
# 1) Load analysis-ready dataset from Script 02
# ------------------------------------------------------------
dat <- readRDS("data_processed/clinical_surv.rds")
cat("Analysis dataset rows:", nrow(dat), "\n\n")

# For KM by stage, restrict to records with stage_group
km_dat <- dat %>% filter(!is.na(stage_group))
cat("Rows available for KM-by-stage:", nrow(km_dat), "\n\n")

if (nrow(km_dat) == 0) {
  stop("No rows with stage_group available for KM analysis.")
}

surv_obj <- Surv(time = km_dat$os_time_days, event = km_dat$os_event)

# ------------------------------------------------------------
# 2) Kaplan–Meier by stage (with p-value and risk table)
# ------------------------------------------------------------
km_fit <- survfit(surv_obj ~ stage_group, data = km_dat)

km_plot <- ggsurvplot(
  km_fit,
  data = km_dat,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Days",
  ylab = "Overall Survival Probability",
  title = "TCGA-BRCA Overall Survival by Stage",
  legend.title = "Stage"
)

ggsave(
  filename = "results/figures/km_by_stage.png",
  plot = km_plot$plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "results/figures/km_by_stage_risktable.png",
  plot = km_plot$table,
  width = 8,
  height = 3.5,
  dpi = 300
)

cat("Saved KM plots to results/figures/\n\n")

# Log-rank test
lr <- survdiff(surv_obj ~ stage_group, data = km_dat)
writeLines(capture.output(lr), "results/tables/logrank_stage.txt")
cat("Saved log-rank test output to results/tables/logrank_stage.txt\n\n")

# ------------------------------------------------------------
# 3) Cox model (clinical-only): age + stage
# ------------------------------------------------------------
cox_dat <- dat %>% filter(!is.na(age_years), !is.na(stage_group))
cat("Rows available for Cox model:", nrow(cox_dat), "\n\n")

if (nrow(cox_dat) == 0 || sum(cox_dat$os_event, na.rm = TRUE) == 0) {
  stop("Insufficient data/events for Cox model fitting.")
}

cox_fit <- coxph(
  Surv(os_time_days, os_event) ~ age_years + stage_group,
  data = cox_dat
)

cox_sum <- summary(cox_fit)

cox_table <- as.data.frame(cox_sum$coefficients) %>%
  tibble::rownames_to_column("term")

write_csv(cox_table, "results/tables/cox_summary.csv")
cat("Saved Cox summary table to results/tables/cox_summary.csv\n\n")

# ------------------------------------------------------------
# 4) Proportional hazards assumption check
# ------------------------------------------------------------
ph <- cox.zph(cox_fit)
writeLines(capture.output(ph), "results/tables/ph_test.txt")
cat("Saved PH test output to results/tables/ph_test.txt\n\n")

png("results/figures/ph_diagnostics.png", width = 1200, height = 900)
plot(ph)
dev.off()

cat("Saved PH diagnostic plot to results/figures/ph_diagnostics.png\n\n")
cat("Survival analysis script completed successfully.\n")
