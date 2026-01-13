# ------------------------------------------------------------
# Script: 05_molecular_signature.R
# Purpose:
#   Optional molecular extension for TCGA-BRCA:
#     1) Download gene expression (HTSeq FPKM)
#     2) Compute a simple proliferation signature
#     3) Merge with clinical survival data
#     4) Fit a Cox model (age + stage + proliferation score)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(survival)
})

dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_file <- file("results/logs/molecular_log.txt", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}, add = TRUE)

cat("Starting optional molecular extension (TCGA-BRCA gene expression)...\n")
cat("This step may take a while depending on download speed.\n\n")

# ------------------------------------------------------------
# 1) Load clinical survival dataset
# ------------------------------------------------------------
clinical <- readRDS("data_processed/clinical_surv.rds")
cat("Loaded clinical survival dataset.\n")
cat("Clinical dataset rows:", nrow(clinical), "\n\n")

pick_first_existing <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  return(NA_character_)
}

patient_id_col <- pick_first_existing(
  clinical,
  c("bcr_patient_barcode", "submitter_id", "case_id", "patient_id")
)

if (is.na(patient_id_col)) {
  stop("Patient ID column not found in clinical dataset.")
}

clinical <- clinical %>%
  mutate(patient_id = str_sub(as.character(.data[[patient_id_col]]), 1, 12))

# ------------------------------------------------------------
# 2) Download gene expression (HTSeq FPKM)
# ------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - FPKM"
)

GDCdownload(query)
se <- GDCprepare(query)
expr <- assay(se)

cat("Expression matrix loaded.\n")
cat("Dimensions:", nrow(expr), "genes x", ncol(expr), "samples\n\n")

# ------------------------------------------------------------
# 3) Build a simple proliferation signature
# ------------------------------------------------------------
prolif_genes <- c("MKI67", "CCNB1", "CDK1", "TOP2A", "BUB1")
gene_map <- rowData(se) %>% as.data.frame()

if (!"gene_name" %in% names(gene_map)) {
  stop("gene_name not found in expression rowData. Check GDCprepare output.")
}

prolif_rows <- gene_map %>%
  filter(gene_name %in% prolif_genes) %>%
  pull(row.names)

if (length(prolif_rows) == 0) {
  stop("No proliferation genes found in expression matrix.")
}

expr_prolif <- expr[prolif_rows, , drop = FALSE]
expr_z <- t(scale(t(expr_prolif)))
prolif_score <- colMeans(expr_z, na.rm = TRUE)

expr_scores <- tibble(
  sample_barcode = colnames(expr),
  patient_id = str_sub(sample_barcode, 1, 12),
  proliferation_score = as.numeric(prolif_score)
) %>%
  group_by(patient_id) %>%
  summarize(proliferation_score = mean(proliferation_score, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 4) Merge with clinical data and fit Cox model
# ------------------------------------------------------------
merged <- clinical %>%
  inner_join(expr_scores, by = "patient_id") %>%
  filter(!is.na(stage_group), !is.na(age_years), !is.na(proliferation_score))

cat("Merged clinical + expression data.\n")
cat("Merged dataset rows:", nrow(merged), "\n\n")

if (nrow(merged) == 0 || sum(merged$os_event, na.rm = TRUE) == 0) {
  stop("Insufficient data/events after merging molecular scores.")
}

cox_fit <- coxph(
  Surv(os_time_days, os_event) ~ age_years + stage_group + proliferation_score,
  data = merged
)

cox_sum <- summary(cox_fit)
cox_table <- as.data.frame(cox_sum$coefficients) %>%
  rownames_to_column("term")

write_csv(cox_table, "results/tables/cox_molecular_summary.csv")
cat("Saved molecular Cox summary to results/tables/cox_molecular_summary.csv\n")
cat("Tip: compare this model to the clinical-only Cox model for context.\n")

cat("Molecular extension completed successfully.\n")