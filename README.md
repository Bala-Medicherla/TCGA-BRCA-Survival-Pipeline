# Real-World Survival Analysis in Breast Cancer Using TCGA Clinical data

## Project overview:
This project presents a real-world survival analysis workflow using publicly available clinical data from The Cancer Genome Atlas(TCGA) breast cancer cohort. The goal is to study overall survival patterns in an observational oncology dataset and to evaluate how well commonly collected clinical variables explain differences in patient outcomes.

The focus of this work is not algorthmic novelty, but rather the careful handling of clinical endpoints , transparent modeling choices, and reproducble analysis. These are important in applied biomedical and health informatics research.

## Problem Statement:
Clinical trial data are generally generated under controlled conditions, but real-world oncology data often reflect patient heterogenity, incomplete follow-up, and varaibility in clinical documentation. As a result, survival analysis using observational clinical datasets requires careful definition of endpoints, explicit handling of censoring, and validation of model assumptions.

This project will address two primary questions:
1. How does overall survival differ across clinically meaningful patient groups in a real world breast cancer cohort?
2. Can a multivariable survival model improve risk stratification beyond simple subgroup comparisons and do the model assumptions hold when applied to real clinical data?

The broader goal is to build a reproducible and interpretable survival analysis pipeline that reflects challenges commonly encountered in real-world biomedical data.

## Dataset:
The analysis uses publicly available clinical data from the **TCGA Breast Invasive Carcinoma(TCGA-BRCA)** cohort.Data are accessed programmatically using the TCGAbiolinks R package.

Only open access, de-idetified clinical variables are used. No protected or controlled access molecular data are required. The primary endpoint for this analysis is **Overall Survival(OS)**

## Methods:
**Endpoint Defnition**:
Overall survival is defined using documented TCGA clinical fields. Survival time is derived from days-to-death or days-to-last-follow-up, and event status is defined as death versus censoring. Records with missing or invalid survival information are excluded using clearly documented rules.Endpoint derivation is treated as a core methodological step rather than an assumption.
   
 **Descriptive Survival Analysis**:
 Kaplan–Meier survival curves are generated to visualize survival patterns across clinically relevant patient groups, such as disease stage.  Differences between groups are assessed using log-rank tests.

 **Multivariable Survival Modeling**:
 A Cox proportional hazards model is fitted using available clinical covariates such as age and stage. Model interpretation focuses on      effect direction, magnitude, and clinical plausibility rather than purely statistical significance.The proportional hazards assumption is explicitly evaluated using diagnostic tests. When violations are identified, alternative model specifications are explored and documented.


 **Model Validation**:
 To avoid over-optimistic conclusions, internal validation is performed using concordance-based metrics. Model performance is assessed using both a train-test split and bootstrap-based out-of-bag validation to provide a more robust estimate of generalization performance

  **Model Validation**:
    To avoid over-optimistic conclusions, internal validation is performed using concordance-based metrics. Model performance is assessed using both a train-test split and bootstrap-based out-of-bag validation to provide a more robust estimate of generalization performance. 

    
 **Missing Data Handling**:
 Patterns of missing data are summarized and reported explicitly. The primary analysis uses a complete-case approach for transparency, with sensitivity summaries provided to assess the potential impact of excluded observations.
    
  **Reproducibility**:
  The analysis is organized as a scripted pipeline with clearly defined steps. To fully pin package versions, initialize `renv` in your environment (`renv::init()`), then snapshot with `renv::snapshot()` so runs can be reproduced later. All figures, tables, and logs are generated programmatically and saved to structured directories.
    The analysis is organized as a scripted pipeline with clearly defined steps. To fully pin package versions, initialize `renv` in your environment (`renv::init()`), then snapshot with `renv::snapshot()` and commit the resulting `renv.lock` so runs can be reproduced later. All figures, tables, and logs are generated programmatically and saved to structured directories.
    
## Key Outputs:
   1. Kaplan–Meier survival plots
   2. Log-rank test results
   3. Cox model summary tables
   4. Proportional hazards diagnostics
   5. Internal validation metrics(train/test + bootstrap OOB C-index)
   6. A rendered analysis report (HTML or PDF)
   7. (Optional) Clinical + molecular Cox model results

## Results & Interpretation
Below are example outputs you can surface directly in the README after running the pipeline.

### Kaplan–Meier curves by stage
![Kaplan–Meier by stage](results/figures/km_by_stage.png)

**Interpretation**
- Later-stage disease shows visibly worse overall survival than early-stage groups, consistent with clinical expectations.
- The log-rank test output (see `results/tables/logrank_stage.txt`) provides a formal test of stage-stratified survival differences.

### Proportional hazards diagnostics
![PH diagnostics](results/figures/ph_diagnostics.png)

**Interpretation**
- The diagnostics indicate whether the proportional hazards assumption holds for key covariates.
- If violations appear, consider stratification or time-varying effects and document the decision.

### Model performance (internal validation)
**Interpretation**
- The train/test C-index (see `results/tables/c_index.txt`) gives a baseline estimate of discrimination.
- The bootstrap out-of-bag C-index (see `results/tables/bootstrap_c_index.txt`) provides a more robust internal validation estimate.

## Optional Molecular Extension (Transcriptome-augmented model):
 To strengthen biological depth, an optional script integrates TCGA-BRCA gene expression data and builds a simple proliferation signature. The signature (average z-score of proliferation genes) is then evaluated alongside age and stage in a Cox model.

  **Script**: `scripts/05_molecular_signature.R`

  **What it does**:
   - Downloads TCGA-BRCA gene expression (STAR - Counts) via TCGAbiolinks
   - Computes a proliferation score from a small gene set (e.g., MKI67, CCNB1, CDK1)
   - Merges expression-derived scores with clinical survival data
   - Fits a Cox model (age + stage + proliferation score)

  **How to run**:
   1. Run the clinical pipeline first: `Rscript run_pipeline.R`
   2. Run the molecular extension: `Rscript scripts/05_molecular_signature.R`

  **Note**: Expression downloads can be large and may take time depending on bandwidth and GDC availability.


## Limitations:
This analysis is based on observational clinical data and is subject to missingness and potential confounding.Treatment exposure details are not modeled in this clinical-only version, and results represent associations rather than causal effects.

These limitations are inherent to real-world clinical datasets and motivate future extensions of this work.

## Next Steps:
Potential extensions include external validation using independent public cohorts, more robust missing data strategies, and controlled integration of additional clinical or molecular features.

## Why this Project matters:
This project reflects the type of applied, methodologically careful analysis commonly required in biomedical informatics and clinical data science. Emphasis is placed on endpoint definition, assumption checking, validation, and reproducibility all of which are central to real-world clinical research workflows.
