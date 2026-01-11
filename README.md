# Real World Survival Analysis in Breast Cancer Using TCGA Clinical data

## Project overview:
 This project presents a real world survival analysis workflow using publicly available clinical data from The Cancer Genome Atlas(TCGA) breast cancer cohort. The goal is to study overall survival patterns in a observational oncology dataset and to evaluate how well commonly collected clinical variables explain differences in patient outcomes.

   The focus of this work is not algortmic novelty, but rather the careful handling of clinical endpoints , transparent modeling choices and reproducble analysis, these are important in applied biomedical and health informatics research.

## Problem Statement:
 Clinical trial data are generally generated under controlled conditions, but real world oncology data often reflect patient heterogenity, incomplete followup and varaibility in clinical documentation. As a result, survival analysis using observational clinical datasets requires careful definition of endpoints, explicit handling of censoring and validation of model assumptions.

  This project will adrress two primary questions:
   1. How does overall survival differ across clinically meaningful patient groups in a real world breast cancer cohort?
   2. Can a multivariable survival model improve risk stratification beyond simple subgroup comparisons and do the model assumptions hold when applied to real clinical data?

The border goal is to build a reproducible and interpretable survival analysis pipeline that reflects challenges commonly encountered in real world biomedical data.

## Dataset:
  The analysis uses publicly available clinical data from the **TCGA Breast Invasive Carcinoma(TCGA-BRCA)** cohort.      Data are accessed programmatically using the TCGAbiolinks R package.

  Only open access, de-idetified clinical variables are used. No protected or controlled access molecular data are      required. The primary endpoint for this analysis is **Overall Survival(OS)**

## Methods:
   **Endpoint Defnition**:
     Overall survival is defined using documented TCGA clinical fields. Survival time is derived from days-to-death or      days-to-last-follow-up, and event status is defined as death versus censoring. Records with missing or invalid         survival information are excluded using clearly documented rules.Endpoint derivation is treated as a core              methodological step rather than an assumption.
   
  **Descriptive Survival Analysis**:
     Kaplan–Meier survival curves are generated to visualize survival patterns across clinically relevant patient           groups, such as disease stage. Differences between groups are assessed using log-rank tests.

  **Multivariable Survival Modeling**:
     A Cox proportional hazards model is fitted using available clinical covariates such as age and stage. Model            interpretation focuses on effect direction, magnitude, and clinical plausibility rather than purely statistical        significance.The proportional hazards assumption is explicitly evaluated using diagnostic tests. When violations       are identified, alternative model specifications are explored and documented.

  **Model Validation**:
    To avoid over-optimistic conclusions, internal validation is performed using concordance-based metrics. Model          performance is assessed using either bootstrap-based validation or a train-test evaluation strategy, depending on      the analysis configuration.
    
  **Missing Data Handling**:
    Patterns of missing data are summarized and reported explicitly. The primary analysis uses a complete-case             approach for transparency, with sensitivity summaries provided to assess the potential impact of excluded              observations.
    
  **Reproducibility**:
    The analysis is organized as a scripted pipeline with clearly defined steps. Package versions are managed using       renv, and all figures, tables, and logs are generated programmatically and saved to structured directories to          ensure reproducibility.
    
## Key Outputs:
   1. Kaplan–Meier survival plots
   2. Log-rank test results
   3. Cox model summary tables
   4. Proportinal hazards diagnostics
   5. Internal validation metrics
   6. A rendered analysis report(HTML or PDF)

## Limitations:
 This analysis is based on observational clinical data and is subject to missingness and potential                      confounding.Treatment exposure details are not modeled in this clinical-only version, and results represent            associations rather than causal effects.

 These limitations are inherent to real-world clinical datasets and motivate future extensions of this work.

## Next Steps:
 Potential extensions include external validation using independent public cohorts, more robust missing data            strategies, and controlled integration of additional clinical or molecular features.

## Why this Project matters:
 This project reflects the type of applied, methodologically careful analysis commonly required in biomedical           informatics and clinical data science. Emphasis is placed on endpoint definition, assumption checking, validation,     and reproducibility all of which are central to real-world clinical research workflows.
