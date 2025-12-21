# Survival analysis
  A reproducible R project exploring survival analysis and machine learning for time-to-event clinical outcomes using simulated data.

## Project overview:
  1. This project explores how time to event clinical outcomes can be analyzed using both traditional survival analysis techniques and machine learning–based methods.
  2. I developed this project as part of my preparation for PhD applications in Biomedical / Health Informatics, with the goal of demonstrating practical skills in clinical data modeling, statistical reasoning, and reproducible research.
  3. The analysis focuses on a simulated clinical dataset and compares classical survival models with more flexible machine learning approaches to understand treatment effects and patient risk.

## Motive:
 Many important clinical outcomes such as survival, disease progression, or time to relapse are censored which implies that the event of interest is not observed for all patients. Standard regression methods are not appropriate for this type of data. With this project I would like to emphasize the topics on:
   1. Practice modeling time-to-event outcomes in a realistic clinical setting.
   2. Compare classical statistical approaches with machine learning methods.
   3. Understand the strengths and limitations of different survival models.
   4. Build a clean, reproducible workflow that reflects real research practice.

## Dataset:
  The dataset used in this project is fully simulated and doesn't have any real data. It was generated to resemble an oncology style cohort.
  The data include:
   1. Demographic variables (age, sex).
   2. Clinical characteristics (disease stage).
   3. Treatment assignment (control vs treatment).
   4. Continuous biomarker measurements
   5. Follow-up time and event indicator (with censoring)

## Analysis Insights:
   **Survival analysis**:
   1. Kaplan–Meier curves to visualize survival differences between treatment groups.
   2. Log-rank test to assess group differences.
   3. Cox proportional hazards modeling to estimate treatment and covariate effects.
   4. Evaluation of proportional hazards assumptions.
   
  **Machine learning**:
   1. Random Survival Forests to capture non-linear effects and interactions.
   2. Variable importance analysis to identify influential predictors.
   3. Penalized Cox regression (Lasso) for feature selection and model stability

  **Model comparison**:
   1. Comparison of models using the concordance index (C-index).
   2. Discussion of interpretability versus predictive flexibility.

## Tools and methods:
  This project was implemented in R, using commonly used packages in clinical and biomedical informatics research:
   1. Survival, survminer.
   2. randomForestSRC.
   3. glmnet.
   4. ggplot2.
      
## What this project demonstrates:
  This project reflects my interest in Biomedical and Health Informatics, particularly:
   1. Working with clinically meaningful, time-to-event outcomes.
   2. Applying statistical and machine learning methods to healthcare data.
   3. Building reproducible and readable analysis pipelines.
   4. Balancing model interpretability with predictive performance.
   Although the data are simulated, the modeling approaches and analytical reasoning closely mirror which are used in real clinical and translational research.
