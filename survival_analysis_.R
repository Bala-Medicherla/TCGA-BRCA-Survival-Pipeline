############################################################
# Survival Analysis + Machine Learning for Time-to-Event Data
############################################################

rm(list=ls())
############################
# SECTION 0: USER SETTINGS 
############################
project_dir <- "/Users/srirammedicherla/Documents/Survival_analysis"

if (nzchar(project_dir)) {
  if (!dir.exists(project_dir)) dir.create(project_dir, recursive = TRUE)
  setwd(project_dir)
}

# Output folder for figures
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir)


############################
# Package setup
############################

# We install missing packages and load everything we need.
packages <- c(
  "survival",
  "survminer",
  "randomForestSRC",
  "glmnet",
  "ggplot2"
)

to_install <- packages[!packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(survival)
library(survminer)
library(randomForestSRC)
library(glmnet)
library(ggplot2)


############################################################
# SECTION 1: Simulate a clinical-style survival dataset
############################################################

# Simulated variables:
# - age, sex, stage, biomarkers, treatment
# Survival outcome:
# - time: observed follow-up time
# - event: 1 if event occurred,
#          0 if censored
#
# Important note:
# - This is a “Cox-style” simulation: covariates influence hazard via a linear predictor.
set.seed(123)
n <- 500

dat <- data.frame(
  id         = 1:n,
  age        = round(rnorm(n, mean = 60, sd = 10)),
  sex        = rbinom(n, 1, 0.45),          # 0 = Female, 1 = Male
  treatment  = rbinom(n, 1, 0.50),          # 0 = Control, 1 = Treatment
  biomarker1 = rnorm(n),
  biomarker2 = rnorm(n),
  stage      = sample(1:4, n, replace = TRUE)
)

# Risk score (linear predictor): higher score -> higher hazard -> shorter survival
lp <- 0.04 * (dat$age - 60) +
  0.5  * dat$sex +
  -0.7 * dat$treatment +
  0.6  * dat$biomarker1 +
  0.3  * (dat$stage - 2)

# Generate event times from exponential baseline hazard with covariate effects
baseline_hazard <- 0.05
u <- runif(n)
true_time <- -log(u) / (baseline_hazard * exp(lp))

# Independent censoring (loss to follow-up / study end)
censor_time <- rexp(n, rate = 0.03)

# Observed time + event indicator
dat$time  <- pmin(true_time, censor_time)
dat$event <- as.integer(true_time <= censor_time)

# Create survival object
surv_obj <- Surv(time = dat$time, event = dat$event)

# Quick sanity checks
message("Simulated dataset created (n = ", n, "). Event rate = ", round(mean(dat$event), 3))

############################################################
# SECTION 2: Kaplan–Meier curves (treatment comparison)
############################################################
fit_km <- survfit(surv_obj ~ treatment, data = dat)

km_plot <- ggsurvplot(
  fit_km,
  data         = dat,
  pval         = TRUE,
  risk.table   = TRUE,
  legend.labs  = c("Control", "Treatment"),
  legend.title = "Group",
  xlab         = "Time",
  ylab         = "Survival probability",
  ggtheme      = theme_minimal()
)

ggsave(
  filename = file.path(fig_dir, "km_curve.png"),
  plot     = km_plot$plot,
  width    = 8,
  height   = 5,
  dpi      = 300
)

############################################################
# SECTION 3: Cox proportional hazards model + PH check
############################################################
cox2 <- coxph(
  surv_obj ~ treatment + age + sex + stage + biomarker1 + biomarker2,
  data = dat
)

# Forest plot (interpretable summary)
cox_forest <- ggforest(cox2, data = dat)

ggsave(
  filename = file.path(fig_dir, "cox_forest.png"),
  plot     = cox_forest,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# Proportional hazards assumption test
ph_test <- cox.zph(cox2)

png(file.path(fig_dir, "cox_ph_test.png"), width = 1200, height = 900, res = 150)
par(mfrow = c(2, 3))
plot(ph_test)
dev.off()
par(mfrow = c(1, 1))

############################################################
# SECTION 4: Random Survival Forest (machine learning)
############################################################
dat_ml <- dat[, c("time", "event", "age", "sex", "treatment",
                  "biomarker1", "biomarker2", "stage")]

set.seed(123)
rsf_fit <- rfsrc(
  Surv(time, event) ~ .,
  data       = dat_ml,
  ntree      = 1000,
  importance = TRUE
)

vip <- data.frame(
  variable   = names(rsf_fit$importance),
  importance = as.numeric(rsf_fit$importance)
)

vip_plot <- ggplot(vip, aes(x = reorder(variable, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  xlab("Variable") +
  ylab("Variable importance") +
  theme_minimal()

ggsave(
  filename = file.path(fig_dir, "rsf_variable_importance.png"),
  plot     = vip_plot,
  width    = 8,
  height   = 5,
  dpi      = 300
)

############################################################
# SECTION 5: Penalized Cox regression (Lasso)
############################################################
# Why Lasso:
# - Helps with feature selection and stability, especially as variables grow.
x <- model.matrix(
  ~ age + sex + treatment + biomarker1 + biomarker2 + factor(stage),
  data = dat
)[, -1]

y <- Surv(dat$time, dat$event)

set.seed(123)
fit_lasso <- cv.glmnet(
  x, y,
  family = "cox",
  alpha  = 1
)

png(file.path(fig_dir, "lasso_cv_plot.png"), width = 1200, height = 900, res = 150)
plot(fit_lasso)
dev.off()


#######
# END
#######
