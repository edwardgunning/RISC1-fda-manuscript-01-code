# ------------------------------------------------------------------------#
# "Modelling"/ Estimation Step
# Fit a series of scalar linear mixed effects models to the FPC scores.
# ------------------------------------------------------------------------#

# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(fda)        # CRAN v5.5.1
library(lme4)       # CRAN v1.1-28


# Read in data: -----------------------------------------------------------
results_path <- here::here("outputs", "results")

basis_transformation_results <- readRDS(
  file.path(results_path, "basis-transform-results.rds"))


# Load helper functions: --------------------------------------------------
source(file = here::here(
  "code", "functions", "BFMM-paper-helper-functions.R"
))


# Obtain data to be used as covariates: -----------------------------------
subject_side_data <- copy(basis_transformation_results$coef_dts[[1]])
subject_side_data[, paste0("lm_coef_", 1:80) := NULL]
head(subject_side_data)



# Extract FPC scores + check: ---------------------------------------------
scores_uncentered <- basis_transformation_results$scores_uncentred
k_retain <- basis_transformation_results$k_retain
stopifnot(ncol(scores_uncentered) == k_retain)


# Construct Data Table with Covariates for Model: -------------------------
covariates_dt <- subject_side_data[, .(ris = factor(retrospective_injury_status, ordered = F),
                                           sex,
                                           side,
                                           speed_cent,
                                           weight_kg_cent,
                                           age_cent,
                                           height_cm_cent,
                                           runner_category,
                                           treadmill, 
                                           dominance, 
                                           bmi_kgm_cent = bmi_kgm_cent,
                                           subject_id = factor(subject_id))]

covariates_dt[, side_is_dominant := ifelse(side == dominance, "dominant", "non-dominant")]
# sense check that each subject should have a n dominant and non dominant side
stopifnot(covariates_dt[, .(count_var = uniqueN(side_is_dominant)), by = subject_id][, all(count_var==2)])


colnames(scores_uncentered) <- paste0("score_", seq_len(k_retain))
lme_dt <- cbind(covariates_dt, scores_uncentered)
head(lme_dt)




# Fixed Effects Formula for the Model: ------------------------------------
model_formula_fixef <- "ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent"

lme_fit <- fit_lme_to_scores(
  lme_df = lme_dt,
  k_retain = k_retain,
  fixef_formula = model_formula_fixef, 
  ranef_formula =  "+ (1|subject_id)",
  REML = TRUE)


# Save results: -----------------------------------------------------------
saveRDS(object = lme_fit, file.path(results_path, "model-fit-results.rds"))
