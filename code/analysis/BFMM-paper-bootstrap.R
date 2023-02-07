# ------------------------------------------------------------------------#
# Bootstrap (and Wald) Intervals for Approximate Fixed-Effects Inference  #
# ------------------------------------------------------------------------#

# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(fda)        # CRAN v5.5.1
library(lme4)       # CRAN v1.1-30

# Load helper functions: --------------------------------------------------
source(here::here("code", "functions", "BFMM-paper-helper-functions.R"))

ncores <- parallel::detectCores()
B <- 2500 # No. of Bootstrap Replicates (Degras, 2017)

# Load results: -----------------------------------------------------------
results_path <- here::here("outputs", "results")

basis_transformation_results <- readRDS(
  file = file.path(results_path, "basis-transform-results.rds"))

lme_fit <- readRDS(file.path(results_path, "model-fit-results.rds"))


# Extract objects from modelling results: ---------------------------------
# Extract things from fitted model:
model_formula_fixef <- lme_fit$fixef_formula
lme_dt <- lme_fit$lme_df
bfpca <- basis_transformation_results$bfpca_obj
bfd_obj <- basis_transformation_results$bfd_obj
k_retain <- basis_transformation_results$k_retain
Psi <- rbind(eval.fd(0:100, bfpca$harmonics[,1]),
             eval.fd(0:100, bfpca$harmonics[,2]))
fixef_coef_point <- lme_fit$fixef_mat
fixef_coef_var <- lme_fit$fixef_var
fixef_names <- names(fixef(lme_fit$full_list[[1]]))
names(fixef_names) <- fixef_names



# Do bootstrap of subjects: -----------------------------------------------
set.seed(1)
# and store time:
boot_time <- system.time(
  boot_result <- bootstrap_of_subjects_coefs(
    lme_df = lme_dt,
    fixef_formula = model_formula_fixef,
    k_retain = k_retain,
    REML = TRUE,
    ranef_formula = "+ (1|subject_id)",
    B = B, # Degras (2017) 
    par_mc = TRUE,
    n_cores = ncores)
)


# Extract results: --------------------------------------------------------

# Functional Point Estimates = coefficient point-estimates x basis functions. 
fixef_fun_point <- fixef_coef_point %*% t(Psi)

# Construct estimated covariance function for each functional fixed effect.
# Diagonal covariance matrix for coefficients reflects that each coefficient
# is modelled independently.
fixef_fun_covar <- lapply(fixef_names, function(x) {
  Psi %*% diag(fixef_coef_var[x, ]) %*% t(Psi)
})

# Estimated standard-error function for each fixed effect.
# I.e., the pointwise standard deviation.
fixef_fun_se <- t(sapply(fixef_fun_covar, function(x) {
  sqrt(diag(x))
}))

# Extract the 
fixef_coef_samples_boot <- lapply(
  fixef_names,
  FUN = function(y) {
    t(sapply(boot_result, function(x) {
      x$fixef[y, ]
    }))
  }
)


fixef_coef_covar_boot <- lapply(fixef_coef_samples_boot, var) 

fixef_fun_se_boot <- lapply(fixef_names, function(x) {
  sqrt(diag(Psi %*% fixef_coef_covar_boot[[x]] %*% t(Psi)))
})


sim_cb_bootstrap <- lapply(X = fixef_names, FUN = function(x) {
  mvn_sim(coef_point_est = fixef_coef_point[x, , drop = TRUE],
          coef_covar_mat = fixef_coef_covar_boot[[x]],
          Psi_basis = Psi,
          N_simulation_mvn = 10000, 
          coverage_level = 0.95)
})

sim_cb_wald <- lapply(X = fixef_names, FUN = function(x) {
  mvn_sim(coef_point_est = fixef_coef_point[x, , drop = TRUE],
          coef_covar_mat = diag(fixef_coef_var[x, , drop = TRUE]),
          Psi_basis = Psi,
          N_simulation_mvn = 10000, 
          coverage_level = 0.95)})



# Put into Data Table for Visualising Results: ----------------------------
parameter_results_df <- purrr::map_dfr(fixef_names, function(x) {
  data.frame(
    dimension = rep(c("hip", "knee"), each = 101),
    t = rep(0:100, times = 2),
    point_est = fixef_fun_point[x,, drop = TRUE],
    se_boot = fixef_fun_se_boot[[x]],
    se_wald = fixef_fun_se[x, ],
    sim_boot_lower = sim_cb_bootstrap[[x]]$lower,
    sim_boot_upper = sim_cb_bootstrap[[x]]$upper,
    sim_wald_lower = sim_cb_wald[[x]]$lower,
    sim_wald_upper = sim_cb_wald[[x]]$upper
  )
}, .id = "parameter")

parameter_results_dt <- data.table(parameter_results_df)

parameter_results_dt[
  , `:=`(
    pw_wald_lower = point_est - 2 * se_wald,
    pw_wald_upper = point_est + 2 * se_wald,
    pw_boot_lower = point_est - 2 * se_boot,
    pw_boot_upper = point_est + 2 * se_boot
  )
]



# -------------------------------------------------------------------------

# Save results: -----------------------------------------------------------
saveRDS(object = list(parameter_results_dt = parameter_results_dt,
                      boot_time = boot_time, 
                      boot_result = boot_result),
  file.path(results_path, "model-fit-results.rds"))
