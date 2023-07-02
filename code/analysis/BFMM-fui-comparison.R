# -------------------------------------------------------------------------
# Compare point estimates and standard errors to 
# fast univatiate inference methods by Cui (2022).
# note I have modified the bootstrap approach in lfosr3s to  use
# pseudo IDs in cluster labels.
# -------------------------------------------------------------------------

# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(ggplot2)    # CRAN v3.3.5
library(fda)        # CRAN v5.5.1
library(multifamm)  # CRAN v0.1.1
library(sparseFLMM) # CRAN v0.4.1
library(xtable)     # CRAN v1.8-4
source("code/functions/FUI-functions/lfosr3s.R")
source("code/functions/FUI-functions/lfosr3s-updated.R")

source(file = here::here(
  "code",
  "functions",
  "functions-helper-smoothing.R"
))
source(here::here(
  "code/functions/rough_fit_mfamm_model.R"
))

data_path <- here::here("data")
subject_side_coef_reg <- readRDS(file.path(data_path, "subject_side_coef_reg.rds"))

# set up basis to represent the data.
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)


# Prepare data for FUI: ---------------------------------------------------
subject_side_eval_reg <- copy(subject_side_coef_reg)
subject_side_eval_reg[, paste0(0:100) := {
  coef_mat <- coef_to_mat(.SD)
  fd_object <- fd(coef = coef_mat, basisobj = bspl80)
  fd_eval <- eval.fd(evalarg = 0:100, fdobj = fd_object)
  return_list <- get_list_of_rows(fd_eval)
  names(return_list) <- 0:100
  return_list},
  .SDcols = paste0("lm_coef_", 1:80)]
subject_side_eval_reg[, paste0("lm_coef_", 1:80) := NULL]
subject_side_eval_reg <- subject_side_eval_reg[location %in% c("Hip", "Knee")]
covariates_df <- data.frame(subject_side_eval_reg[, .(location = location, ris = factor(retrospective_injury_status, ordered = FALSE), sex, speed_cent, weight_kg_cent, height_cm_cent, age_cent, subject_id)])
covariates_df$Y <- as.matrix(subject_side_eval_reg[, paste0(0:100)])

fui_hip_analytic_time <- system.time(
  fui_hip_analytic <- lfosr3s(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id), 
                     data = covariates_df[covariates_df$location == "Hip",], 
                     family = "gaussian",
                     var = TRUE,
                     analytic = TRUE)
)["elapsed"]


fui_knee_analytic_time <- system.time(
  fui_knee_analytic <- lfosr3s(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id), 
                   data = covariates_df[covariates_df$location == "Knee",], 
                   family = "gaussian",
                   var = TRUE,
                   analytic = TRUE))["elapsed"]


fui_hip_boot_time <- system.time(
  fui_hip_boot <- lfosr3s_updated(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id),
                     data = covariates_df[covariates_df$location == "Hip",],
                     family = "gaussian",
                     var = TRUE,
                     analytic = FALSE,
                     parallel = TRUE))["elapsed"]


fui_knee_boot_time <- system.time(
  fui_knee_boot <- lfosr3s_updated(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id),
                             data = covariates_df[covariates_df$location == "Knee",],
                             family = "gaussian",
                             var = TRUE,
                             analytic = FALSE,
                             parallel = TRUE))["elapsed"]
fui_hip_analytic_time


fui_df <-data.frame(
  Approach = rep(c("Analytic", "Bootstrap"), each = 2),
  Dimension = rep(c("Hip", "Knee"), times = 2),
  time = c(fui_hip_analytic_time,
                     fui_knee_analytic_time,
                     fui_hip_boot_time, 
                     fui_knee_boot_time))
colnames(fui_df)[3] <- "Time (secs)"
bold <- function(x) {
  paste0("{\\bfseries ", x, "}") 
}
fui_table <- xtable(fui_df, 
                     digits =  c(0, 0, 0, 2), 
                     label = "tab:fui-comp-time",
                     caption = "Computation time for fitting separate FUI models to the hip and knee data using bootstrap and analytic approaches.")
align(fui_table)[1] <- "l"
print(fui_table, 
      file = here::here("outputs", "tables", "fui-comp-time.tex"),
      sanitize.text.function = function(x){x},
      sanitize.colnames.function = bold,
      include.rownames = FALSE,
      booktabs = TRUE)


t <- 0:100
fui_results_df <- data.frame()
for(i in 1:9) {
  param <- rownames(fui_hip_analytic$betaHat)[i]
  fui_point_est <- fui_hip_analytic$betaHat[i,]
  
  se_boot <- sqrt(diag((1/1.2) * fui_hip_boot$betaHat.var[,,i]))
  se_analytic <- sqrt(diag((1/1) * fui_hip_analytic$betaHat.var[,,i]))
  
  df_i <- data.frame(t = t,
             dimension = "hip",
             parameter = param,
             fui_point_est = fui_point_est,
             pw_fui_lower_boot = fui_point_est - 2 * se_boot,
             pw_fui_upper_boot = fui_point_est + 2 * se_boot,
             pw_fui_lower_analytic = fui_point_est - 2 * se_analytic,
             pw_fui_upper_analytic = fui_point_est + 2 * se_analytic)
  fui_results_df <- rbind(fui_results_df, df_i)
}
for(i in 1:9) {
  param <- rownames(fui_knee_analytic$betaHat)[i]
  fui_point_est <- fui_knee_analytic$betaHat[i,]
  
  se_boot <- sqrt(diag((1/1.2) * fui_knee_boot$betaHat.var[,,i]))
  se_analytic <- sqrt(diag((1/1) * fui_knee_analytic$betaHat.var[,,i]))
  
  df_i <- data.frame(t = t,
                     dimension = "knee",
                     parameter = param,
                     fui_point_est = fui_point_est,
                     pw_fui_lower_boot = fui_point_est - 2 * se_boot,
                     pw_fui_upper_boot = fui_point_est + 2 * se_boot,
                     pw_fui_lower_analytic = fui_point_est - 2 * se_analytic,
                     pw_fui_upper_analytic = fui_point_est + 2 * se_analytic)
  fui_results_df <- rbind(fui_results_df, df_i)
}

saveRDS(object = fui_results_df, file = here::here("outputs", "results", "fui_comparison.rds"))


