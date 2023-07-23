# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(ggplot2)    # CRAN v3.3.5
library(fda)        # CRAN v5.5.1
library(multifamm)  # CRAN v0.1.1
library(sparseFLMM) # CRAN v0.4.1
source("code/functions/FUI-functions/lfosr3s.R") # put your own paths to fui functions.
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


fui_hip_time_standard <- system.time(
  fui_hip_standard <- lfosr3s(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id), 
                     data = covariates_df[covariates_df$location == "Hip",], 
                     family = "gaussian",
                     var = TRUE,
                     analytic = FALSE,
                     parallel = TRUE)
)["elapsed"]

fui_knee_time_standard <- system.time(fui_knee_standard <- lfosr3s(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id), 
                                                 data = covariates_df[covariates_df$location == "Knee",], 
                                                 family = "gaussian",
                                                 var = TRUE,
                                                 analytic = FALSE,
                                                 parallel = TRUE))["elapsed"]


fui_hip_time_updated <- system.time(
  fui_hip_updated <- lfosr3s_updated(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id), 
                              data = covariates_df[covariates_df$location == "Hip",], 
                              family = "gaussian",
                              var = TRUE,
                              analytic = FALSE,
                              parallel = TRUE)
)["elapsed"]

fui_knee_time_updated <- system.time(fui_knee_updated <- lfosr3s_updated(formula = Y ~ ris + sex + speed_cent +  weight_kg_cent + height_cm_cent + age_cent + (1|subject_id), 
                                                          data = covariates_df[covariates_df$location == "Knee",], 
                                                          family = "gaussian",
                                                          var = TRUE,
                                                          analytic = FALSE,
                                                          parallel = TRUE))["elapsed"]



i<-1
par(mfrow = c(3, 3))
for(i in 1:9) {
  se_standard <- sqrt(diag(fui_hip_standard$betaHat.var[,,i]))
  est_standard <- fui_hip_standard$betaHat[i,]
  se_updated <- sqrt(diag(fui_hip_updated$betaHat.var[,,i]))
  est_updated <- fui_hip_updated$betaHat[i,]
  
  plot(est_standard, type = "l", 
       ylim = range(est_standard + 2 * se_standard,
                    est_standard - 2 * se_standard,
                    est_updated + 2 * se_updated,
                    est_updated - 2 * se_updated))
  lines(est_standard - 2 * se_standard, lty = 2)
  lines(est_standard + 2 * se_standard, lty = 2)
  title(paste("SE Ratio =", mean(se_updated / se_standard)))
  lines(est_updated, col = 2)
  lines(est_updated - 2 * se_updated, lty = 2, col = 2)
  lines(est_updated + 2 * se_updated, lty = 2, col = 2)
  abline(h = 0)
}

for(i in 1:9) {
  se_standard <- sqrt(diag(fui_knee_standard$betaHat.var[,,i]))
  est_standard <- fui_knee_standard$betaHat[i,]
  se_updated <- sqrt(diag(fui_knee_updated$betaHat.var[,,i]))
  est_updated <- fui_knee_updated$betaHat[i,]
  
  plot(est_standard, type = "l", 
       ylim = range(est_standard + 2 * se_standard,
                    est_standard - 2 * se_standard,
                    est_updated + 2 * se_updated,
                    est_updated - 2 * se_updated))
  title(paste("SE Ratio =", mean(se_updated / se_standard)))
  lines(est_standard - 2 * se_standard, lty = 2)
  lines(est_standard + 2 * se_standard, lty = 2)
  lines(est_updated, col = 2)
  lines(est_updated - 2 * se_updated, lty = 2, col = 2)
  lines(est_updated + 2 * se_updated, lty = 2, col = 2)
  abline(h = 0)
}


