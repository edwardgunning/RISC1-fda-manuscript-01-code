fit_multifamm <- function(data, 
                          bf_covs = c(5, 5),
                          mfpc_cutoff = 0.95) {
  gc() # collect garbage help memory:
  try(mfamm_model <- multiFAMM(data = data,
                                       covariate = TRUE,
                                       fRI_B = TRUE,
                                       fRI_C = FALSE,
                                       nested = FALSE,
                                       num_covariates = 8,
                                       min_grid = 0, 
                                       max_grid = 100,
                                       d_grid = 101,
                                       bf_mean = 26, # min{35, n/4} (Ruppert)
                                       m_mean = c(2, 2),
                                       bf_covariates = 26, # min{35, n/4} (Ruppert)
                                       my_grid = 0:100,
                                       covariate_form = rep("by", 8), 
                                       interaction = FALSE,
                                       which_interaction = matrix(NA),
                                       bf_covs = bf_covs,
                                       m_covs = list(c(2, 3), c(2, 3)),
                                       var_level = 1,
                                       mfpc_cutoff = mfpc_cutoff,
                                       mfpc_cut_method = "total_var",
                                       final_method = "bam"))
  
  newdata_hip <- mfamm_model$data[subject_long == 1 & n_long == 1 & dim == "Hip"]
  newdata_hip[, c(paste0("covariate.", 1:8),
                  paste0("dimHip.", 1:8))] <- 1
  terms_hip <- predict(mfamm_model$model, newdata = newdata_hip, type = "terms", se.fit = TRUE,
                       terms = c("dim", "s(t):dimHip",
                                 paste0("s(t):dimHip.", 1:8))
  )
  
  head(terms_hip$fit)
  mfamm_hip_est_df <- data.frame(t = rep(0:100, times = 9))
  mfamm_hip_est_df$param <- rep(c("intercept",
                                  "speed_cent_coefficient",
                                  "ris_effect[2]", 
                                  "ris_effect[3]",
                                  "ris_effect[4]",
                                  "sex_effect[2]",
                                  "age_cent",
                                  "weight_cent",
                                  "height_cent"),
                                each = 101)
  mfamm_hip_est_df$estimate <- mfamm_hip_est_df$lower <- mfamm_hip_est_df$upper <- NA
  
  # Intercept:
  intercept_est <- terms_hip$fit[, 1] + terms_hip$fit[, 2]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "intercept", "estimate"] <- intercept_est
  intercept_lower <- (terms_hip$fit[, 1] - 2 * terms_hip$se.fit[, 1]) + (terms_hip$fit[, 2] - 2 * terms_hip$se.fit[, 2])
  intercept_upper <- (terms_hip$fit[, 1] + 2 * terms_hip$se.fit[, 1]) + (terms_hip$fit[, 2] + 2 * terms_hip$se.fit[, 2])
  mfamm_hip_est_df[mfamm_hip_est_df$param == "intercept", "lower"] <- intercept_lower
  mfamm_hip_est_df[mfamm_hip_est_df$param == "intercept", "upper"] <- intercept_upper
  
  
  # Speed:
  mfamm_hip_est_df[mfamm_hip_est_df$param == "speed_cent_coefficient", "estimate"] <- terms_hip$fit[, "s(t):dimHip.2"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "speed_cent_coefficient", "lower"] <- terms_hip$fit[, "s(t):dimHip.2"] - 2 * terms_hip$se.fit[, "s(t):dimHip.2"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "speed_cent_coefficient", "upper"] <- terms_hip$fit[, "s(t):dimHip.2"] + 2 * terms_hip$se.fit[, "s(t):dimHip.2"]
  
  # Sex
  head(mfamm_hip_est_df)
  mfamm_hip_est_df[mfamm_hip_est_df$param == "sex_effect[2]", "estimate"] <- terms_hip$fit[, "s(t):dimHip.1"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "sex_effect[2]", "lower"] <- terms_hip$fit[, "s(t):dimHip.1"] - 2 * terms_hip$se.fit[, "s(t):dimHip.1"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "sex_effect[2]", "upper"] <- terms_hip$fit[, "s(t):dimHip.1"] + 2 * terms_hip$se.fit[, "s(t):dimHip.1"]
  
  
  # ris_effect[2]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[2]", "estimate"] <- terms_hip$fit[, "s(t):dimHip.3"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[2]", "lower"] <- terms_hip$fit[, "s(t):dimHip.3"] - 2 * terms_hip$se.fit[, "s(t):dimHip.3"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[2]", "upper"] <- terms_hip$fit[, "s(t):dimHip.3"] + 2 * terms_hip$se.fit[, "s(t):dimHip.3"]
  
  # ris_effect[3]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[3]", "estimate"] <- terms_hip$fit[, "s(t):dimHip.4"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[3]", "lower"] <- terms_hip$fit[, "s(t):dimHip.4"] - 2 * terms_hip$se.fit[, "s(t):dimHip.4"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[3]", "upper"] <- terms_hip$fit[, "s(t):dimHip.4"] + 2 * terms_hip$se.fit[, "s(t):dimHip.4"]
  
  # ris_effect[4]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[4]", "estimate"] <- terms_hip$fit[, "s(t):dimHip.5"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[4]", "lower"] <- terms_hip$fit[, "s(t):dimHip.5"] - 2 * terms_hip$se.fit[, "s(t):dimHip.5"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "ris_effect[4]", "upper"] <- terms_hip$fit[, "s(t):dimHip.5"] + 2 * terms_hip$se.fit[, "s(t):dimHip.5"]
  
  
  # age
  mfamm_hip_est_df[mfamm_hip_est_df$param == "age_cent", "estimate"] <- terms_hip$fit[, "s(t):dimHip.6"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "age_cent", "lower"] <- terms_hip$fit[, "s(t):dimHip.6"] - 2 * terms_hip$se.fit[, "s(t):dimHip.6"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "age_cent", "upper"] <- terms_hip$fit[, "s(t):dimHip.6"] + 2 * terms_hip$se.fit[, "s(t):dimHip.6"]
  
  # weight
  mfamm_hip_est_df[mfamm_hip_est_df$param == "weight_cent", "estimate"] <- terms_hip$fit[, "s(t):dimHip.7"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "weight_cent", "lower"] <- terms_hip$fit[, "s(t):dimHip.7"] - 2 * terms_hip$se.fit[, "s(t):dimHip.7"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "weight_cent", "upper"] <- terms_hip$fit[, "s(t):dimHip.7"] + 2 * terms_hip$se.fit[, "s(t):dimHip.7"]
  
  # height
  mfamm_hip_est_df[mfamm_hip_est_df$param == "height_cent", "estimate"] <- terms_hip$fit[, "s(t):dimHip.8"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "height_cent", "lower"] <- terms_hip$fit[, "s(t):dimHip.8"] - 2 * terms_hip$se.fit[, "s(t):dimHip.8"]
  mfamm_hip_est_df[mfamm_hip_est_df$param == "height_cent", "upper"] <- terms_hip$fit[, "s(t):dimHip.8"] + 2 * terms_hip$se.fit[, "s(t):dimHip.8"]
  
  
  newdata_knee <- mfamm_model$data[subject_long == 1 & n_long == 1 & dim == "Knee"]
  
  newdata_knee[, c(paste0("covariate.", 1:5),
                   paste0("dimKnee.", 1:5))] <- 1
  terms_knee <- predict(mfamm_model$model,
                        newdata_knee, type = "terms",
                        se.fit = TRUE,
                        terms = c("dim", "s(t):dimKnee",
                                  paste0("s(t):dimKnee.", 1:8))
  )
  
  terms_knee$fit
  terms_knee
  mfamm_knee_est_df <- data.frame(t = rep(0:100, times = 9))
  mfamm_knee_est_df$param <- rep(c("intercept",
                                   "speed_cent_coefficient",
                                   "ris_effect[2]",
                                   "ris_effect[3]",
                                   "ris_effect[4]",
                                   "sex_effect[2]",
                                   "age_cent",
                                   "weight_cent",
                                   "height_cent"),
                                 each = 101)
  mfamm_knee_est_df$estimate <- mfamm_knee_est_df$lower <- mfamm_knee_est_df$upper <- NA
  
  
  # Intercept:
  intercept_est <- terms_knee$fit[, 1] + terms_knee$fit[, 2]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "intercept", "estimate"] <- intercept_est
  intercept_lower <- (terms_knee$fit[, 1] - 2 * terms_knee$se.fit[, 1]) + (terms_knee$fit[, 2] - 2 * terms_knee$se.fit[, 2])
  intercept_upper <- (terms_knee$fit[, 1] + 2 * terms_knee$se.fit[, 1]) + (terms_knee$fit[, 2] + 2 * terms_knee$se.fit[, 2])
  mfamm_knee_est_df[mfamm_knee_est_df$param == "intercept", "lower"] <- intercept_lower
  mfamm_knee_est_df[mfamm_knee_est_df$param == "intercept", "upper"] <- intercept_upper
  
  # Speed:
  
  mfamm_knee_est_df[mfamm_knee_est_df$param == "speed_cent_coefficient", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.2"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "speed_cent_coefficient", "lower"] <- terms_knee$fit[, "s(t):dimKnee.2"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.2"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "speed_cent_coefficient", "upper"] <- terms_knee$fit[, "s(t):dimKnee.2"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.2"]
  
  # Sex
  head(mfamm_knee_est_df)
  mfamm_knee_est_df[mfamm_knee_est_df$param == "sex_effect[2]", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.1"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "sex_effect[2]", "lower"] <- terms_knee$fit[, "s(t):dimKnee.1"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.1"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "sex_effect[2]", "upper"] <- terms_knee$fit[, "s(t):dimKnee.1"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.1"]
  
  
  # ris_effect[2]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[2]", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.3"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[2]", "lower"] <- terms_knee$fit[, "s(t):dimKnee.3"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.3"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[2]", "upper"] <- terms_knee$fit[, "s(t):dimKnee.3"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.3"]
  
  # ris_effect[3]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[3]", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.4"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[3]", "lower"] <- terms_knee$fit[, "s(t):dimKnee.4"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.4"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[3]", "upper"] <- terms_knee$fit[, "s(t):dimKnee.4"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.4"]
  
  # ris_effect[4]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[4]", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.5"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[4]", "lower"] <- terms_knee$fit[, "s(t):dimKnee.5"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.5"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "ris_effect[4]", "upper"] <- terms_knee$fit[, "s(t):dimKnee.5"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.5"]
  
  # age
  mfamm_knee_est_df[mfamm_knee_est_df$param == "age_cent", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.6"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "age_cent", "lower"] <- terms_knee$fit[, "s(t):dimKnee.6"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.6"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "age_cent", "upper"] <- terms_knee$fit[, "s(t):dimKnee.6"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.6"]
  
  # weight
  mfamm_knee_est_df[mfamm_knee_est_df$param == "weight_cent", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.7"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "weight_cent", "lower"] <- terms_knee$fit[, "s(t):dimKnee.7"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.7"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "weight_cent", "upper"] <- terms_knee$fit[, "s(t):dimKnee.7"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.7"]
  
  # height
  mfamm_knee_est_df[mfamm_knee_est_df$param == "height_cent", "estimate"] <- terms_knee$fit[, "s(t):dimKnee.8"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "height_cent", "lower"] <- terms_knee$fit[, "s(t):dimKnee.8"] - 2 * terms_knee$se.fit[, "s(t):dimKnee.8"]
  mfamm_knee_est_df[mfamm_knee_est_df$param == "height_cent", "upper"] <- terms_knee$fit[, "s(t):dimKnee.8"] + 2 * terms_knee$se.fit[, "s(t):dimKnee.8"]
  
  
  mfamm_knee_est_df$dim <- "knee"
  mfamm_hip_est_df$dim <- "hip"
  
  rbind(mfamm_knee_est_df, mfamm_hip_est_df)
}