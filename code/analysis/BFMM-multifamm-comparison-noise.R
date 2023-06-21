# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(ggplot2)    # CRAN v3.3.5
library(fda)        # CRAN v5.5.1
library(multifamm)  # CRAN v0.1.1
library(sparseFLMM) # CRAN v0.4.1


source(file = here::here(
  "code",
  "functions",
  "functions-helper-smoothing.R"
))

data_path <- here::here("data")
subject_side_coef_reg <- readRDS(file.path(data_path, "subject_side_coef_reg.rds"))

# set up basis to represent the data.
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)

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

rm(subject_side_coef_reg)


# Prepare the data for fitting the multiFAMM ------------------------------
# reshape to long form to fit the multiFAMM
subject_side_eval_reg_lng <- melt.data.table(subject_side_eval_reg,
                                             measure.vars = paste0(0:100),
                                             variable.name = "t",
                                             variable.factor = FALSE,
                                             value.name = "y_vec")

subject_side_eval_reg_lng[, t := as.numeric(t)]

names(subject_side_eval_reg_lng)

mfamm_data <- subject_side_eval_reg_lng[location %in% c("Hip", "Knee"),
                                        c("subject_id", 
                                          "side",
                                          "location",
                                          "t", 
                                          "y_vec",
                                          "retrospective_injury_status",
                                          "sex",
                                          "speed_cent",
                                          "age_cent",
                                          "weight_kg_cent",
                                          "height_cm_cent")]



mfamm_data[, sex := fifelse(sex == "female", 1, 0)]
mfamm_data[, ris2 := fifelse(retrospective_injury_status == "injured_greater_than_2_yr", 1, 0)]
mfamm_data[, ris3 := fifelse(retrospective_injury_status == "injured_1_to_2_yr", 1, 0)]
mfamm_data[, ris4 := fifelse(retrospective_injury_status == "injured_less_than_1_yr", 1, 0)]
mfamm_data[, retrospective_injury_status := NULL]
mfamm_data[, subject_id := as.integer(as.factor(subject_id))]
mfamm_data[, n_long := as.integer(factor(interaction(side, subject_id)))]

setnames(mfamm_data,
         old = c("subject_id", "location", "sex", "speed_cent", 
                 "ris2", "ris3", "ris4", "age_cent",
                 "weight_kg_cent", "height_cm_cent"),
         new = c("subject_long", "dim",
                 paste0("covariate.", 1:8)))

mfamm_data[, dim := factor(dim)]
stopifnot(is.numeric(mfamm_data$y_vec))
stopifnot(is.numeric(mfamm_data$t))
stopifnot(is.factor(mfamm_data$dim))
stopifnot(is.integer(mfamm_data$subject_long))
stopifnot(is.integer(mfamm_data$n_long))
stopifnot(purrr::map_lgl(mfamm_data[, paste0("covariate.", 1:8)],
                         is.numeric))
setorderv(mfamm_data, c("dim", "subject_long", "n_long", "t"))
mfamm_data$y_vec <- mfamm_data$y_vec + rnorm(n = length(mfamm_data$y_vec), mean = 0, sd = 1)



system.time(mfamm_model <- multiFAMM(data = mfamm_data,
                                     covariate = TRUE,
                                     fRI_B = TRUE,
                                     fRI_C = FALSE,
                                     nested = FALSE,
                                     num_covariates = 8,
                                     min_grid = 0, 
                                     max_grid = 100,
                                     d_grid = 101,
                                     bf_mean = 35,
                                     m_mean = c(2, 2),
                                     bf_covariates = 35,
                                     my_grid = 0:100,
                                     covariate_form = rep("by", 8), 
                                     interaction = FALSE,
                                     which_interaction = matrix(NA),
                                     bf_covs = c(10, 10),
                                     m_covs = list(c(2, 3), c(2, 3)),
                                     var_level = 1,
                                     mfpc_cutoff = 0.95,
                                     mfpc_cut_method = "total_var",
                                     final_method = "bam"))
mfamm_model$mfpc$B
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



# -------------------------------------------------------------------------

ggplot(data = mfamm_hip_est_df) + 
  facet_wrap(~ param, scales = "free_y") +
  aes(x = t) +
  geom_line(aes(y = estimate)) +
  geom_line(aes(y = lower), lty = 2) +
  geom_line(aes(y = upper), lty = 2) +
  geom_hline(yintercept  = 0)




# -------------------------------------------------------------------------

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

ggplot(data = mfamm_knee_est_df) + 
  facet_wrap(~ param, scales = "free_y") +
  aes(x = t) +
  geom_line(aes(y = estimate)) +
  geom_line(aes(y = lower), lty = 2) +
  geom_line(aes(y = upper), lty = 2) +
  geom_hline(yintercept  = 0)


