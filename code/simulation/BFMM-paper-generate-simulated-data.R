# -------------------------------------------------------------------------
# Functions to generate the bivariate functional data.
# in two scenarios in the paper!
# -------------------------------------------------------------------------

results_path <- here::here("outputs", "results")
sim_data_list <- readRDS(file.path(results_path, "BFMM-simulation-parameters.rds"))
efuns_s2 <- readRDS(file.path(results_path, "BFMM-simulation-efuns.rds"))

# Scenario 1 & 2: Empirical Fixed Effects Functions Used
generate_data_scenario_1 <- function(N = 280, J = 2) {
  
  # generate non-functional part of dataset
  subject_id <- factor(rep(seq_len(length.out = N), each = J))
  side <- rep(seq_len(J), times = N)
  df <- data.frame(subject_id, side)
  covariate_speed <- rnorm(n = N, sd = sim_data_list$speed_sd)
  df$speed <- rep(covariate_speed, each = J)
  covariate_sex <- factor(ifelse(runif(N) > 0.5, "male", "female"), levels = c("male", "female"))
  df$sex <- rep(covariate_sex, each = J)
  X <- model.matrix( ~ sex + speed, data = df) # don't need intercept
  Z <- model.matrix(~ subject_id - 1, data = df)
  
  # Functional Part
  U_star <- mvtnorm::rmvnorm(n = N, sigma = diag(sim_data_list$q_vec))
  E_star <- mvtnorm::rmvnorm(n = N * J, sigma = diag(sim_data_list$s_vec))
  # Scalar model
  Y_star <- X %*% sim_data_list$B_empirical_scores + Z %*% U_star + E_star
  Y_uc <- Y_star %*% t(sim_data_list$Phi) # multiply by basis
  Y <- sweep(Y_uc, MARGIN = 2, STATS = sim_data_list$mean_eval_vec,FUN = "+", check.margin = TRUE) # add back on mean
  list(df = df,
       Y = Y,
       Y_star = Y_star)
 }



generate_data_scenario_2 <- function(N = 280, J = 2) {
  
  # generate non-functional part of dataset
  subject_id <- factor(rep(seq_len(length.out = N), each = J))
  side <- rep(seq_len(length.out = J), times = N)
  df <- data.frame(subject_id, side)
  covariate_speed <- rnorm(n = N, sd = sim_data_list$speed_sd) # draw speed from a normal dist w/ mean zero and empirical sd
  df$speed <- rep(covariate_speed, each = J)
  covariate_sex <- factor(ifelse(runif(N) > 0.5, "male", "female"), levels = c("male", "female")) # binary sex variable
  df$sex <- rep(covariate_sex, each = J)
  X <- model.matrix( ~ sex + speed, data = df) # don't need intercept
  Z <- model.matrix(~ subject_id - 1, data = df)

  
  # Functional Part
  U_star <- mvtnorm::rmvnorm(n = N, sigma = diag(sim_data_list$q_vec))
  E_star <- mvtnorm::rmvnorm(n = N * J, sigma = diag(sim_data_list$s_vec))

  # Scalar model
  Y_uc <- X %*% sim_data_list$B_empirical_scores %*% t(sim_data_list$Phi) + 
    Z %*% U_star %*% t(efuns_s2$efuns_U) +
    E_star %*% t(efuns_s2$efuns_E)
  Y <- sweep(Y_uc, MARGIN = 2, STATS = sim_data_list$mean_eval_vec,FUN = "+", check.margin = TRUE) # add back on mean
  list(df = df,
       Y = Y)
}


