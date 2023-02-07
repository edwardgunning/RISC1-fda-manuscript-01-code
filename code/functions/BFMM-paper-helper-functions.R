fit_lme_to_scores <- function(lme_df,
                              k_retain,
                              ranef_formula = "+ (1|subject_id)",
                              fixef_formula,
                              REML = TRUE) {
  require(lme4) 
  stopifnot(is.logical(REML))
  stopifnot(is.data.frame(lme_df))
  stopifnot(all(paste0("score_", seq_len(k_retain)) %in% names(lme_df)))
  
  lme_list <- lapply(seq_len(k_retain), function(k) {
    formula_k <- formula(
      paste0("score_", k, " ~ ", fixef_formula, ranef_formula))
    lmer(formula = formula_k, data = lme_df, REML = REML)
  })
  
  fixef_mat <- sapply(lme_list, fixef)
  fixef_var <- sapply(lme_list, function(x) {
    diag(vcov(x))
  })
  rownames(fixef_var) <- rownames(fixef_mat)
  q_vec <- sapply(lme_list, function(x) VarCorr(x)$subject_id[1])
  s_vec <- sapply(lme_list, function(x) {
    attr(VarCorr(x), "sc")^2
  })
  
  list(full_list = lme_list,
       fixef_formula = fixef_formula,
       lme_df = lme_df,
       fixef_mat = fixef_mat,
       fixef_var = fixef_var,
       q_vec = q_vec, 
       s_vec = s_vec)
}



# -------------------------------------------------------------------------
# Optimised version which only returns fixed effects
# parameters for bootstrap
fit_lme_to_scores_boot <- function(lme_df,
                              k_retain,
                              ranef_formula = "+ (1|subject_id)",
                              fixef_formula,
                              REML = TRUE) {
  require(lme4) 
  stopifnot(is.logical(REML))
  stopifnot(is.data.frame(lme_df))
  stopifnot(all(paste0("score_", seq_len(k_retain)) %in% names(lme_df)))
  
  lme_list <- lapply(seq_len(k_retain), function(k) {
    formula_k <- formula(
      paste0("score_", k, " ~ ", fixef_formula, ranef_formula))
    lmer(formula = formula_k, data = lme_df, REML = REML)
  })
  
  fixef_mat <- sapply(lme_list, fixef)
  fixef_var <- sapply(lme_list, function(x) {
    diag(vcov(x))
  })
  rownames(fixef_var) <- rownames(fixef_mat)
  q_vec <- sapply(lme_list, function(x) VarCorr(x)$subject_id[1])
  s_vec <- sapply(lme_list, function(x) {
    attr(VarCorr(x), "sc")^2
  })
  
  list(fixef_mat = fixef_mat,
       q_vec = q_vec, 
       s_vec = s_vec)
}
  

# -------------------------------------------------------------------------

bootstrap_of_subjects_coefs <- function(lme_df,
                                        k_retain,
                                        fixef_formula,
                                        ranef_formula = "+ (1|subject_id)",
                                        REML = TRUE,
                                        B = 1000,
                                        par_mc = TRUE,
                                        n_cores = 4) {
  require(parallel)
  # checks input
  stopifnot("subject_id" %in% names(lme_df))
  stopifnot(paste0("score_", seq_len(k_retain)) %in% names(lme_df))
  
  
  # Generate resampled indices before doing parralelisation to avoid RNG problems
  resampled_ids_list <- lapply(seq_len(B), function(b) {
    sample(unique(lme_df$subject_id), replace = TRUE)
  })
  
  # do bootstrapping in parralell
  if(!par_mc) {
    bootstrap_list <- lapply(resampled_ids_list, function(x) {
      df_b <- purrr::map_dfr(.x = seq_along(x), .f = function(i) {
        lme_df[lme_df$subject_id == x[i], ]
      }, .id = "subject_id_b")
      df_b$subject_id <- factor(df_b$subject_id_b)
      lme_fit_b <- fit_lme_to_scores_boot(lme_df = df_b, 
                        k_retain = k_retain,
                        REML = REML,
                        fixef_formula = fixef_formula,
                        ranef_formula = ranef_formula)
      list(fixef = lme_fit_b$fixef_mat, q_vec = lme_fit_b$q, s_vec = lme_fit_b$s)
    })
  }
  
  if(par_mc) {
    bootstrap_list <- mclapply(resampled_ids_list, function(x) {
      df_b <- purrr::map_dfr(.x = seq_along(x), .f = function(i) {
        lme_df[lme_df$subject_id == x[i], ]
      }, .id = "subject_id_b")
      df_b$subject_id <- factor(df_b$subject_id_b)
      lme_fit_b <- fit_lme_to_scores_boot(lme_df = df_b, 
                        k_retain = k_retain,
                        REML = REML,
                        fixef_formula = fixef_formula,
                        ranef_formula = ranef_formula)
      list(fixef = lme_fit_b$fixef_mat, q_vec = lme_fit_b$q, s_vec = lme_fit_b$s)
    }, mc.silent = FALSE, mc.cores = n_cores)
  }
  
  return(bootstrap_list)
}


# -------------------------------------------------------------------------

mvn_sim <- function(coef_point_est, # Vector of coefficients for the parameter function (length = k)
                    coef_covar_mat, # Estimated covariance matrix for coef_point_est (dimension = k x k)
                    Psi_basis, # Basis to which the coefficieients in coef_point_est correspond (dimension = npts x k)
                    N_simulation_mvn = 10000, # Number of samples to draw from MVN distribution.
                    coverage_level = 0.95) { # coverage level for the two-sided CI. 
  
  k <- length(coef_point_est)
  
  # Checks:
  stopifnot(is.vector(coef_point_est))
  stopifnot(k == ncol(Psi_basis))
  stopifnot(dim(coef_covar_mat) == c(k, k))
  # Check this is a proper covariance matrix being supplied:
  stopifnot(isSymmetric(coef_covar_mat))
  stopifnot(matrixcalc::is.positive.definite(coef_covar_mat))
  stopifnot((coverage_level > 0) & (coverage_level < 1))
  
  npts <- nrow(Psi_basis)
  
  # Get pointwise point estimate:
  fun_point_est <- (coef_point_est %*% t(Psi_basis))[1, ]
  
  # and pointwise standard error:
  fun_se <- sqrt(diag(
    Psi_basis %*% coef_covar_mat %*% t(Psi_basis)
  ))
  
  # Simulate from distribution of coefficients:
  coefs_samples <- mvtnorm::rmvnorm(n = N_simulation_mvn,
                                    mean = coef_point_est, 
                                    sigma = coef_covar_mat)
  # Convert simulated vectors to functions # by multiplying them by basis functions:
  fun_samples <- coefs_samples %*% t(Psi_basis)
  
  # Center samples on fixed effects point estimate
  fun_samples_cent <- sweep(fun_samples,
                            MARGIN = 2, 
                            STATS = fun_point_est,
                            FUN = "-",
                            check.margin = TRUE)
  
  # Now also scale these by pointwise standard error,
  # to give stanardised curves:
  fun_samples_stand <- sweep(fun_samples_cent,
                             MARGIN = 2, 
                             STATS = fun_se,
                             FUN = "/",
                             check.margin = TRUE)
  
  # calculate the maximum T statistic
  maxT_dist <- apply(X = fun_samples_stand,
                     MARGIN = 1,
                     FUN = function(x) max(abs(x)))
  
  # Obtain quantile of distribution:
  maxT_q <- quantile(maxT_dist, probs = coverage_level)
  
  
  if(maxT_q < 2) warning("Something is weird -- Simultaneous intervals narrower than pointwise!")
  
  lower <- fun_point_est - maxT_q * fun_se
  upper <- fun_point_est + maxT_q * fun_se
  
  list(lower = lower, upper = upper, q = maxT_q, fun_se = fun_se)
}




