results_path <- here::here("outputs", "results")
bootstrap_results <- readRDS(file = file.path(results_path, "bootstrap-results.rds"))

lme_fit <- readRDS(file.path(results_path, "model-fit-results.rds"))
q_vec <- lme_fit$q_vec
s_vec <- lme_fit$s_vec

icc_hat <- sum(q_vec) / sum(s_vec, q_vec)

icc_dist_boot <- sapply(bootstrap_results$boot_result, function(x) {
  sum(x$q_vec) / sum(x$s_vec, x$q_vec)
})

icc_hat_se <- sd(icc_dist_boot)
icc_ci_normal <- icc_hat + c(-2, 2) * icc_hat_se
icc_ci_quantile <- quantile(icc_dist_boot, probs = c(0.025, 0.975))

icc_confint_df <- data.frame(
  method = c("normal", "quantile"), 
  lower = round(c(icc_ci_normal[1], icc_ci_quantile[1]), 2),
  upper = round(c(icc_ci_normal[2], icc_ci_quantile[2]), 2))

data.table::fwrite(file = file.path(results_path, "icc_confint.csv"),
                   x = icc_confint_df)
