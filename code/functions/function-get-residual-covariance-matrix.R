get_residual_covariance_matrix <- function(x) {
  require(nlme)
  require(mixedup)
  stopifnot(class(x) == "lme")
  resid_het_var <- extract_het_var(x, digits = 16)
  resid_het_sd <- sqrt(resid_het_var)
  resid_cor_matrix <- as.matrix(extract_cor_structure(x, digits = 16))
  diag(resid_het_sd) %*% resid_cor_matrix %*% diag(resid_het_sd)
}
