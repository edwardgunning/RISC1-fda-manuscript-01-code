require(fda)      # CRAN v5.5.1
# Wrapper function top do bfpca on simulated data -------------------------
do_bfpca <- function(Y, argvals, pc_var_cutoff) {
  Y_hip <- Y[, 1:101]
  Y_knee <- Y[, 102:202]
  Y_array <- array(data = NA, dim = c(101, N_sub * J_rep, 2))
  Y_array[,, 1] <- t(Y_hip)
  Y_array[,, 2] <- t(Y_knee)
  Y_bfd <- Data2fd(argvals = argvals, y = Y_array)
  Y_basis <- Y_bfd$basis
  Y_bfpca <- pca.fd(Y_bfd, harmfdPar = fdPar(Y_basis))
  k_retain <- min(which(cumsum(Y_bfpca$values/ sum(Y_bfpca$values)) > pc_var_cutoff))
  bfpca <- pca.fd(fdobj = Y_bfd,
                  nharm = k_retain, 
                  harmfdPar = fdPar(Y_basis))
  psi_hat_array <- eval.fd(evalarg = argvals, fdobj = bfpca$harmonics)
  Psi_hat <- rbind(psi_hat_array[,, 1], psi_hat_array[,,2])
  scores <- apply(bfpca$scores, c(1, 2), sum)

  stopifnot(ncol(scores) == k_retain)
  stopifnot(ncol(Psi_hat) == k_retain)
  stopifnot(nrow(scores) == nrow(Y))
  stopifnot(nrow(Psi_hat) == 2 * length(argvals))
  # Return 
  list(k_retain = k_retain,
       scores = scores,
       Psi_hat = Psi_hat,
       bfpca = bfpca)
}

grid_to_bivariate_fd <- function(Y, argvals, basis_obj) {
  Y_hip <- Y[, 1:101]
  Y_knee <- Y[, 102:202]
  Y_array <- array(data = NA, dim = c(101, N_sub * J_rep, 2))
  Y_array[,, 1] <- t(Y_hip)
  Y_array[,, 2] <- t(Y_knee)
  Data2fd(argvals = argvals, y = Y_array, basisobj = basis_obj)
}


