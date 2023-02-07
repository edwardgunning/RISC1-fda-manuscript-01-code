# -------------------------------------------------------------------------
# Test Multilevel FPCA Covariance Estimation for MULTIVARIATE functional
# data:
# Construct bivariate Eigenfunctions by the "weighting" method
# of Happ & Greven (2018)
# -------------------------------------------------------------------------
source(here::here("code", 
                  "functions",
                  "functions-unstructured-covariance.R"))
library(Matrix)    # CRAN v1.4-0
library(funData)   # CRAN v1.3-8
library(fda)       # CRAN v5.5.1

# ------------------------------------------------------------------------#
# Generate data from a univariate functional mixed model
# ------------------------------------------------------------------------#

N_tot <- 300 # no. of subjects
J <- 2 # no. of replicate curves per subject
sampling_points <- 0:100

# generate dataset
subject_id <- rep(seq_len(N_tot), each = J)
replicate <- rep(seq_len(J), times = N_tot)
df <- data.frame(factor(subject_id), factor(replicate))

# design matrices for random effects
Z1 <- sparse.model.matrix(~ - 1 + factor(subject_id), data = df)

M <- 10 # number of FPCs, 
# used to generate the functional random intercepts and curve-level random intercepts
# Generate Eigenfunctions. 
# Use different bases for each level.


# -------------------------------------------------------------------------
# generate bivariate basis functions to be used in simulation:
sim_bfundata_U <- funData::simMultiFunData(argvals = list(list(0:100), list(0:100)), 
                                           eFunType = list("Fourier", "Poly"),
                                           M = list(M, M),
                                           N = 1,
                                           eValType = "linear",
                                           type = "weighted")

efuns_U_hip <- fda::eval.fd(0:100, funData::funData2fd(sim_bfundata_U$trueFuns[[1]]))
efuns_U_knee <- fda::eval.fd(0:100, funData::funData2fd(sim_bfundata_U$trueFuns[[2]]))
efuns_U <- rbind(efuns_U_hip, efuns_U_knee)


sim_bfundata_E <- funData::simMultiFunData(argvals = list(list(0:100), list(0:100)), 
                                           eFunType = list("Wiener", "Poly"),
                                           M = list(M, M),
                                           N = 1,
                                           # ignoreDeg = list(1:2, 0),
                                           eValType = "linear",
                                           type = "weighted")

efuns_poly <- sim_bfundata_E$trueFuns
efuns_E_hip <- fda::eval.fd(0:100, funData::funData2fd(efuns_poly[[1]]))
efuns_E_knee <- fda::eval.fd(0:100, funData::funData2fd(efuns_poly[[2]]))
efuns_E <- rbind(efuns_E_hip, efuns_E_knee)
# -------------------------------------------------------------------------


# Generate eigenvalues
scale_wb <- 0.5 # ratio of subject-level evals and error term evals
eval1 <- 0.5 ^ (seq_len(M) - 1) # generate evals as 2^{m-1}
eval2 <- scale_wb * ( 0.5 ^ (seq_len(M) - 1) )


# -------------------------------------------------------------------------

# draw random effects scores.
set.seed(1996) 


N_sim <- 2000
K_U_11_res <- K_U_12_res <- K_U_21_res <- K_U_22_res <- array(NA,dim = c(101, 101, N_sim))
K_E_11_res <- K_E_12_res <- K_E_21_res <- K_E_22_res <- array(NA,dim  = c(101, 101, N_sim))
for(i in seq_len(N_sim)) {
  print(paste0("Iteration:", i))
  # -------------------------------------------------------------------------
  # center and scale them so covariance estimates will be almost exact.
  U <- mvtnorm::rmvnorm(n = N_tot,  sigma = diag(eval1)) # generate random scores for subject-level random intercepts
  E <- mvtnorm::rmvnorm(n = N_tot * J,  sigma = diag(eval2)) # generate random scores for curve-level random intercepts
  Y_grid <- Z1 %*% U %*% t(efuns_U) + E %*% t(efuns_E)
  bspl35 <- create.bspline.basis(rangeval = c(0, 100), nbasis = 80, norder = 4)
  Y_hip_fd <- Data2fd(argvals = sampling_points, y = as.matrix(t(Y_grid[, 1:101])), basisobj = bspl35)
  Y_knee_fd <- Data2fd(argvals = sampling_points, y = as.matrix(t(Y_grid[, 102:202])), basisobj = bspl35)
  test <- cov_unstruc_mlfpca_bi_fd(fd_obj_list = list(Y_hip_fd, Y_knee_fd),
                                   id_vec = subject_id)
  # -------------------------------------------------------------------------
  K_U_11_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_U_bifd_list$K_U_11)
  K_U_12_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_U_bifd_list$K_U_12)
  K_U_21_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_U_bifd_list$K_U_21)
  K_U_22_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_U_bifd_list$K_U_22)
  # -------------------------------------------------------------------------
  K_E_11_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_E_bifd_list$K_E_11)
  K_E_12_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_E_bifd_list$K_E_12)
  K_E_21_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_E_bifd_list$K_E_21)
  K_E_22_res[,,i] <- eval.bifd(sevalarg = 0:100,
                               tevalarg = 0:100,
                               bifd = test$K_E_bifd_list$K_E_22)
}



# ------------------------------------------------------------------------#
# Get truth: --------------------------------------------------------------
K_U_11_true <- efuns_U[1:101, ] %*% diag(eval1) %*% t(efuns_U[1:101, ])
K_U_12_true <- efuns_U[1:101, ] %*% diag(eval1) %*% t(efuns_U[102:202, ])
K_U_21_true <- efuns_U[102:202, ] %*% diag(eval1) %*% t(efuns_U[1:101, ])
K_U_22_true <- efuns_U[102:202, ] %*% diag(eval1) %*% t(efuns_U[102:202, ])

K_E_11_true <- efuns_E[1:101, ] %*% diag(eval2) %*% t(efuns_E[1:101, ])
K_E_12_true <- efuns_E[1:101, ] %*% diag(eval2) %*% t(efuns_E[102:202, ])
K_E_21_true <- efuns_E[102:202, ] %*% diag(eval2) %*% t(efuns_E[1:101, ])
K_E_22_true <- efuns_E[102:202, ] %*% diag(eval2) %*% t(efuns_E[102:202, ])
# ------------------------------------------------------------------------#


# -------------------------------------------------------------------------
# True Parameter vs E[Estimate]:
filled.contour(K_U_11_true)
filled.contour(apply(K_U_11_res, c(1, 2), mean))

filled.contour(K_U_12_true)
filled.contour(apply(K_U_12_res, c(1, 2), mean))

filled.contour(K_U_21_true)
filled.contour(apply(K_U_21_res, c(1, 2), mean))

filled.contour(K_U_22_true)
filled.contour(apply(K_U_22_res, c(1, 2), mean))


# -------------------------------------------------------------------------
filled.contour(K_E_11_true)
filled.contour(apply(K_E_11_res, c(1, 2), mean))

filled.contour(K_E_12_true)
filled.contour(apply(K_E_12_res, c(1, 2), mean))

filled.contour(K_E_21_true)
filled.contour(apply(K_E_21_res, c(1, 2), mean))

filled.contour(K_E_22_true)
filled.contour(apply(K_E_22_res, c(1, 2), mean))


# all look good.

# Can look at bias, too.
# In U: -------------------------------------------------------------------
K_U_11_bias <- apply(sweep(K_U_11_res,
                           MARGIN = c(1, 2),
                           STATS = K_U_11_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_U_11_bias)

K_U_12_bias <- apply(sweep(K_U_12_res,
                           MARGIN = c(1, 2),
                           STATS = K_U_12_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_U_12_bias)
mean(abs(K_U_12_bias))

K_U_21_bias <- apply(sweep(K_U_21_res,
                           MARGIN = c(1, 2),
                           STATS = K_U_21_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_U_21_bias)
mean(abs(K_U_21_bias))

K_U_22_bias <- apply(sweep(K_U_22_res,
                           MARGIN = c(1, 2),
                           STATS = K_U_22_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_U_22_bias)
mean(abs(K_U_11_bias))

# In E: -------------------------------------------------------------------
K_E_11_bias <- apply(sweep(K_E_11_res,
                           MARGIN = c(1, 2),
                           STATS = K_E_11_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_E_11_bias)

K_E_12_bias <- apply(sweep(K_E_12_res,
                           MARGIN = c(1, 2),
                           STATS = K_E_12_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_E_12_bias)
mean(abs(K_E_12_bias))

K_E_21_bias <- apply(sweep(K_E_21_res,
                           MARGIN = c(1, 2),
                           STATS = K_E_21_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_E_21_bias)
mean(abs(K_E_21_bias))

K_E_22_bias <- apply(sweep(K_E_22_res,
                           MARGIN = c(1, 2),
                           STATS = K_E_22_true,
                           FUN = "-"), c(1,2), mean)
filled.contour(K_E_22_bias)
mean(abs(K_E_11_bias))

# Bias is overall small... 


