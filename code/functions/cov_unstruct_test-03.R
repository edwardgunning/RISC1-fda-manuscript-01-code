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

N_tot <- 10000 # no. of subjects
J <- 10 # no. of replicate curves per subject
sampling_points <- 0:100

# generate dataset
subject_id <- rep(seq_len(N_tot), each = J)
replicate <- rep(seq_len(J), times = N_tot)
df <- data.frame(factor(subject_id), factor(replicate))

# design matrices for random effects
Z1 <- sparse.model.matrix(~ - 1 + factor(subject_id), data = df)

M <- 3 # number of FPCs, used to generate the functional random intercepts and curve-level random intercepts
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


# draw random effects scores.
set.seed(1996) 
# center and scale them so covariance estimates will be almost exact.
U <- mvtnorm::rmvnorm(n = N_tot,  sigma = diag(eval1)) # generate random scores for subject-level random intercepts
U <- apply(U, 2, scale, center = TRUE, scale = TRUE)
U <- sweep(U, MARGIN = 2, STATS = sqrt(eval1), FUN = "*")


E <- mvtnorm::rmvnorm(n = N_tot * J,  sigma = diag(eval2)) # generate random scores for curve-level random intercepts
E <- apply(E, 2, scale, center = TRUE, scale =TRUE)
E <- sweep(E, MARGIN = 2, STATS = sqrt(eval2), FUN = "*")


# Generate data according to model
Y_grid <- Z1 %*% U %*% t(efuns_U) + E %*% t(efuns_E)

# plot a subset of the data
par(mfrow = c(1, 2))
matplot(t(Y_grid)[, 1:5], type = "l")
matplot(t(Y_grid)[, 15:20], type = "l")




# -------------------------------------------------------------------------

# Put data onto basis to make an fd object:
bspl35 <- create.bspline.basis(rangeval = c(0, 100), nbasis = 80, norder = 4)
Y_hip_fd <- Data2fd(argvals = sampling_points, y = as.matrix(t(Y_grid[, 1:101])), basisobj = bspl35)
Y_knee_fd <- Data2fd(argvals = sampling_points, y = as.matrix(t(Y_grid[, 102:202])), basisobj = bspl35)

test <- cov_unstruc_mlfpca_bi_fd(fd_obj_list = list(Y_hip_fd, Y_knee_fd),
                                 id_vec = subject_id)



# Testing U: ---------------------------------------------------------------
test_U_list <- test$K_U_bifd_list
K_U_11_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_U_bifd_list$K_U_11)
K_U_11_true <- efuns_U[1:101, ] %*% diag(eval1) %*% t(efuns_U[1:101, ])
zlims_U_11 <- range(K_U_11_test, K_U_11_true)
filled.contour(K_U_11_test, zlim = zlims_U_11)
filled.contour(K_U_11_true, zlim = zlims_U_11)
#filled.contour(K_U_11_test - K_U_11_true)
mean((K_U_11_test - K_U_11_true)^2)


K_U_12_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_U_bifd_list$K_U_12)
K_U_12_true <- efuns_U[1:101, ] %*% diag(eval1) %*% t(efuns_U[102:202, ])
zlims_U_12 <- range(K_U_12_test, K_U_12_true)
filled.contour(K_U_12_test, zlim = zlims_U_12)
filled.contour(K_U_12_true, zlim = zlims_U_12)
#filled.contour(K_U_11_test - K_U_11_true)
mean((K_U_12_test - K_U_12_true)^2)


K_U_21_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_U_bifd_list$K_U_21)
K_U_21_true <- efuns_U[102:202, ] %*% diag(eval1) %*% t(efuns_U[1:101, ])
zlims_U_21 <- range(K_U_21_test, K_U_21_true)
filled.contour(K_U_21_test, zlim = zlims_U_21)
filled.contour(K_U_21_true, zlim = zlims_U_21)
#filled.contour(K_U_11_test - K_U_11_true)
mean((K_U_21_test - K_U_21_true)^2)


K_U_22_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_U_bifd_list$K_U_22)
K_U_22_true <- efuns_U[102:202, ] %*% diag(eval1) %*% t(efuns_U[102:202, ])
zlims_U_22 <- range(K_U_22_test, K_U_22_true)
filled.contour(K_U_22_test, zlim = zlims_U_22)
filled.contour(K_U_22_true, zlim = zlims_U_22)
#filled.contour(K_U_11_test - K_U_11_true)
mean((K_U_22_test - K_U_22_true)^2)



# Testing E: --------------------------------------------------------------
test_E_list <- test$K_E_bifd_list
K_E_11_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_E_bifd_list$K_E_11)
K_E_11_true <- efuns_E[1:101, ] %*% diag(eval2) %*% t(efuns_E[1:101, ])
zlims_E_11 <- range(K_E_11_test, K_E_11_true)
filled.contour(K_E_11_test, zlim = zlims_E_11)
filled.contour(K_E_11_true, zlim = zlims_E_11)
# filled.contour(K_E_11_test - K_E_11_true)
mean((K_E_11_test - K_E_11_true)^2)


K_E_12_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_E_bifd_list$K_E_12)
K_E_12_true <- efuns_E[1:101, ] %*% diag(eval2) %*% t(efuns_E[102:202, ])
zlims_E_12 <- range(K_E_12_test, K_E_12_true)
filled.contour(K_E_12_test, zlim = zlims_E_12)
filled.contour(K_E_12_true, zlim = zlims_E_12)
#filled.contour(K_E_11_test - K_E_11_true)
mean((K_E_12_test - K_E_12_true)^2)


K_E_21_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_E_bifd_list$K_E_21)
K_E_21_true <- efuns_E[102:202, ] %*% diag(eval2) %*% t(efuns_E[1:101, ])
zlims_E_21 <- range(K_E_21_test, K_E_21_true)
filled.contour(K_E_21_test, zlim = zlims_E_21)
filled.contour(K_E_21_true, zlim = zlims_E_21)
#filled.contour(K_E_11_test - K_E_11_true)
mean((K_E_21_test - K_E_21_true)^2)


K_E_22_test <- eval.bifd(sevalarg = 0:100, tevalarg = 0:100, bifd = test$K_E_bifd_list$K_E_22)
K_E_22_true <- efuns_E[102:202, ] %*% diag(eval2) %*% t(efuns_E[102:202, ])
zlims_E_22 <- range(K_E_22_test, K_E_22_true)
filled.contour(K_E_22_test, zlim = zlims_E_22)
filled.contour(K_E_22_true, zlim = zlims_E_22)
#filled.contour(K_E_11_test - K_E_11_true)
mean((K_E_22_test - K_E_22_true)^2)

