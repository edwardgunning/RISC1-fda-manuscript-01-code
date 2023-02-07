# -------------------------------------------------------------------------
# Test Multilevel FPCA Covariance Estimation for UNIVARIATE functional data:
# -------------------------------------------------------------------------
source(here::here("code", #path to functions script here 
                  "functions", 
                  "functions-unstructured-covariance.R"))
library(Matrix)    # CRAN v1.4-0
library(funData)   # CRAN v1.3-8
library(fda)       # CRAN v5.5.1

# ------------------------------------------------------------------------#
# Generate data from a univariate functional mixed model
# ------------------------------------------------------------------------#

N_tot <- 5000 # no. of subjects
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
Phi <- eval.fd(sampling_points, funData::funData2fd(funData::eFun(argvals = sampling_points, M = 3, ignoreDeg = c(1,2), type = "PolyHigh")))
Psi <- eval.fd(sampling_points, funData::funData2fd(funData::eFun(argvals = sampling_points, M = 3, type = "Fourier")))

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
Y_grid <- Z1 %*% U %*% t(Phi) + E %*% t(Psi)

# plot a subset of the data
par(mfrow = c(1, 2))
matplot(t(Y_grid)[, 1:5], type = "l")
matplot(t(Y_grid)[, 15:20], type = "l")

# Put data onto basis to make an fd object:
bspl35 <- create.bspline.basis(rangeval = c(0, 100), nbasis = 35, norder = 5)
Y_fd <- Data2fd(argvals = sampling_points, y = as.matrix(t(Y_grid)), basisobj = bspl35)
Y_coef <- t(Y_fd$coefs)



# Test our function -------------------------------------------------------
test <- cov_unstruc_mlfpca_uni_fd(fd_obj = Y_fd, id_vec = subject_id)
str(test)

# ------------------------------------------------------------------------#
# Compare results (visually)
# ------------------------------------------------------------------------#

K_U_bifd <- test$K_U_bifd
K_U_est <- eval.bifd(sampling_points, sampling_points, K_U_bifd)
K_U_true <- Phi %*% diag(eval1) %*% t(Phi)
zlims_U <- range(K_U_est, K_U_true)
filled.contour(K_U_est, zlim = zlims_U)
filled.contour(K_U_true, zlim = zlims_U)

K_E_bifd <- test$K_E_bifd
K_E_est <- eval.bifd(sampling_points, sampling_points, K_E_bifd)
K_E_true <- Psi %*% diag(eval2) %*% t(Psi)
zlims_E <- range(K_E_est, K_E_true)
filled.contour(K_E_est, zlim = zlims_E)
filled.contour(K_E_true, zlim = zlims_E)

# 3-d plots
par(mfrow = c(1, 2))
persp(K_U_est, zlim = zlims_U)
persp(K_U_true, zlim = zlims_U)

persp(K_E_est, zlim = zlims_E)
persp(K_E_true, zlim = zlims_E)

# RMSE:
mean((K_U_est - K_U_true)^2)
mean((K_E_est - K_E_true)^2)
# all small
