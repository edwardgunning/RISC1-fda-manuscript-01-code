# ------------------------------------------------------------------------#
# Generate data for simulation scenario.
# Generate by a subset model fit to our data.
# This script gets the parameters for the model.
# (We keep it very simple)
# ------------------------------------------------------------------------#

# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(fda)        # CRAN v5.5.1
library(lme4)       # CRAN v1.1-28

# Load some helper functions for functional data manipulation: ------------
source(file = here::here("code","functions","functions-helper-smoothing.R"))
source(file = here::here("code","functions","BFMM-paper-helper-functions.R"))
# And helper function for projecting mean function back onto the FPCs:
source(here::here("code","functions","function-project-mean-onto-fpcs.R"))

# Paths -------------------------------------------------------------------
data_path <- here::here("data")
results_path <- here::here("outputs", "results")

# Read in data ------------------------------------------------------------
# Simplified dataset stored as an rds object.
# functional data is represented in terms of 80 basis coefficients:
subject_side_coef_reg <- readRDS(file.path(data_path, "subject_side_coef_reg.rds"))

# set up basis to represent the data.
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)


# Do mfpca ----------------------------------------------------------------
# To do this we need to put the coefficcients into a 3-d array
# of 
# N_curve x N_coef x N_dimensions
# (i.e., hip and knee coefs are in slices of array)
# First, lets split them into separate FDA objects and do exploratory plots:
subject_side_coef_hip <- subject_side_coef_reg[location== "Hip"]
subject_side_coef_knee <- subject_side_coef_reg[location== "Knee"]
stopifnot(subject_side_coef_hip$subject_id == subject_side_coef_knee$subject_id)
stopifnot(subject_side_coef_hip$side == subject_side_coef_knee$side)


subject_side_fd_knee <- fd(
  coef = coef_to_mat(subject_side_coef_knee[, paste0("lm_coef_", 1:80)]),
  basisobj = bspl80)

subject_side_fd_hip <- fd(
  coef = coef_to_mat(subject_side_coef_hip[, paste0("lm_coef_", 1:80)]),
  basisobj = bspl80)

par(mfrow = c(1, 2))
plot(subject_side_fd_hip)
title("Hip")
plot(subject_side_fd_knee)
title("Knee")
dev.off()


# use `fda` package which requires a multivariate fd object
# so make array of coefs
coef_array <- array(data = NA, dim = c(dim(subject_side_fd_hip$coefs), 2))
coef_array[,, 1] <- subject_side_fd_hip$coefs
coef_array[,, 2] <- subject_side_fd_knee$coefs
bfd_obj <- fd(coef = coef_array, basisobj = bspl80)

# Plot to check:
par(mfrow = c(1, 2))
plot(bfd_obj[, 1])
title("Hip")
plot(bfd_obj[, 2])
title("Knee")
dev.off()

# Set up functional parameter object for the FPCA:
# - use same basis 80 bsplines
# - use no smoothing (lambda = 0)
harm_fd_par <- fdPar(fdobj = bspl80, Lfdobj = int2Lfd(2), lambda = 0)
bfpca <- pca.fd(fdobj = bfd_obj, nharm = 50, harmfdPar = harm_fd_par)


# We want to retain enough fpcas s.t. > 99 of var is retained to be used in simulation
k_retain <- min(which(cumsum(bfpca$varprop) > 0.99))
colours <- c(rep("darkgreen", k_retain), rep("red4", 50 - k_retain))

# Re-run bfpca onlt retianing the specified number of fpcs.
bfpca <- pca.fd(fdobj = bfd_obj,
                nharm = k_retain,
                harmfdPar = harm_fd_par)

# The fpc scores are seperated into parts projected onto the
# hip and knee dimensions in each slice of the array.
# To get overall score, we need to sum over slices.
scores <- apply(bfpca$scores,
                MARGIN = c(1, 2),
                FUN = sum)

mean_eval <- eval.fd(evalarg = 0:100, mean.fd(bfd_obj))
mean_eval_vec <- c(mean_eval[,1,1], mean_eval[,1,2])


# Evaluate basis functions on a grid:
phi <- eval.fd(evalarg = 0:100,
               fdobj = bfpca$harmonics)
Phi <- rbind(phi[,,1], phi[,,2])

scores <- apply(bfpca$scores,
                        MARGIN = c(1, 2),
                        FUN = sum)



# Model for scores -----------------------------------------------
# Prepare data for fitting linear mixed effects (lme) models:
# First put scores and covariates in to a dataset:
covariates_dt <- subject_side_coef_hip[, .(sex,
                                           speed_cent,
                                           subject_id = factor(subject_id))]

colnames(scores) <- paste0("score_", seq_len(ncol(scores)))
lme_dt <- cbind(covariates_dt, scores)



# Fixed Effects Formula for the Model: ------------------------------------
model_formula_fixef <- "sex + speed_cent"

lme_fit <- fit_lme_to_scores(
  lme_df = lme_dt,
  k_retain = k_retain,
  fixef_formula = model_formula_fixef, 
  ranef_formula =  "+ (1|subject_id)",
  REML = TRUE)

# Extract Fixed Effects:
B <- lme_fit$fixef_mat

# Extract Random-Intercept Variances:
q_vec <- lme_fit$q_vec

# Extract Error Variances:
s_vec <- lme_fit$s_vec

# For siulation data, only take first ten inds to keep simulation start
# (hence subsetting [, 1:10])
sim_data_list <- list(
  k = k_retain,
  mean_fd = bfpca$meanfd,
  mean_eval_vec = mean_eval_vec,
  B_empirical_scores = B, # don't take intercept (2:3), work with mean-centered data
  Phi = Phi,
  arg_vals = 0:100,
  q_vec = q_vec,
  s_vec = s_vec,
  speed_sd = sd(unique(covariates_dt$speed_cent)) # standard deviation of speed variable used to generate
)

saveRDS(object = sim_data_list,
        file = file.path(results_path, "BFMM-simulation-parameters.rds"))
