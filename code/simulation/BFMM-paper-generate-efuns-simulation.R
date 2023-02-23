# ------------------------------------------------------------------------#
# Generate eigenfunctions for second scenario in the simulation.
# ------------------------------------------------------------------------#

# ------------------------------------------------------------------------#
# Use "splitting" method from Happ and Greven to generate (different)
# bivariate bases for the level of the random intercept and random error,
# respectively
# ------------------------------------------------------------------------#

# ------------------------------------------------------------------------#
# Explanation:
# type = "split": The basis functions of an underlying 'big' orthonormal 
# basis are split in M parts, translated and possibly reflected. This yields
# an orthonormal basis of multivariate functions with M elements. This option
# is implemented only for one-dimensional domains.
# ------------------------------------------------------------------------#

# ------------------------------------------------------------------------#
# Reference
# Happ, C., & Greven, S. (2018). Multivariate Functional Principal Component
# Analysis for Data Observed on Different (Dimensional) Domains. Journal of
# the American Statistical Association, 113(522), 649–659.
# https://doi.org/10.1080/01621459.2016.1273115
# ------------------------------------------------------------------------#

results_path <- here::here("outputs", "results")
k <- readRDS(file.path(results_path, "BFMM-simulation-parameters.rds"))$k
set.seed(2) # because sign in front of eignefunctions in each dimensoon can flip

# Basis functions for U ---------------------------------------------------
sim_bfundata_U <- funData::simMultiFunData(
  argvals = replicate(n = 2, expr = 0:100, simplify = FALSE), # domains of bivariate functional data.
  eFunType = "Fourier", # Fourier basis for level 1
  M = k, # Number of basis functions 
  N = 1, # NA because only taking basis functions
  eValType = "linear", # NA because only taking basis functions
  type = "split") # ('splitting' method, see Happ and Greven)

# And extract:
efuns_fourier <- sim_bfundata_U$trueFuns
efuns_U_hip <- fda::eval.fd(0:100, funData::funData2fd(efuns_fourier[[1]]))
efuns_U_knee <- fda::eval.fd(0:100, funData::funData2fd(efuns_fourier[[2]]))
efuns_U <-  -1 * rbind(efuns_U_hip, efuns_U_knee) # so consistent with initial sim


# Basis functions for E: --------------------------------------------------
set.seed(1996)
sim_bfundata_E <- funData::simMultiFunData(argvals = replicate(n = 2, expr = 0:100, simplify = FALSE),
                                           eFunType = "Poly", # Orthogonal polynomial basis for level 2
                                           M = k,
                                           N = 1,
                                           eValType = "linear",
                                           type = "split")
efuns_poly <- sim_bfundata_E$trueFuns
efuns_E_hip <- fda::eval.fd(0:100, funData::funData2fd(efuns_poly[[1]]))
efuns_E_knee <- fda::eval.fd(0:100, funData::funData2fd(efuns_poly[[2]]))
efuns_E <- rbind(efuns_E_hip, efuns_E_knee)


saveRDS(object = list(efuns_U = efuns_U, efuns_E = efuns_E),
        file = file.path(results_path, "BFMM-simulation-efuns.rds"))


