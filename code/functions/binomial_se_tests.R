# -------------------------------------------------------------------------
# Simple tests to make sure binomial SE function works as expected with
# vectorised inputs.
# -------------------------------------------------------------------------

# Load function: ----------------------------------------------------------
functions_path <- here::here("code", "functions")
source(file.path(functions_path, "binomial_se.R"))

# Test 1 - Vectorised p and n: --------------------------------------------
phat_vec <- c(0.95, 0.95, 0.9)
n_vec <- c(475, 1900, 900)
expected_results_vec <- c(0.01, 0.005, 0.01)
# test:
testthat::expect_equal(object = binomial_se(phat = phat_vec, n = n_vec),
                       expected_results_vec)



# Test 2: Vector p scalar n -----------------------------------------------
# (this is what used for simulation results)
phat_vec <- c(0.975, 0.95, 0.9)
n_fixed <- 500
expected_results_vec <- c( # (to 8 decimal places)
  0.00698212,
  0.00974679,
  0.01341641)

testthat::expect_equal(object = binomial_se(phat = phat_vec, n = n_fixed),
                       expected_results_vec,
                       tolerance = 10^-8)



