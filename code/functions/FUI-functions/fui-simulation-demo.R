################################################################################
## Load packages
################################################################################
rm(list = ls())
library(refund) # CRAN v0.1-26 
library(dplyr)  # CRAN v1.0.8
source("code/functions/FUI-functions/lfosr3s.R") # put your own paths to FUI functions here
source("code/functions/FUI-functions/lfosr3s-updated.R")
source("code/functions/FUI-functions/lfosrsim.R")

################################################################################
## Set parameters
################################################################################
I <- 50 ## number of subjects
J <- 5 ## mean number of observations per subject
L <- 50 ## dimension of the functional domain
SNR_B <- 0.5 ## relative importance of random effects
SNR_sigma <- 1 ## signal-to-noise ratio
family <- "gaussian" ## family of longitudinal functional data

## specify true fixed and random effects functions
grid <- seq(0, 1, length = L)
beta_true <- matrix(NA, 2, L)
scenario <- 1 ## indicate the scenario to use
if(scenario == 1){
  beta_true[1,] = -0.15 - 0.1*sin(2*grid*pi) - 0.1*cos(2*grid*pi)
  beta_true[2,] = dnorm(grid, .6, .15)/20
}else if(scenario == 2){
  beta_true[1,] <- 0.53 + 0.06*sin(3*grid*pi) - 0.03*cos(6.5*grid*pi)
  beta_true[2,] <- dnorm(grid, 0.2, .1)/60 + dnorm(grid, 0.35, .1)/200 - 
    dnorm(grid, 0.65, .06)/250 + dnorm(grid, 1, .07)/60
}
rownames(beta_true) <- c("Intercept", "x")

psi_true <- matrix(NA, 2, L)
psi_true[1,] <- (1.5 - sin(2*grid*pi) - cos(2*grid*pi) )
psi_true[1,] <- psi_true[1,] / sqrt(sum(psi_true[1,]^2))
psi_true[2,] <- sin(4*grid*pi)
psi_true[2,] <- psi_true[2,] / sqrt(sum(psi_true[2,]^2))


# First let's do a single simulated dataset: -----------------------------
set.seed(2023)
data <- lfosrsim(family, I, J, L, beta_true, psi_true, SNR_B = SNR_B, SNR_sigma = SNR_sigma)
################################################################################
## Implement different estimation methods
################################################################################
## fit the lfosr3s model
fit_lfosr3s_standard <- lfosr3s(formula = Y ~ X + (1 | ID), data = data, family = family, 
                       var = TRUE, analytic = FALSE, parallel = TRUE)


fit_lfosr3s_updated <- lfosr3s_updated(formula = Y ~ X + (1 | ID), data = data, family = family, 
                       var = TRUE, analytic = FALSE, parallel = TRUE)
par(mfrow = c(1, 2))
for(i in 1:2) {
  se_standard <- sqrt(diag(fit_lfosr3s_standard$betaHat.var[,,i]))
  est_standard <- fit_lfosr3s_standard$betaHat[i,]
  se_updated <- sqrt(diag(fit_lfosr3s_updated$betaHat.var[,,i]))
  est_updated <- fit_lfosr3s_updated$betaHat[i,]
  
  plot(est_standard, type = "l", 
       ylim = range(est_standard + 2 * se_standard,
                    est_standard - 2 * se_standard,
                    est_updated + 2 * se_updated,
                    est_updated - 2 * se_updated))
  lines(est_standard - 2 * se_standard, lty = 2)
  lines(est_standard + 2 * se_standard, lty = 2)
  title(paste("mean SE Ratio =", mean(se_updated / se_standard)))
  lines(est_updated, col = 2)
  lines(est_updated - 2 * se_updated, lty = 2, col = 2)
  lines(est_updated + 2 * se_updated, lty = 2, col = 2)
  abline(h = 0)
}


# Now let's do it repeatedly ----------------------------------------------
N_sim <- 10 # small no. of iterations because doing this on laptop.
SE_ratio <- matrix(NA, nrow = N_sim, ncol = 2)
for(i in seq_len(N_sim)) {
  print(paste0("Iteration", i))
  data <- lfosrsim(family, I, J, L, beta_true, psi_true, SNR_B = SNR_B, SNR_sigma = SNR_sigma)
  ################################################################################
  ## Implement different estimation methods
  ################################################################################
  ## fit the lfosr3s model
  fit_lfosr3s_standard <- lfosr3s(formula = Y ~ X + (1 | ID),
                                  data = data,
                                  family = family, 
                                  var = TRUE, 
                                  analytic = FALSE, 
                                  parallel = TRUE)
  

  fit_lfosr3s_updated <- lfosr3s_updated(formula = Y ~ X + (1 | ID),
                                         data = data,
                                         family = family, 
                                         var = TRUE,
                                         analytic = FALSE,
                                         parallel = TRUE)
  # take average SE ratio between two approaches averaged over the functional domain
  # ratio = updated / standard so expect top be > 1 
  for(j in 1:2) {
    SE_ratio[i, j] <- mean(sqrt(diag(fit_lfosr3s_updated$betaHat.var[,,j])) / sqrt(diag(fit_lfosr3s_standard$betaHat.var[,,j])))
  }
  
}
par(mfrow = c(1, 1))
boxplot(SE_ratio)
