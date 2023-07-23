## Simulation study comparing FUI (implemented in the "lfosr3s" function) and FAMM (implemented in the "pffr" function)

################################################################################
## Load packages
################################################################################
rm(list = ls())
library(refund)
library(dplyr)
source("./lfosrsim.R") ## required package: mvtnorm
source("./lfosr3s.R") ## required packages: lme4, refund, dplyr, mgcv, progress, mvtnorm, parallel


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


################################################################################
## Do simulations on a local laptop (results in the manuscript are from JHPCE)
################################################################################
nsim <- 10
sim_local <- list() ## store simulation results

for(iter in 1:nsim){
  
  ################################################################################
  ## Generate simulated data
  ################################################################################
  set.seed(iter)
  data <- lfosrsim(family, I, J, L, beta_true, psi_true, SNR_B = SNR_B, SNR_sigma = SNR_sigma)
  
  ################################################################################
  ## Implement different estimation methods
  ################################################################################
  ## fit the lfosr3s model
  ptm <- proc.time()
  fit_lfosr3s <- lfosr3s(formula = Y ~ X + (1 | ID), data = data, family = family, 
                         var = TRUE, analytic = TRUE, parallel = FALSE, silent = TRUE)
  lfosr3stime <- (proc.time() - ptm)[3]
  
  ## fit the pffr model
  ### we use "bam" for faster computation, and change k = 15 in bs.yindex to capture curvatures in both scenarios
  ptm <- proc.time()
  fit_pffr <- pffr(formula = Y ~ X + s(ID, bs = "re"), data = data, family = family, algorithm = "bam",
                   bs.yindex = list(bs = "ps", k = 15, m = c(2, 1)))
  pffrtime <- (proc.time() - ptm)[3]
  
  ################################################################################
  ## Organize simulation results
  ################################################################################
  ## lfosr3s
  MISE_lfosr3s <- rowMeans((fit_lfosr3s$betaHat - beta_true)^2)
  cover_lfosr3s <- rep(NA, nrow(beta_true))
  cover_lfosr3s_pw <- matrix(FALSE, nrow(beta_true), L) 
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(fit_lfosr3s$betaHat[p,]+fit_lfosr3s$qn[p]*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) > beta_true[p,])
    cover_lower <- which(fit_lfosr3s$betaHat[p,]-fit_lfosr3s$qn[p]*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) < beta_true[p,])
    cover_lfosr3s[p] <- length(intersect(cover_lower, cover_upper)) == L
    cover_upper_pw <- which(fit_lfosr3s$betaHat[p,]+1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) > beta_true[p,])
    cover_lower_pw <- which(fit_lfosr3s$betaHat[p,]-1.96*sqrt(diag(fit_lfosr3s$betaHat.var[,,p])) < beta_true[p,])
    cover_lfosr3s_pw[p,intersect(cover_lower_pw, cover_upper_pw)] <- TRUE
  }
  
  ## pffr
  coef_pffr <- coef(fit_pffr, n1 = L)
  betaHat_pffr <- betaHat_pffr.var <- beta_true
  betaHat_pffr[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$value)
  betaHat_pffr[2,] <- as.vector(coef_pffr$smterms$`X(yindex)`$value)
  betaHat_pffr.var[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$se^2)
  betaHat_pffr.var[2,] <- as.vector(coef_pffr$smterms$`X(yindex)`$se^2)
  MISE_pffr <- rowMeans((betaHat_pffr - beta_true)^2)
  cover_pffr <- matrix(FALSE, nrow(beta_true), L)
  for(p in 1:length(cover_lfosr3s)){
    cover_upper <- which(betaHat_pffr[p,]+1.96*sqrt(betaHat_pffr.var[p,]) > beta_true[p,])
    cover_lower <- which(betaHat_pffr[p,]-1.96*sqrt(betaHat_pffr.var[p,]) < beta_true[p,])
    cover_pffr[p,intersect(cover_lower, cover_upper)] <- TRUE
  }
  ## organize pffr results
  fit_pffr_cleaned <- list(betaHat = betaHat_pffr, betaHat.var = betaHat_pffr.var)
  
  ## results of a single simulation
  sim_result <- list(MISE_lfosr3s = MISE_lfosr3s, MISE_pffr = MISE_pffr,
                     cover_lfosr3s = cover_lfosr3s, cover_lfosr3s_pw = cover_lfosr3s_pw,
                     cover_pffr = cover_pffr, time_lfosr3s = lfosr3stime, time_pffr = pffrtime,
                     I = I, J = J, L = L, SNR_B = SNR_B, SNR_sigma = SNR_sigma,
                     family = family, scenario = scenario)
  
  sim_local[[iter]] <- sim_result
  print(iter)
  
}


################################################################################
## Obtain MISE, coverage, and computing time
################################################################################
## MISE
MISE_lfosr3s <- lapply(sim_local, '[[', 1) %>% bind_rows()
MISE_pffr <- lapply(sim_local, '[[', 2) %>% bind_rows()
colMeans(MISE_lfosr3s)
colMeans(MISE_pffr)


## coverage
### joint
cover_lfosr3s <- t(lapply(sim_local, '[[', 3) %>% bind_cols())
rownames(cover_lfosr3s) <- 1:nsim
### pointwise
cover_lfosr3s_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_lfosr3s_pw[i,,1] <- sim_local[[i]]$cover_lfosr3s_pw[1,]
  cover_lfosr3s_pw[i,,2] <- sim_local[[i]]$cover_lfosr3s_pw[2,]
}
cover_pffr_pw <- array(NA, dim = c(nsim, L, p))
for(i in 1:nsim){
  cover_pffr_pw[i,,1] <- sim_local[[i]]$cover_pffr[1,]
  cover_pffr_pw[i,,2] <- sim_local[[i]]$cover_pffr[2,]
}
### print results
print(paste0("lfosr3s coverage (joint CI): ", colMeans(cover_lfosr3s)[2]))
print(paste0("lfosr3s coverage (pointwise CI): ", apply(cover_lfosr3s_pw, 3, mean)[2]))
print(paste0("pffr coverage (pointwise CI): ", apply(cover_pffr_pw, 3, mean)[2]))


## time
time_lfosr3s <- lapply(sim_local, '[[', 6) %>% bind_rows()
print(apply(time_lfosr3s, 2, median))
time_pffr <- lapply(sim_local, '[[', 7) %>% bind_rows()
print(apply(time_pffr, 2, median))
