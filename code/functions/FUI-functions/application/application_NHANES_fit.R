# Implement FUI in the NHANES study.
# For fast implementation, parallelization is necessary

rm(list=ls())
## Check for packages needed to run analyses/install the rnhanesdata package.
pckgs <- c("devtools","mgcv","lme4","glmmTMB","tidyverse","parallel","doParallel","foreach")
sapply(pckgs, function(x) if(!require(x,character.only=TRUE,quietly=TRUE)) {
  install.packages(x)
  require(x, character.only=TRUE)
})
rm(list=c("pckgs"))
## Install the rnhanesdata package and dependencies.
## This may take a few minutes because of the size of the data package.
if(!require("rnhanesdata")){
  install_github("andrew-leroux/rnhanesdata",build_vignettes = FALSE)
  require("rnhanesdata")
}
## load activity count and wear/non-wear flag data
data("PAXINTEN_C"); data("PAXINTEN_D")
data("Flags_C"); data("Flags_D")
data("Covariate_C"); data("Covariate_D")

## re-code activity counts which are considered "non-wear" to be 0
## this doesn't impact much data, most estimated non-wear times correspond to 0 counts anyway
PAXINTEN_C[,paste0("MIN",1:1440)] <- PAXINTEN_C[,paste0("MIN",1:1440)]*Flags_C[,paste0("MIN",1:1440)]
PAXINTEN_D[,paste0("MIN",1:1440)] <- PAXINTEN_D[,paste0("MIN",1:1440)]*Flags_D[,paste0("MIN",1:1440)]


## mege 2003-2004 and 2005-2006 waves' data
PAXINTEN <- bind_rows(PAXINTEN_C, PAXINTEN_D)
Flags    <- bind_rows(Flags_C, Flags_D)
Covariate <- bind_rows(Covariate_C, Covariate_D) %>% 
  mutate(Age = RIDAGEEX/12) %>% 
  dplyr::select(SEQN, Gender, BMI_cat, Age)

## clear up the workspace
rm(PAXINTEN_C, PAXINTEN_D, Flags_C, Flags_D, Covariate_C, Covariate_D)

## 1) subset to good days of data (>= 10 hours estiamted wear time) and good quality data indicators (PAXCAL/PAXSTAT)
## 2) create day number of longitudinal information
## 3) re-code day of the week as factor variable (for ease of use in regressions)
PAXINTEN$nmins_wear <- rowSums(Flags[,paste0("MIN",1:1440)], na.rm=TRUE)
PAXINTEN <- 
  PAXINTEN %>% 
  mutate(good_day = as.numeric(nmins_wear >= 600 & PAXCAL %in% 1 & PAXSTAT %in% 1)) %>% 
  group_by(SEQN) %>% 
  mutate(day = 1:n(),
         n_good_days = sum(good_day)) %>% 
  ungroup() %>% 
  filter(good_day == 1, n_good_days >= 3) %>% 
  dplyr::select(-PAXCAL, -PAXSTAT) %>% 
  mutate(dow_fac = factor(WEEKDAY, levels=1:7, labels=c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")))



## merge PA and covariate data
df_fit <- left_join(PAXINTEN, Covariate, by="SEQN")
## add PA to the data frame as a matrix
df_fit$X <- I(as.matrix(df_fit[,paste0("MIN",1:1440)]))
## replace any NAs with 0
df_fit$X[is.na(df_fit$X)] <- 0
## get rid of redundant columns
df_fit   <- 
  df_fit %>% 
  dplyr::select(-one_of(paste0("MIN",1:1440))) 
## create active/inactive indicators
df_fit$Y <- I(matrix(as.numeric(df_fit$X >= 100),ncol=1440,nrow=nrow(df_fit$X),byrow=FALSE))
## clear up the workspace
rm(PAXINTEN, Flags, Covariate)
gc()


### bin the active/inactive indicators
# binning window
tlen <- 1
nt   <- floor(1440/tlen)
tind <- seq(0,1,len=nt)
# create a list of indices for binning into tlen minute windows
inx_col_ls       <- split(1:1440, rep(1:nt,each=tlen))
df_fit$Y_bin     <- I(sapply(inx_col_ls, function(x) rowSums(df_fit$Y[,x,drop=FALSE])))
df_fit$X_bin     <- I(sapply(inx_col_ls, function(x) rowSums(df_fit$X[,x,drop=FALSE])))
# subset to just ages 18-30 to avoid issues associated with potential non-linear effects across age 
df_fit <- 
  df_fit %>% 
  filter(Age < 30 & Age > 18 & !is.na(BMI_cat)) 

## functional domain 
sind <- 1:ncol(df_fit$Y_bin)

## create empty matrices/arrays for holding results
uid <- unique(df_fit$SEQN)
nid <- length(uid)
coef_fixed_mat_BB  <- coef_fixed_mat_B  <- matrix(NA, ncol=10, nrow=length(sind))     # fixed effects
coef_random_mat_BB <- coef_random_mat_B <- array(NA, dim=c(length(sind), nid, 2))    # random effects (here only 2 per subject b/c we fit a random intercept + slope model)
eta_mat_BB <- eta_mat_B <-
  p_mat_BB <- p_mat_B <-matrix(NA, nrow= nrow(df_fit), ncol=length(sind)) # linear predictors

## scale and center continuous variables (helps with convergence issues)
df_fit$Age_sc <- scale(df_fit$Age)
df_fit$day_sc <- scale(df_fit$day)



## fit the models
time_st <- Sys.time()
expit <- function(x) 1/(1+exp(-x))

## set up cluster for parallel computing
cl <- makeCluster(14,setup_strategy = "sequential")
registerDoParallel(cl)

time_start <- Sys.time()
results <- foreach(s = sind, .packages=c("glmmTMB", "tidyverse")) %dopar% {
  df_fit_s <- df_fit %>% dplyr::select(-X,-Y,-Y_bin) %>% data.frame("Y" = df_fit$Y_bin[,s])
  
  df_fit_s$SEQN_fac <- as.factor(df_fit_s$SEQN)
  df_fit_s$nfail <- tlen - df_fit_s$Y
  
  time_st <- Sys.time()
  fit_B <- try(glmmTMB::glmmTMB(cbind(Y, nfail) ~ Gender + Age_sc + dow_fac +  (day_sc|SEQN),
                       data=df_fit_s, family=binomial))
  time_ed <- Sys.time()
  time_ed-time_st

  if(!inherits(fit_B, "try-error")){
    eta_hat_B <- predict(fit_B, type="link")
    p_hat_B   <- predict(fit_B, type="response")
    
    ret <- list(eta_hat = eta_hat_B,
                p_hat = p_hat_B,
                coef_fixed = fit_B$fit$par[1:10],
                coef_random = as.matrix(ranef(fit_B)$cond$SEQN))
    return(ret)
  } else {
    return(NULL)
  }
}
time_end <- Sys.time()
time_end-time_start
stopCluster(cl)

## Save the result into an RDS file. This is just the point estimation result.
## The bootstrap results loaded in "application_NHANES_make_figures.R" have different formats due to parallelization.
write_rds(list("results" = results, "time_st" = time_start, time_end = time_end), file=outfile)

## derive point estimators from step 1
beta_tilde <- t(vapply(results, function(x) as.vector(x$coef_fixed), numeric(10)))
beta_tilde_sc <- beta_hat
beta_tilde_sc[,3] <- beta_hat_sc[,3]/sd(df_fit$Age) ## scale back age
beta_tilde_sc[,10] <- beta_hat_sc[,10]/sd(df_fit$day) ## scale back day

## derive smooth estimators from step 2
beta_hat <- apply(beta_tilde_sc, 2, function(x) gam(x ~ s(tind, bs="cr", k = 15))$fitted.values)
colnames(beta_hat) <- c("Intercept", "Gender (Female)", "Age", 
                        "DoW:Mon", "DoW:Tue", "DoW:Wed", "DoW:Thu", "DoW:Fri", "DoW:Sat", "Day")

