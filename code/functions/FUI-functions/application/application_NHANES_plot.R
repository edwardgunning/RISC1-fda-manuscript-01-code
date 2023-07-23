## This script makes Figure 6 of the paper using the bootstrap results obtained from the cluster.
## The NHANES application code is modified and fully parallelized on the cluster. Therefore, the loaded bootstrap 
## results below have different format from that saved in the "application_NHANES_fit.R" as described.
## The bootstrap results are too large to share in the supplementary material.

rm(list=ls())
## Check for (and install if necessary) packages needed to create figures
pckgs <- c("tidyverse", "stringr", "mgcv", "mvtnorm", "gridExtra")
sapply(pckgs, function(x) if(!require(x,character.only=TRUE,quietly=TRUE)) {
        install.packages(x)
        require(x, character.only=TRUE)
})
rm(list=c("pckgs"))

## Directory for where bootstrap results are saved
data_path <- "./results/"
## set the directory for saving figures
figure_path <- "./figures"

## Load the bootstrap results
## the structure of "df" is an array where the first dimension is the bootstrap sample
## and the 2nd/3rd dimensions are the pth coefficient and functional domain (minute of the day), respectively. 
## Each entry is a pointwise coefficient estimate.
## Also, the first index of the array (i.e. df[1,,]) are the estimated coefficients from the original data.
df <- read_rds(file.path(data_path, "merged_bootstrap_NHANES_small.rds"))

## use the first 100 bootstrap samples which we have full data for inference
nboot    <- 100
inx_comp <- which(apply(df, 1, function(x) !any(is.na(x))))
P        <- dim(df)[2]
df_sub   <- df[inx_comp[1:(nboot+1)],,]

## smooth the bootstrapped coefficients using 15 cubic regresison spline basis functions (therefore k=16)
## and REML smoothing parameter selection
df_sub_sm <- array(NA, dim=dim(df_sub))
tind <- seq(0,1,len=1440)
for(i in 1:(nboot+1)){
        for(p in 1:P){
                y_ip <- as.vector(df_sub[i,p,])
                fit_i_p <- gam(y_ip ~ s(tind, bs="cr", k=16), method="REML")                
                df_sub_sm[i,p,] <- fit_i_p$fitted.values
        }
        if(i %% 10 == 0) print(i)
}


## Get pointwise CIs
## get estimated coefficients from the original data
beta_hat <- df_sub_sm[1,,]

## set inflation paramter \phi for the pointwise CIs, use \phi = 1.1
phi_pw <- 1.1

## get SD over the bootstrap sample -- IGNORE first index df_sub_sm[1,,] 
## because these are the in-sample estimates
SD_pointwise <- apply(df_sub_sm[-1,,], c(2,3), sd, na.rm=TRUE)
CI_LB_pointwise <- df_sub_sm[1,,] - 2*phi_pw*SD_pointwise
CI_UB_pointwise <- df_sub_sm[1,,] + 2*phi_pw*SD_pointwise


## set the random seed for obtaining joint CIs
set.seed(10101)
## get joint CIs using "nbs" samples
nbs <- 1000
xi_arr <- array(NA, dim=c(nbs, P, 1440))
u_mat <- matrix(NA, nbs, P)
for(p in 1:P){  
        ## get variance/sd estimates again ignoring the first index 
        sigma_p <- var(df_sub_sm[-1,p,], na.rm = T)
        sd_p    <- sqrt(diag(sigma_p))
        ## sample coefficients from multivariate normal
        xi_arr[,p,] <- rmvnorm(n=nbs,mean=beta_hat[p,], sigma=sigma_p)
        ## get the Z-scores along the functional domain
        u_mat[,p] <- apply(xi_arr[,p,], 1, function(x) max( abs(x-beta_hat[p,]) / sd_p ))
        print(paste0("Variable = ", p))
}
## set alpha level
alpha <- 0.95
## obtain u_n for each coefficient
q_vec <- apply(u_mat, 2, quantile, alpha)
## set inflation paramter \phi for the joint CIs, again use \phi = 1.1
phi_jt <- 1.1
## set up empty matrices for holding upper and lower bounds
CI_LB_joint <- CI_UB_joint <- matrix(NA, nrow(CI_LB_pointwise), ncol(CI_LB_pointwise))
## loop over coefficients and get the CIs
for(p in 1:P){
        CI_LB_joint[p,] <- df_sub_sm[1,p,] - phi_jt*q_vec[p]*SD_pointwise[p,]
        CI_UB_joint[p,] <- df_sub_sm[1,p,] + phi_jt*q_vec[p]*SD_pointwise[p,]
}

## set the ticks and labels for the x-axis
xlab_inx <- (c(1,6,12,18,23)*60 + 1)/(1440)
xlab     <- c("01:00","06:00","12:00","18:00","23:00")
## create vector of labels for each coefficient
## note that the order of the coefficients in the model fit is NOT the plotting order
labels <- c("Intercept","Gender","Age",
            "Day-of-week: Monday","Day-of-week: Tuesday","Day-of-week: Wednesday","Day-of-week: Thursday","Day-of-week: Friday","Day-of-week: Saturday",
            "Day")
## set the indices for custom plotting order
inx_plt <- c(1, 10, 2, 3, 4:9)
## set global text size
textsize <- 2

## visualize the result using ggplot2
plist <- list()
## loop over coefficients to plot
for(inx in seq_along(inx_plt)){
        ## get current index
        i <- inx_plt[inx]
        ## set the y-axis limits
        ## we use custom limits for "age" (i = 3) and "day" (i = 10) due to scaling
        ## as well as baseline (i = 1)
        ylims <- c(-2.3,3.2)
        if(i %in% c(3,10)){
                ylims <- c(-0.2,0.3)
        }
        if(i %in% 1){
                ylims <- c(-11,1)
        } 
        ## add in multiplicative constant to multiply estimates and CIs by for the "age" and "day" coefficients
        ## c = 1/sd(variable) from the data
        ## copied and pasted here to avoid re-processing the data in this script
        c <- 1
        if(i %in% 3){
                c <- (1/3.758175)
        } 
        if(i %in% 10){
                c <- (1/1.95203)
        }
        
        ## organize the results into a plot format
        beta.hat.plt <- data.frame(s = tind, beta = beta_hat[i,]*c,
                                   lower = CI_LB_pointwise[i,]*c,
                                   upper = CI_UB_pointwise[i,]*c,
                                   lower.joint = CI_LB_joint[i,]*c,
                                   upper.joint = CI_UB_joint[i,]*c)
        
        p.CI <- ggplot() +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
                geom_ribbon(aes(x = s, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray30", alpha = 0.2) +
                geom_ribbon(aes(x = s, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
                geom_line(aes(x = s, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
                scale_colour_manual(name="", values=c("Estimate"="blue3")) +
                labs(title = labels[i]) +
                geom_hline(yintercept=0, color = "black", lwd=textsize/4, linetype = "dotted") +
                ylim(ylims) +
                scale_x_continuous(name = "Time of Day (s)", breaks = xlab_inx, labels = xlab)
        
        if(inx == 1){ ## add legend for the first panel
                p.CI <- p.CI + labs(y = bquote(beta[.(inx-1)](s))) +
                        theme(legend.title=element_blank(),
                              legend.position = c(0.02, 1.02),
                              legend.justification = c("left", "top"),
                              legend.box.just = "right",
                              legend.margin = margin(0, 0, 0, 0),
                              legend.background = element_rect(fill=alpha('white', 0)))
        }
        if(inx %in% c(2:4)){ ## add corresponding ylabs
                p.CI <- p.CI + labs(y = bquote(beta[.(inx-1)](s))) +
                        theme(legend.position = "none")
        } 
        if(inx %in% 5:10){ ## add corresponding ylabs
                p.CI <- p.CI + labs(y = bquote(f[.(inx-4)](s))) +
                        theme(legend.position = "none")
        }
        
        plist[[inx]] <- p.CI
}
## make the plot
do.call("grid.arrange", c(plist, nrow = 2))

