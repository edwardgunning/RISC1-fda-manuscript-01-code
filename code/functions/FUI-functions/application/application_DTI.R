## Implement FUI in the DTI study.
## The dataset "DTI_full" is not publicly available since it includes sensitive age and scan date information

################################################################################
## load packages
################################################################################
rm(list = ls())
library(refund)
source("./lfosr3s.R")
library(ggplot2)
library(gridExtra)


################################################################################
## load dataset
################################################################################
data("DTI") ## a puclicly available dataset with limited variables
load("./DTI_full.RDA") ## a restricted dataset


################################################################################
## do FUI
################################################################################
fit_DTI <- lfosr3s(formula = cca ~ case + scan_dates + sex + 
                     age_baseline + (1 + scan_dates | ID), data = DTI_full, 
                   family = "gaussian", var = TRUE, analytic = TRUE)


################################################################################
## plot the results
################################################################################
## Figure 1
plot.DTI <- list()
for(id in c(2017, 2085)){
  subject.id <- id
  sample <- DTI[which(DTI$ID == subject.id),]
  sample.long <- data.frame(ID = subject.id, 
                            visit = as.factor(rep(unique(sample$visit), each = ncol(sample$cca))), 
                            s = rep(1:ncol(sample$cca), length(unique(sample$visit))),
                            FA.cca = c(t(sample$cca)))
  
  plot.DTI[[as.character(id)]] <- ggplot(sample.long) +
    geom_line(aes(x = s, y = FA.cca, color = visit, group = visit)) +
    labs(x = "Tract location", y = "Fractional anisotropy", 
         title = paste0("ID ", unique(sample.long$ID)) ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ylim(1.01*min(sample.long$FA.cca)-0.01*max(sample.long$FA.cca), 
         1.01*max(sample.long$FA.cca)-0.01*min(sample.long$FA.cca) ) +
    scale_x_continuous(breaks = c(1, 24, 47, 70, 93)) +
    scale_color_manual(values=c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#377EB8", "#FB9A99"))
  
}
do.call("grid.arrange", c(plot.DTI, ncol = 2))


## Figure 4
plot.FUI <- function(r, name = NULL){
  var_name <- c("Intercept", "Case", "Scan Date (yr)", "Sex", "Age at Baseline (yr)")
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_DTI$betaHat), length.out = ncol(fit_DTI$betaHat)), 
                             beta = fit_DTI$betaHat[r,],
                             lower = fit_DTI$betaHat[r,] - 2*sqrt(diag(fit_DTI$betaHat.var[,,r])),
                             upper = fit_DTI$betaHat[r,] + 2*sqrt(diag(fit_DTI$betaHat.var[,,r])),
                             lower.joint = fit_DTI$betaHat[r,] - fit_DTI$qn[r]*sqrt(diag(fit_DTI$betaHat.var[,,r])),
                             upper.joint = fit_DTI$betaHat[r,] + fit_DTI$qn[r]*sqrt(diag(fit_DTI$betaHat.var[,,r])))
  p.CI <- ggplot() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray30", alpha = 0.2) +
    geom_ribbon(aes(x = s, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) +
    scale_x_continuous(breaks = c(1, 24, 47, 70, 93))
  
  if(r == 1){
    p.CI <- p.CI + labs(x = "Tract location", y = expression(beta[0](s)), title = var_name[r]) +
      theme(legend.title=element_blank(),
            legend.position = c(0.15, 0.99),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill=alpha('white', 0)))
  }else{
    p.CI <- p.CI + labs(x = "Tract location", y = bquote(paste(beta[.(r-1)], "(s)")), 
                        title = var_name[r]) +
      theme(legend.position = "none") +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted")
  }
  if(!is.null(name)){
    p.CI <- p.CI + labs(title = name)
  }
  return(p.CI)
}

plot.f4 <- list()
for(prdtr in 1:nrow(fit_DTI$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr)
}

do.call("grid.arrange", c(plot.f4, nrow = 1))

