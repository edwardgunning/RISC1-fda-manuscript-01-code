# ------------------------------------------------------------------------#
# function to calculate binomial SE.
# Useful for calculating the Monte-Carlo Standard Errors for 
# confidence interval coverage as described in:

# Morris, T. P., White, I. R., & Crowther, M. J. (2019). 
# Using simulation studies to evaluate statistical methods.
# Statistics in Medicine, 38(11), 2074–2102. https://doi.org/10.1002/sim.8086

binomial_se <- function(phat, n) {
  sqrt(
    (phat  * (1- phat)) / n
  )
}