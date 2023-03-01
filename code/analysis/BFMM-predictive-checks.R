# ------------------------------------------------------------------------#
# Plots of model BLUp predictions. Compare to data to make sure 
# individual predictions from the model are reasonable.
# ------------------------------------------------------------------------#

# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(fda)        # CRAN v5.5.1
library(lme4)       # CRAN v1.1-30

# Load results: -----------------------------------------------------------
results_path <- here::here("outputs", "results")

basis_transformation_results <- readRDS(
  file = file.path(results_path, "basis-transform-results.rds"))

lme_fit <- readRDS(file.path(results_path, "model-fit-results.rds"))

# -------------------------------------------------------------------------
inds <- 1:6
predictions_scores <- sapply(lme_fit$full_list, function(x) {
  predict(x)[inds]
})
bfpca <- basis_transformation_results$bfpca_obj
bfd_obj <- basis_transformation_results$bfd_obj
Psi <- rbind(eval.fd(0:100, bfpca$harmonics[,1]),
             eval.fd(0:100, bfpca$harmonics[,2]))
mean_eval_vec <- basis_transformation_results$mean_eval_vec
predictions <- predictions_scores %*% t(Psi)
predictions <- sweep(predictions, MARGIN = 2, mean_eval_vec, FUN = "+", check.margin = TRUE)

matplot(t(predictions)[1:101, ], type = "l", col = rep(1:3, each = 2), lty=1)
lines(bfd_obj[inds,1], col = rep(1:3, each = 2), lty = 2)
matplot(t(predictions)[102:202, ], type = "l", col = rep(1:3, each = 2), lty = 1)
lines(bfd_obj[inds,2], col = rep(1:3, each = 2), lty = 2)


# -------------------------------------------------------------------------
inds <- 13:18
predictions_scores <- sapply(lme_fit$full_list, function(x) {
  predict(x)[inds]
})
bfpca <- basis_transformation_results$bfpca_obj
bfd_obj <- basis_transformation_results$bfd_obj
Psi <- rbind(eval.fd(0:100, bfpca$harmonics[,1]),
             eval.fd(0:100, bfpca$harmonics[,2]))
mean_eval_vec <- basis_transformation_results$mean_eval_vec
predictions <- predictions_scores %*% t(Psi)
predictions <- sweep(predictions, MARGIN = 2, mean_eval_vec, FUN = "+", check.margin = TRUE)

matplot(t(predictions)[1:101, ], type = "l", col = rep(1:3, each = 2), lty=1)
lines(bfd_obj[inds,1], col = rep(1:3, each = 2), lty = 2)
matplot(t(predictions)[102:202, ], type = "l", col = rep(1:3, each = 2), lty = 1)
lines(bfd_obj[inds,2], col = rep(1:3, each = 2), lty = 2)


# -------------------------------------------------------------------------

inds <- 19:24
predictions_scores <- sapply(lme_fit$full_list, function(x) {
  predict(x)[inds]
})
bfpca <- basis_transformation_results$bfpca_obj
bfd_obj <- basis_transformation_results$bfd_obj
Psi <- rbind(eval.fd(0:100, bfpca$harmonics[,1]),
             eval.fd(0:100, bfpca$harmonics[,2]))
mean_eval_vec <- basis_transformation_results$mean_eval_vec
predictions <- predictions_scores %*% t(Psi)
predictions <- sweep(predictions, MARGIN = 2, mean_eval_vec, FUN = "+", check.margin = TRUE)

matplot(t(predictions)[1:101, ], type = "l", col = rep(1:3, each = 2), lty=1)
lines(bfd_obj[inds,1], col = rep(1:3, each = 2), lty = 2)
matplot(t(predictions)[102:202, ], type = "l", col = rep(1:3, each = 2), lty = 1)
lines(bfd_obj[inds,2], col = rep(1:3, each = 2), lty = 2)


# -------------------------------------------------------------------------
# all look good!
