# ------------------------------------------------------------------------#
# Bivariate FPCA basis transformation step:
# ------------------------------------------------------------------------#

# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(ggplot2)    # CRAN v3.3.5
library(fda)        # CRAN v5.5.1
library(GGally)     # CRAN v2.1.2
library(tikzDevice) # CRAN v0.12.3.1
# "R version 4.1.2 (2021-11-01)"

# Settings for ggplot() figures: ------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning()

# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Load some helper functions for functional data manipulation: ------------
source(here::here("code", "functions", "functions-helper-smoothing.R"))

# And helper function for projecting mean function back onto the FPCs:
source(here::here("code","functions","function-project-mean-onto-fpcs.R"))

# Path to save the ouputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")
results_path <- here::here("outputs", "results")

# Read in data ------------------------------------------------------------
# Simplified dataset stored as an rds object.
# functional data is represented in terms of 80 basis coefficients:
data_path <- here::here("data")
subject_side_coef_reg <- readRDS(file.path(data_path, "subject_side_coef_reg.rds"))

# set up basis to represent the data.
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)

# -------------------------------------------------------------------------
# Do mfpca ----------------------------------------------------------------
# To do this we need to put the coefficcients into a 3-d array
# of 
# N_curve x N_coef x N_dimensions
# (i.e., hip and knee coefs are in slices of array)
# First, lets split them into separate FDA objects and do exploratory plots:


## Split coefficients into hip and knee: ----------------------------------
subject_side_coef_hip <- subject_side_coef_reg[location== "Hip"]
subject_side_coef_knee <- subject_side_coef_reg[location== "Knee"]

# check that there is biulateral observations from all subjects:
stopifnot(subject_side_coef_hip[, .(bilateral = .N == 2),, by = subject_id][, all(bilateral)])
stopifnot(subject_side_coef_knee[, .(bilateral = .N == 2),, by = subject_id][, all(bilateral)])

# Print number of subjects being analysed:
if(subject_side_coef_hip[, uniqueN(subject_id)] == subject_side_coef_knee[, uniqueN(subject_id)]) {
  print(paste("Number of Subjects:", subject_side_coef_hip[, uniqueN(subject_id)]))
} else stop()

# and check that the data that are split are in the same ordering:
stopifnot(subject_side_coef_hip$subject_id == subject_side_coef_knee$subject_id)
stopifnot(subject_side_coef_hip$side == subject_side_coef_knee$side)



#                   # Create two univariate fd objects  --------------------------------------
# This is just for visualisation purposes:
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



## Plot a subset of adjacent B-Spline coefficients: -----------------------
# Here we show that adjacent coefficients are highly correlated and it
# would not make sense to model them independently, hence the need
# for the second-stage transformation:
coef_plot_subset <- subject_side_coef_hip[location == "Hip", paste0("lm_coef_", 1:5)]
names(coef_plot_subset) <- paste0("Coef. ", 1:5)
ggpairs(coef_plot_subset) +
  theme(strip.text = element_text())



## Put data in 3-d  array to make a multivariate FDA objects: -------------
# 3-d array with:
# rows = basis coefficients
# columns = observations
# slices = dimensions (hip and knee)
# but instead let's use `fda` package which requires a multivariate fd object
# so make array of coefs
coef_array <- array(data = NA, dim = c(dim(subject_side_fd_hip$coefs), 2))
coef_array[,, 1] <- subject_side_fd_hip$coefs
coef_array[,, 2] <- subject_side_fd_knee$coefs


## Make bfd object --------------------------------------------------------
bfd_obj <- fd(coef = coef_array, basisobj = bspl80)
bfd_obj$fdnames$reps <- paste(subject_side_coef_hip$subject_id, subject_side_coef_hip$side,sep = "_")
bfd_obj$fdnames$funs <- c("Hip", "Knee")
# Plot to check:
par(mfrow = c(1, 2))
plot(bfd_obj)


# Set up functional parameter object for the FPCA:
# - use same basis 80 bsplines
# - use no smoothing (lambda = 0)
harm_fd_par <- fdPar(fdobj = bspl80, Lfdobj = int2Lfd(2), lambda = 0)
bfpca <- pca.fd(fdobj = bfd_obj, nharm = 50, harmfdPar = harm_fd_par)

var_cutoff <- 0.9999
# We want to retain enough fpcas s.t. > 99.99% of var is retained (Zhu et al.).
k_retain <- min(which(cumsum(bfpca$varprop) > var_cutoff)) # This is \widetilde{K} in the paper:
colours <- c(rep("darkgreen", k_retain), rep("red4", 50 - k_retain))



# Plot retained PCs -------------------------------------------------------
plot(bfpca$values[seq_len(50)],
     type = "b",
     xlab = "FPC Number",
     ylab = expression(lambda[k]),
     col = colours, pch = 20, cex = 0.8)
legend("topright", legend = c("Retained", "Not Retained"), col = c("darkgreen", "red4"), pch = 20)
plot(cumsum(bfpca$varprop),
     type = "b", xlab = "FPC Number",
     ylab = "Cumulative Prop. of Variance Explained",
     col = colours, pch = 20, cex = 0.6)
abline(h = c(var_cutoff), lty = 2, col = 1)
abline(v = k_retain, lty = 3, col = 1)
legend( "bottomright", legend = c("PVE = 0.9999", paste("K =", k_retain)), lty = c(2, 3), cex = 0.9)
# dev.off()


# Make publish-able plot: -------------------------------------------------
bfpc_dt <- data.table(
  fpc_number = seq_len(50),
  eigen_value = bfpca$values[seq_len(50)],
  var_pc = 100 * bfpca$varprop[seq_len(50)],
  retained = c(rep("Retained", k_retain), rep("Not Retained", 50 - k_retain))
)
bfpc_dt[, var_pc_cum := cumsum(var_pc)]

p1 <- ggplot(data = bfpc_dt) +
  aes(x = fpc_number, y = eigen_value, colour = retained) +
  geom_line() +
  theme(legend.title = element_blank()) +
  geom_point(size = 0.3) +
  labs(y = "Eigenvalue ($\\lambda_k$)",
       x= "FPC Number ($k$)",
       title = "\\textbf{(a)} Scree Plot")

p2 <- ggplot(data = bfpc_dt) +
  geom_hline(yintercept = 100 *c(0.95, 0.99), lty = 2:3) +
  aes(x = fpc_number, y = var_pc_cum, colour = retained) +
  geom_line() +
  theme(legend.title = element_blank()) +
  geom_point(size = 0.3) +
  labs(x= "FPC Number ($k$)",
       y = "$(\\%)$ of Variance Explained",
       title = "\\textbf{(b)} Explained Variance")

(variance_explained_plot <- ggpubr::ggarrange(plotlist = list(p1, p2),
                                              ncol = 2, nrow = 1,
                                              common.legend = TRUE, 
                                              legend = "bottom")) 


# Do bfpca after choosing \tilde{K}: --------------------------------------
bfpca <- pca.fd(fdobj = bfd_obj,
                nharm = k_retain,
                harmfdPar = harm_fd_par)
bfpca_eval_array <- eval.fd(0:100, bfpca$harmonics)
bfpca_eval_mat <- rbind(bfpca_eval_array[,,1], bfpca_eval_array[,,2])



# Project mean back onto FPCs: --------------------------------------------
mean_scores <- project_mean_onto_fpcs(pca.fd_obj = bfpca)
mean_scores_vec <- apply(mean_scores,c(1,2),sum)[1, ]
mean_tilde <- (matrix(mean_scores_vec,nrow = 1, ncol = k_retain) %*% t(bfpca_eval_mat))[1, ]
mean_eval_array <- eval.fd(0:100, bfpca$meanfd)  
mean_eval_mat <- c(mean_eval_array[,,1], mean_eval_array[,,2])



# Check mean is reconstructed ok: -----------------------------------------
par(mfrow = c(1, 1))
plot(mean_eval_mat, type = "l")
lines(mean_tilde, type = "l", col = 2)



# Add mean scores to matrix of scores: ------------------------------------
# Think of this as ``un-centering" the FPCs:
scores_centred <- apply(bfpca$scores,
                MARGIN = c(1, 2),
                FUN = sum)

scores_uncentred <- sweep(scores_centred,
                        MARGIN = c(2), 
                        STATS = mean_scores_vec, # add on mean vector to each row!
                        FUN = "+",
                        check.margin = TRUE)



# Check reconstruction of indicidual functions: ---------------------------
reconstructed_functions <- scores_uncentred %*% t(bfpca_eval_mat)
true_functions <- rbind(eval.fd(0:100, subject_side_fd_hip),
                     eval.fd(0:100, subject_side_fd_knee))

# Look at random sample for plot:
set.seed(1) # so this sample will be the same every time!
sample_inds <- sample(seq_len(288),  5)
reconstruction_sample <- t(reconstructed_functions)[, sample_inds]
true_sample <- true_functions[, sample_inds]


colnames(true_sample) <- colnames(reconstruction_sample) <- paste0("rep", sample_inds)
reconstruction_dt <- data.table(
  t = rep(0:100, 2),
  dimension = rep(c("Hip", "Knee"), each = 101),
  reconstruction_sample
)
truth_dt <- data.table(
  t = rep(0:100, 2),
  dimension = rep(c("Hip", "Knee"), each = 101),
  true_sample
)

reconstruction_dt$reconstruct <- "Reconstructed"
truth_dt$reconstruct <- "True"

plot_reconstruction_dt <- rbind(
  reconstruction_dt, truth_dt
)

plot_reconstruction_dt_lng <- melt.data.table(
  data = plot_reconstruction_dt, 
  id.vars = c("t", "dimension", "reconstruct"),
  value.name = "angle",
  variable.name = "rep",
  paste0("rep", sample_inds),
  variable.factor = FALSE, 
  value.factor = FALSE, 
  verbose = TRUE
)
plot_reconstruction_dt_lng[, rep:= factor(stringr::str_remove(rep, "rep"))]



(reconstruction_plot <- ggplot(data = plot_reconstruction_dt_lng) +
  aes(x = t,
      y = angle,
      colour = rep, 
      lty = reconstruct,
      alpha = reconstruct,
      group = interaction(reconstruct, rep)) +
  facet_wrap( ~ dimension, scales = "free_y") +
  scale_alpha_manual(values = c(0.5,1)) +
  geom_line(lwd = 0.75) +
  labs(x = "Normalised Time ($\\%$ of Stride)", y = "Angle ($^{\\circ}$)",
       title = "\\textbf{(c)} Sample Reconstructions") +
  guides(colour ="none") +
  theme(legend.title = element_blank(),
        legend.margin = margin(b = -1),
        legend.position = "bottom"))




 # Create publishable combined plot: ------------------------------
(basis_transform_plot <- ggpubr::ggarrange(variance_explained_plot, 
                  reconstruction_plot))
tikz(file.path(plots_path, "basis-transform-plot.tex"),
     width = 1.25 * doc_width_inches, 
     height = 1.25 * (doc_width_inches/3))
print(basis_transform_plot)
dev.off()


# Store results: ----------------------------------------------------------

basis_transform_results <- list(
  coef_dts = list(subject_side_coef_hip, subject_side_coef_knee),
  bfd_obj = bfd_obj,
  bfpca_obj = bfpca,
  scores_centred = scores_centred,
  k_retain = k_retain,
  var_cutoff = var_cutoff,
  scores_uncentred = scores_uncentred
)

saveRDS(object = basis_transform_results,
        file = file.path(results_path,"basis-transform-results.rds"))

