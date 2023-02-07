# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0
library(fda)        # CRAN v5.5.1
library(ggplot2)    # CRAN v3.4.0

# Helper functions needed -------------------------------------------------
source(here::here("code", "functions", "functions-helper-smoothing.R"))
source(here::here("code", "functions", "functions-unstructured-covariance.R"))

# Settings for ggplot() figures: ------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning() # set theme
theme_update(panel.grid.major = element_blank(),
             legend.title = element_text(hjust = 0.5),
             legend.key.size = unit(0.95,"line"),
             legend.text = element_text(size = 8, hjust = 0.5))
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937



# Path to save the ouputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")
results_path <- here::here("outputs", "results")



# Import Results: ---------------------------------------------------------
basis_transformation_results <- readRDS(file.path(results_path, "basis-transform-results.rds"))
lme_fit <- readRDS(file = file.path(results_path, "model-fit-results.rds"))


model_formula_fixef <- lme_fit$fixef_formula
bfpca <- basis_transformation_results$bfpca_obj
bfd_obj<- basis_transformation_results$bfd_obj
k_retain <- basis_transformation_results$k_retain
lme_dt <- lme_fit$lme_df


# Extract Covariance Parameters -------------------------------------------
q <- lme_fit$q_vec
Q_star <- diag(q)

s <- lme_fit$s_vec
S_star <- diag(s)

Psi <- rbind(eval.fd(0:100, bfpca$harmonics[,1]), eval.fd(0:100, bfpca$harmonics[,2]))

Q <- Psi %*% Q_star %*% t(Psi)
S <- Psi %*% S_star %*% t(Psi)



# Unstructured covariance estimates: --------------------------------------

# Need to do a least squares fir to "center the data" around fixed effects
# do this under working independence assumption.
model_matrix <- model.matrix(as.formula(paste("~", model_formula_fixef)), data = lme_dt)
x_fd_list <- get_list_of_columns(model_matrix)
names(x_fd_list) <- colnames(model_matrix)
beta_list <- replicate(n = length(x_fd_list), expr = fdPar(fdobj = bfd_obj$basis, Lfdobj = 2, lambda = 0), simplify = FALSE)
freg_hip <- fRegress(y = bfd_obj[,1], xfdlist = x_fd_list, betalist = beta_list)
freg_knee <- fRegress(y = bfd_obj[,2], xfdlist = x_fd_list, betalist = beta_list)



# Get unstructured estimates: ---------------------------------------------
Y_centered_hip <- (freg_hip$yfdobj - freg_hip$yhatfdobj)
Y_centered_knee <- (freg_knee$yfdobj - freg_knee$yhatfdobj)
unsturctured_cov <- cov_unstruc_mlfpca_bi_fd(fd_obj_list = list(
  Y_centered_hip, Y_centered_knee
), id_vec =  as.integer(as.factor(lme_dt$subject_id)))


Q_unstruc <- S_unstruc <- matrix(NA, nrow = 202, ncol = 202)
Q_unstruc[1:101, 1:101] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_U_bifd_list$K_U_11)
Q_unstruc[1:101, 102:202] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_U_bifd_list$K_U_12)
Q_unstruc[102:202, 1:101] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_U_bifd_list$K_U_21)
Q_unstruc[102:202, 102:202] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_U_bifd_list$K_U_22)

S_unstruc[1:101, 1:101] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_E_bifd_list$K_E_11)
S_unstruc[1:101, 102:202] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_E_bifd_list$K_E_12)
S_unstruc[102:202, 1:101] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_E_bifd_list$K_E_21)
S_unstruc[102:202, 102:202] <- eval.bifd(0:100, 0:100, unsturctured_cov$K_E_bifd_list$K_E_22)



# Quick and dirty plots! --------------------------------------------------
par(mfrow = c(1, 2))
persp(S[102:202, 102:202], zlim = range(
  S[102:202, 102:202], S_unstruc[102:202, 102:202]
))
persp(S_unstruc[102:202, 102:202])

persp(Q[102:202, 102:202], zlim = range(
  Q[102:202, 102:202], Q_unstruc[102:202, 102:202]
))
persp(Q_unstruc[102:202, 102:202])




# Reshape data for Publish-able plots -------------------------------------
## Q ----------------------------------------------------------------------
Q_dt <- data.table(
  t = rep(0:100, times = 2),
  dim1 = rep(c("Hip", "Knee"), each = 101),
  Q
)
names(Q_dt)[- c(1:2)] <- paste(Q_dt$t, Q_dt$dim1, sep = "_")
Q_dt_lng <- melt.data.table(data = Q_dt, 
                            id.vars = c("t", "dim1"),
                            measure.vars =  paste(Q_dt$t, Q_dt$dim1, sep = "_"),
                            variable.name = "s_dim2",
                            value.name = "Q_st", 
                            variable.factor = FALSE, 
                            value.factor = FALSE,
                            verbose = TRUE, na.rm = FALSE)
Q_dt_lng[, s := stringr::str_extract(s_dim2, pattern = "\\d{1,3}_")]
Q_dt_lng[, s := as.numeric(stringr::str_remove(s, "_"))]
Q_dt_lng[, dim2 := stringr::str_remove(s_dim2, pattern = "\\d{1,3}_")]
Q_dt_lng[, dim_comb := paste(dim1, dim2, sep = " - ")]


Q_dt_unstruc <- data.table(
  t = rep(0:100, times = 2),
  dim1 = rep(c("Hip", "Knee"), each = 101),
  Q_unstruc
)
names(Q_dt_unstruc)[- c(1:2)] <- paste(Q_dt_unstruc$t, Q_dt_unstruc$dim1, sep = "_")
Q_dt_unstruc_lng <- melt.data.table(data = Q_dt_unstruc, 
                            id.vars = c("t", "dim1"),
                            measure.vars =  paste(Q_dt_unstruc$t, Q_dt_unstruc$dim1, sep = "_"),
                            variable.name = "s_dim2",
                            value.name = "Q_st", 
                            variable.factor = FALSE, 
                            value.factor = FALSE,
                            verbose = TRUE, na.rm = FALSE)
Q_dt_unstruc_lng[, s := stringr::str_extract(s_dim2, pattern = "\\d{1,3}_")]
Q_dt_unstruc_lng[, s := as.numeric(stringr::str_remove(s, "_"))]
Q_dt_unstruc_lng[, dim2 := stringr::str_remove(s_dim2, pattern = "\\d{1,3}_")]
Q_dt_unstruc_lng[, dim_comb := paste(dim1, dim2, sep = " - ")]

zlim_Q <- range(Q_dt_lng$Q_st, Q_dt_unstruc_lng$Q_st)
expand_lim <- c(0.02, 0.02)
breaks_Q <- seq(zlim_Q[1], zlim_Q[2], length.out = 13)

Q_plot <- ggplot(Q_dt_lng) +
  aes(x = s, y = t, z = Q_st) +
  scale_x_continuous(expand = expand_lim) +
  scale_y_continuous(expand = expand_lim) +
  facet_wrap(~ dim_comb) +
  geom_contour_filled(breaks = breaks_Q) +
  labs(x = "$t_1$", y = "$t_2$", fill = "$\\mathbf{Q}(t_1, t_2)$", title = "Model Estimate")


Q_plot_unstruc <- ggplot(Q_dt_unstruc_lng) +
  aes(x = s, y = t, z = Q_st) +
  facet_wrap(~ dim_comb) +
  scale_x_continuous(expand = expand_lim) +
  scale_y_continuous(expand = expand_lim) +
  geom_contour_filled(breaks = breaks_Q) +
  labs(x = "$t_1$", y = "$t_2$", fill = "$\\mathbf{Q}(t_1, t_2)$", title = "Unstructured Estimate")



Q_combined <- ggarrange(Q_plot, Q_plot_unstruc, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")



## S ----------------------------------------------------------------------

S_dt <- data.table(
  t = rep(0:100, times = 2),
  dim1 = rep(c("Hip", "Knee"), each = 101),
  S
)
names(S_dt)[- c(1:2)] <- paste(S_dt$t, S_dt$dim1, sep = "_")
S_dt_lng <- melt.data.table(data = S_dt, 
                            id.vars = c("t", "dim1"),
                            measure.vars =  paste(S_dt$t, S_dt$dim1, sep = "_"),
                            variable.name = "s_dim2",
                            value.name = "S_st", 
                            variable.factor = FALSE, 
                            value.factor = FALSE,
                            verbose = TRUE, na.rm = FALSE)
S_dt_lng[, s := stringr::str_extract(s_dim2, pattern = "\\d{1,3}_")]
S_dt_lng[, s := as.numeric(stringr::str_remove(s, "_"))]
S_dt_lng[, dim2 := stringr::str_remove(s_dim2, pattern = "\\d{1,3}_")]
S_dt_lng[, dim_comb := paste(dim1, dim2, sep = " - ")]


S_dt_unstruc <- data.table(
  t = rep(0:100, times = 2),
  dim1 = rep(c("Hip", "Knee"), each = 101),
  S_unstruc
)
names(S_dt_unstruc)[- c(1:2)] <- paste(S_dt_unstruc$t, S_dt_unstruc$dim1, sep = "_")
S_dt_unstruc_lng <- melt.data.table(data = S_dt_unstruc, 
                                    id.vars = c("t", "dim1"),
                                    measure.vars =  paste(S_dt_unstruc$t, S_dt_unstruc$dim1, sep = "_"),
                                    variable.name = "s_dim2",
                                    value.name = "S_st", 
                                    variable.factor = FALSE, 
                                    value.factor = FALSE,
                                    verbose = TRUE, na.rm = FALSE)
S_dt_unstruc_lng[, s := stringr::str_extract(s_dim2, pattern = "\\d{1,3}_")]
S_dt_unstruc_lng[, s := as.numeric(stringr::str_remove(s, "_"))]
S_dt_unstruc_lng[, dim2 := stringr::str_remove(s_dim2, pattern = "\\d{1,3}_")]
S_dt_unstruc_lng[, dim_comb := paste(dim1, dim2, sep = " - ")]

zlim_S <- range(S_dt_lng$S_st, S_dt_unstruc_lng$S_st)
expand_lim <- c(0.02, 0.02)
breaks_S <- seq(zlim_S[1], zlim_S[2], length.out = 13)

S_plot <- ggplot(S_dt_lng) +
  aes(x = s, y = t, z = S_st) +
  scale_x_continuous(expand = expand_lim) +
  scale_y_continuous(expand = expand_lim) +
  facet_wrap(~ dim_comb) +
  geom_contour_filled(breaks = breaks_S) +
  labs(x = "$t_1$", y = "$t_2$", fill = "$\\mathbf{S}(t_1, t_2)$", title = "Model Estimate")



S_plot_unstruc <- ggplot(S_dt_unstruc_lng) +
  aes(x = s, y = t, z = S_st) +
  facet_wrap(~ dim_comb) +
  scale_x_continuous(expand = expand_lim) +
  scale_y_continuous(expand = expand_lim) +
  geom_contour_filled(breaks = breaks_S) +
  labs(x = "$t_1$", y = "$t_2$", fill = "$\\mathbf{S}(t_1, t_2)$", title = "Unstructured Estimate")



S_combined <- ggarrange(S_plot, S_plot_unstruc, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")




# Create plot and save: ---------------------------------------------------
full_plot <- ggarrange(Q_combined, S_combined, ncol = 1, nrow = 2)
full_plot

tikz(file.path(plots_path, "covariance-reconstruction-plot.tex"),
     width = 1.25 * doc_width_inches, 
     height = 1.1 * (doc_width_inches))
print(full_plot)
dev.off()



# Save results for use later: ---------------------------------------------
saveRDS(list(S = S,
             S_unstruc = S_unstruc,
             Q = Q,
             Q_unstruc = Q_unstruc), 
        file = file.path(results_path, "covariance-results.rds"))
