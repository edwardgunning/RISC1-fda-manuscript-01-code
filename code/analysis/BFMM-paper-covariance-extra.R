# ------------------------------------------------------------------------#
# Extra analysis on reconstructing covariance function.
# Use methods of Fieuws, S., & Verbeke, G. (2006) to estimate
# unstructured S^* and GLASSO to selet which off-diagonal elements
# ------------------------------------------------------------------------#


# Packages ----------------------------------------------------------------
library(nlme)       # CRAN v3.1-155
library(mixedup)    # [github::m-clark/mixedup] v0.3.9
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0
library(fda)        # CRAN v5.5.1
library(ggplot2)    # CRAN v3.3.5
library(glasso)     # CRAN v1.11 
library(parallel)

ncores <- detectCores()


# Settings for ggplot() figures: ------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning() # set theme
theme_update(panel.grid.major = element_blank(),
             legend.title = element_text(hjust = 0.5),
             legend.key.size = unit(0.95,"line"),
             legend.text = element_text(size = 8, hjust = 0.5))
# and rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937



# -------------------------------------------------------------------------
# Load helper functions ---------------------------------------------------
source(here::here("code", "functions", "functions-helper-smoothing.R"))
source(here::here("code", "functions", "functions-unstructured-covariance.R"))
source(here::here("code", "functions", "function-get-residual-covariance-matrix.R"))

# Path to save the outputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")
results_path <- here::here("outputs", "results")


# Load results of covariance analysis: ------------------------------------
covariance_results <- readRDS(file.path(results_path, "covariance-results.rds"))

S <- covariance_results$S
Q <- covariance_results$Q
S_unstruc <- covariance_results$S_unstruc
Q_unstruc <- covariance_results$Q_unstruc


# Load results of basis transformation ------------------------------------
basis_transformation_results <- readRDS(file.path(results_path, "basis-transform-results.rds"))
bfpca <- basis_transformation_results$bfpca_obj
bfd_obj<- basis_transformation_results$bfd_obj
k_retain <- basis_transformation_results$k_retain

# Load results of mixed modelling: ----------------------------------------
lme_fit <- readRDS(file.path(results_path, "model-fit-results.rds"))
model_formula_fixef <- lme_fit$fixef_formula
lme_dt <- lme_fit$lme_df


# Extract results ---------------------------------------------------------
q <- lme_fit$q_vec
Q_star <- diag(q)

s <- lme_fit$s_vec
S_star <- diag(s)

Psi <- rbind(eval.fd(0:100, bfpca$harmonics[,1]), 
             eval.fd(0:100, bfpca$harmonics[,2]))



# Fit model with unstructured $S^*$ ---------------------------------------
# use the algorithm of:

# Fieuws, S., & Verbeke, G. (2006). 
# Pairwise Fitting of Mixed Models for the Joint Modeling of Multivariate Longitudinal Profiles.
# Biometrics, 62(2), 424–431. https://doi.org/10.1111/j.1541-0420.2006.00507.x

# which fits all combinations of pairwise mixed moels

# need to go from wide to long so that each score is an observation.
# but need to add an id for each observation
lme_dt$obs <- seq_len(nrow(lme_dt))
id_vars <- c(
  "ris", "sex", "side", "speed_cent", 
  "weight_kg_cent", "age_cent", 
  "height_cm_cent", "runner_category", 
  "treadmill", "dominance", 
  "bmi_kgm_cent", "subject_id", 
  "side_is_dominant", "obs")

lme_dt_long <- melt.data.table(data = lme_dt,
                               verbose = TRUE,
                               measure.vars = paste0("score_", 1:k_retain),
                               id.vars = id_vars,
                               variable.factor = TRUE,
                               variable.name = "dimension", 
                               value.name = "score",
                               value.factor = FALSE)
lme_dt_long[, dimension := stringr::str_remove(dimension, "score_")]
lme_dt_long[, dimension := factor(dimension, levels = seq_len(k_retain))]



# List of all possible pairs of scores:
dim_list <- combn(seq_len(k_retain), m = 2, simplify = FALSE)
head(dim_list)

# fit all possible pairs of bivariate models in parallel.
# (note we have had to specify the model manually, there should be a way
# to specify model using model_formula_fixed)
system.time(model_list <- parallel::mclapply(dim_list, function(d) {
  lme_dt_long_subset <- lme_dt_long[dimension %in% d]
  lme(score ~ - 1 + dimension + ris:dimension + sex:dimension + speed_cent:dimension +
        weight_kg_cent:dimension + height_cm_cent:dimension + age_cent:dimension,
      random = list(subject_id = pdDiag(form = ~ 0 + dimension)),
      weights = varIdent(form = ~ 1 | dimension),
      correlation = corSymm(form = ~ 1 | subject_id/obs),
      data = lme_dt_long_subset,
      method = "REML",
      control = lmeControl(msMaxIter = 1000, msMaxEval = 1000))
}, mc.cores = ncores))


# Extract estimated parameters --------------------------------------------
residual_df_pairwise <- purrr::map_dfr(model_list, 
                                       .f = function(x) {
                                         residual_covariance_matrix_pair <- get_residual_covariance_matrix(x)
                                         dimensions_pairs <- unique(x$data$dimension)
                                         terms_names_pairs <- c(as.character(dimensions_pairs), paste(dimensions_pairs, collapse = ","))
                                         data.frame(parameter = paste0("s_", terms_names_pairs),
                                                    value = c(diag(residual_covariance_matrix_pair), 
                                                              residual_covariance_matrix_pair[1,2]))
                                         },
                                       .id = "Model")

residual_dt_pairwise <- data.table(residual_df_pairwise)
# Because some parameters are estimated twice, we need to average them.
residual_dt_estimates <- residual_dt_pairwise[, .(value = mean(value)), by = parameter]

# and put back into matrix to give S^*
S_star_full <- matrix(data = NA, nrow = k_retain, ncol = k_retain)
for(i in 1:k_retain) {
  S_star_full[i, i] <- residual_dt_estimates[parameter == paste0("s_", i), value]
}

for(i in seq_len(k_retain - 1)) {
  for(j in (i+1):(k_retain)) {
    S_star_full[i, j] <- S_star_full[j, i] <- residual_dt_estimates[parameter == paste0("s_", i, ",", j), value]
  }
}



# Plot of covariance and correlation matrices: ----------------------------

S_star_full_sample <- S_star_full[1:8, 1:8] # only look at first 8:

# reshape data for correlation and covariance plots:
colnames(S_star_full_sample) <- rownames(S_star_full_sample) <- paste0("FPC ", 1:8)
S_star_corr_sample <- cov2cor(S_star_full_sample)

# covariance reshaping:
S_star_full_sample_dt <- data.table(Var1 = paste0("FPC ", 1:8),
                                    S_star_full_sample)
S_star_full_sample_dt_lng <- melt.data.table(data = S_star_full_sample_dt,
                                             id.vars = c("Var1"),
                                             measure.vars = paste0("FPC ", 1:8),
                                             variable.name = "Var2",
                                             value.name = "cov",
                                             variable.factor = FALSE,
                                             value.factor = FALSE,
                                             verbose = TRUE
                                             )
S_star_full_sample_dt_lng[, Var1 := factor(Var1, levels = paste0("FPC ", 1:8))]
S_star_full_sample_dt_lng[, Var2 := factor(Var2, levels = paste0("FPC ", 1:8))]

# correlation reshaping
diag(S_star_corr_sample) <- NA # want to fill diagonal with NAs on plot.
S_star_corr_sample_dt <- data.table(Var1 = paste0("FPC ", 1:8),
                                    S_star_corr_sample)

S_star_corr_sample_dt_lng <- melt.data.table(data = S_star_corr_sample_dt,
                                             id.vars = c("Var1"),
                                             measure.vars = paste0("FPC ", 1:8),
                                             variable.name = "Var2",
                                             value.name = "cor",
                                             variable.factor = FALSE,
                                             value.factor = FALSE,
                                             verbose = TRUE
)
S_star_corr_sample_dt_lng[, Var1 := factor(Var1, levels = paste0("FPC ", 1:8))]
S_star_corr_sample_dt_lng[, Var2 := factor(Var2, levels = paste0("FPC ", 1:8))]

S_cor_plot <- ggplot(data = S_star_corr_sample_dt_lng,
       aes(x=Var1, y=Var2, fill=cor)) + 
  geom_tile(linejoin = "round", ) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="$s_{kk^{'}} / s_{kk^{'}}$", na.value = "lightgrey",
                       guide = guide_colorbar(label = TRUE,
                                              draw.ulim = TRUE, 
                                              draw.llim = TRUE,
                                              # here comes the code change:
                                              frame.colour = "grey",
                                              ticks = TRUE, 
                                              nbin = 10,
                                              label.position = "bottom",
                                              barwidth = 6,
                                              barheight = 0.6, 
                                              direction = 'horizontal'))+
  geom_text(aes(Var2, Var1, label = round(cor, 2)), color = "black", size = 2) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.75),
        panel.border = element_rect(colour = "grey"),
        axis.title = element_blank(), 
        legend.title = element_blank(),
        legend.margin = margin(t = 0),
        legend.position = "bottom")  +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Correlation Matrix")

S_cov_plot <- ggplot(data = S_star_full_sample_dt_lng,
                     aes(x=Var1, y=Var2, fill=cov)) + 
  geom_tile(linejoin = "round") +
  scale_fill_gradient2(low = "white", high = "dodgerblue3", mid = "white", 
                       midpoint = median(S_star_full_sample_dt_lng$cov),
                       limit = range(S_star_full_sample_dt_lng$cov), space = "Lab", 
                       name="$s_{kk^{'}}$:",
                       guide = guide_colorbar(label = TRUE,
                                              draw.ulim = TRUE, 
                                              draw.llim = TRUE,
                                              # here comes the code change:
                                              frame.colour = "grey",
                                              ticks = TRUE, 
                                              nbin = 10,
                                              label.position = "bottom",
                                              barwidth = 6,
                                              barheight = 0.6, 
                                              direction = 'horizontal')) +
  geom_text(aes(Var2, Var1, label = round(cov, 2)), color = "black", size = 2) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.75),
        panel.border = element_rect(colour = "grey"),
        axis.title = element_blank(), 
        legend.margin = margin(t = 0),
        legend.title = element_blank(),
        legend.position = "bottom") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Covariance Matrix")


ggarrange(S_cov_plot, S_cor_plot)

tikz(file.path(plots_path, "S-star-cov-and-cor-plot.tex"),
     width = 0.8 * doc_width_inches, 
     height = 0.5 * doc_width_inches)
ggarrange(S_cov_plot, S_cor_plot)
dev.off()



# Look at Simpson's Paradox Phenom: -----------------------------------------
# (not included in paper)
# overall, scores are uncorrelated:
ggplot(lme_dt) +
  aes(x = score_1, y = score_3, colour = subject_id) +
  geom_point() +
  theme(legend.position = "none")+
  geom_line(aes(group = subject_id))

# Correlated after being centered around subject's mean!
lme_dt_simpson <- copy(lme_dt)
lme_dt_simpson[, score_1_mc := score_1 - mean(score_1), by = subject_id]
lme_dt_simpson[, score_3_mc := score_3 - mean(score_3), by = subject_id]

ggplot(lme_dt_simpson) +
  aes(x = score_1_mc, y = score_3_mc, colour = subject_id) +
  geom_point() +
  theme(legend.position = "none")


cor(lme_dt_simpson$score_1, lme_dt_simpson$score_3)
cor(lme_dt_simpson$score_3_mc, lme_dt_simpson$score_1_mc)




# Making a Sparse S^* with the Graphical Lasso ----------------------------

# Use the glassopath()
# Going to apply Graphical Lasso to S^*
# Arguments:
# 
# s	= Covariance matrix:p by p matrix (symmetric): 
# - We plug in S^*
# 
# 
# rholist	= Vector of non-negative regularization parameters for
# the lasso. Should be increasing from smallest to largest;
# actual path is computed from largest to smallest value of rho).
# If NULL, 10 values in a (hopefully reasonable) range are used. 
# Note that the same parameter rholist[j] is used for all entries
# of the inverse covariance matrix; different penalties for 
# different entries are not allowed: 
# - We'll do trial and error for a range of values.
#
# penalize.diagonal = Should diagonal of inverse covariance be
# penalized? Default TRUE:
# - We set this to FALSE!

# -------------------------------------------------------------------------
glasso_path <- glasso::glassopath(s = S_star_full,
                                  rholist = seq(0, 200, by = 1),
                                  penalize.diagonal = FALSE)
# check no errors have been flagged.
stopifnot(glasso_path$errflag == 0)
length(seq(0, 100, by = 1))
# w	contains:
# 'Estimated covariance matrices, an array of dimension
# (nrow(s),ncol(n), length(rholist))'
# sowe will look at this,

# count the number of non-zero elements in the graphical lasso.
non_zero <- apply(glasso_path$w, MARGIN = 3, FUN = function(x) {
  sum(x[lower.tri(x)] != 0)
})

par(mfrow = c(1, 1))
plot(seq(0, 200, by = 1),
     non_zero,
     type = "b",
     pch = 20,
     xlab = "Tuning Parameter $\\rho$",
     ylab = "Non-Zero Elements")

ise_vec <- vector("numeric", 201)
for(i in seq_len(201)) {
  S_star_regularised <- glasso_path$w[,, i]
  S_star_tmp <- S_star_full
  S_star_tmp[S_star_regularised == 0] <- 0
  S_regularised <- Psi %*% S_star_tmp %*% t(Psi)
  ise_vec[i] <- mean((S_regularised - S_unstruc)^2)
}

# Join back on the error from the diagonal s:
non_zero <- c(non_zero, 0)
ise_vec <- c(ise_vec, mean((S-S_unstruc)^2))

# quick and dirty plot
plot(non_zero, ise_vec, xlim = c(0, 40), type = "b", pch = 20)


# Create Publishable Plot -------------------------------------------------

tuning_param_dt <- data.table(
  non_zero_elements = non_zero[-length(non_zero)],
  tuning_parameter = seq(0, 200, by = 1)
)

ise_dt <- data.table(
  non_zero_elements = non_zero,
  ise = ise_vec
)


non_zero_label <- "No. of Non-Zero Elements in $\\mathbf{S^*}$"

theme_gunning() # set theme

p1 <- ggplot(data = tuning_param_dt) +
  aes(x = tuning_parameter, y = non_zero_elements) +
  geom_point(size = 0.5) +
  geom_line() +
  labs(x = "Graphical LASSO Tuning Parameter ($\\rho$)",
       y = non_zero_label,
       title = "Variable Selection") +
  ylim(c(0, 200))



p2 <- ggplot(data = ise_dt) +
  aes(x = non_zero_elements, y = ise) +
  geom_point(size = 0.5) +
  geom_line() +
  xlim(c(0, 45)) +
  labs(x = non_zero_label,
       y = "(Scaled) Reconstruction Error of $\\mathbf{S}(t_1, t_2)$",
       title = "Reconstruction Error")

tikz(file.path(plots_path, "glasso-selection-plot.tex"),
     width = 1 * doc_width_inches, 
     height = 0.5 * doc_width_inches)
ggarrange(p1, p2, ncol = 2, nrow = 1)
dev.off()




# Now, final plot: --------------------------------------------------------
# ------------------------------------------------------------------------#


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







# Sparse reconstruction ---------------------------------------------------

ind <- min(which(non_zero == 3))
S_star_regularised <- glasso_path$w[,, ind]

S_star_sparse <- S_star_full
S_star_sparse[S_star_regularised == 0] <- 0

S_sparse_fun <- Psi %*% S_star_sparse %*% t(Psi)

S_dt_sparse <- data.table(
  t = rep(0:100, times = 2),
  dim1 = rep(c("Hip", "Knee"), each = 101),
  S_sparse_fun
)
names(S_dt_sparse)[- c(1:2)] <- paste(S_dt_sparse$t, S_dt_sparse$dim1, sep = "_")
S_dt_sparse_lng <- melt.data.table(data = S_dt_sparse, 
                                    id.vars = c("t", "dim1"),
                                    measure.vars =  paste(S_dt_sparse$t, S_dt_sparse$dim1, sep = "_"),
                                    variable.name = "s_dim2",
                                    value.name = "S_st", 
                                    variable.factor = FALSE, 
                                    value.factor = FALSE,
                                    verbose = TRUE, na.rm = FALSE)
S_dt_sparse_lng[, s := stringr::str_extract(s_dim2, pattern = "\\d{1,3}_")]
S_dt_sparse_lng[, s := as.numeric(stringr::str_remove(s, "_"))]
S_dt_sparse_lng[, dim2 := stringr::str_remove(s_dim2, pattern = "\\d{1,3}_")]
S_dt_sparse_lng[, dim_comb := paste(dim1, dim2, sep = " - ")]

zlim_S <- range(S_dt_lng$S_st, S_dt_unstruc_lng$S_st, S_dt_sparse_lng$S_st)
expand_lim <- c(0.02, 0.02)
breaks_S <- seq(zlim_S[1], zlim_S[2], length.out = 13)

theme_update(panel.grid.major = element_blank(),
             legend.title = element_text(hjust = 0.5),
             legend.key.size = unit(0.95,"line"),
             legend.text = element_text(size = 8, hjust = 0.5))

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


S_plot_sparse <- ggplot(S_dt_sparse_lng) +
  aes(x = s, y = t, z = S_st) +
  facet_wrap(~ dim_comb) +
  scale_x_continuous(expand = expand_lim) +
  scale_y_continuous(expand = expand_lim) +
  geom_contour_filled(breaks = breaks_S) +
  labs(x = "$t_1$", y = "$t_2$", fill = "$\\mathbf{S}(t_1, t_2)$", title = "Sparse Estimate")



S_combined <- ggarrange(S_plot, S_plot_unstruc, S_plot_sparse, nrow = 1, ncol = 3, common.legend = TRUE, legend = "bottom")

tikz(file.path(plots_path, "covariance-reconstruction-plot-extra.tex"),
     width = 1.25 * doc_width_inches * (3/2), 
     height = 0.75 * (doc_width_inches))
print(S_combined)
dev.off()


