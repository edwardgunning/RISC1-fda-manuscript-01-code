library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0
library(fda)        # CRAN v5.5.1
library(ggplot2)    # CRAN v3.4.0
library(lme4)       # CRAN v1.1-30
source(here::here("code", "BFMM-paper-helper-functions.R"))

# Settings for ggplot() figures: ------------------------------------------
theme_set(theme_bw())
theme_update(strip.text = element_text(size = 9),
             text = element_text(size = 9),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             axis.title = element_text(size = 9),
             legend.title = element_text(hjust = 0.5),
             legend.key.size = unit(0.95,"line"),
             legend.text = element_text(size = 8, hjust = 0.5),
             plot.title = element_text(hjust = 0.5, size = 9.5, face = "bold"))

# Path to save the ouputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "plots", "subject-side")



# -------------------------------------------------------------------------

basis_transformation_results <- readRDS(
  file = here::here("outputs",
                    "full-data",
                    "RDS-objects", 
                    "basis-transform-results.rds"))


lme_fit <- readRDS(
  file = here::here("outputs",
                    "full-data",
                    "RDS-objects", 
                    "model-fit-results.rds")
)
# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937


# Extract things from fitted model:
model_formula_fixef <- lme_fit$fixef_formula
lme_dt <- lme_fit$lme_df
bfpca <- basis_transformation_results$bfpca_obj
bfd_obj<- basis_transformation_results$bfd_obj
k_retain <- basis_transformation_results$k_retain
Psi <- rbind(eval.fd(0:100, bfpca$harmonics[,1]),
             eval.fd(0:100, bfpca$harmonics[,2]))
fixef_coef_point <- lme_fit$fixef_mat
fixef_coef_var <- lme_fit$fixef_var
fixef_names <- names(fixef(lme_fit$full_list[[1]]))
names(fixef_names) <- fixef_names

# --------------------------------------------------------------
set.seed(1)
boot_result <- bootstrap_of_subjects_coefs(lme_df = lme_dt,
                            fixef_formula = model_formula_fixef,
                            k_retain = k_retain,
                            REML = TRUE,
                            B = 2500, # Degras (2017) 
                            par_mc = TRUE,
                            n_cores = 4)

# -------------------------------------------------------------------------

fixef_fun_point <- fixef_coef_point %*% t(Psi)

fixef_fun_covar <- lapply(fixef_names, function(x) {
  Psi %*% diag(fixef_coef_var[x, ]) %*% t(Psi)
})

fixef_fun_se <- t(sapply(fixef_fun_covar, function(x) {
  sqrt(diag(x))
}))

fixef_coef_samples_boot <- lapply(
  fixef_names,
  FUN = function(y) {
    t(sapply(boot_result, function(x) {
      x$fixef[y, ]
    }))
  }
)


fixef_coef_covar_boot <- lapply(fixef_coef_samples_boot, var) 

fixef_fun_se_boot <- lapply(fixef_names, function(x) {
  sqrt(diag(Psi %*% fixef_coef_covar_boot[[x]] %*% t(Psi)))
})


sim_cb_bootstrap <- lapply(X = fixef_names, FUN = function(x) {
  mvn_sim(coef_point_est = fixef_coef_point[x, , drop = TRUE],
          coef_covar_mat = fixef_coef_covar_boot[[x]],
          Psi_basis = Psi,
          N_simulation_mvn = 10000, 
          coverage_level = 0.95)
})

sim_cb_wald <- lapply(X = fixef_names, FUN = function(x) {
  mvn_sim(coef_point_est = fixef_coef_point[x, , drop = TRUE],
          coef_covar_mat = diag(fixef_coef_var[x, , drop = TRUE]),
          Psi_basis = Psi,
          N_simulation_mvn = 10000, 
          coverage_level = 0.95)})

parameter_results_df <- purrr::map_dfr(fixef_names, function(x) {
  data.frame(
    dimension = rep(c("hip", "knee"), each = 101),
    t = rep(0:100, times = 2),
    point_est = fixef_fun_point[x,, drop = TRUE],
    se_boot = fixef_fun_se_boot[[x]],
    se_wald = fixef_fun_se[x, ],
    sim_boot_lower = sim_cb_bootstrap[[x]]$lower,
    sim_boot_upper = sim_cb_bootstrap[[x]]$upper,
    sim_wald_lower = sim_cb_wald[[x]]$lower,
    sim_wald_upper = sim_cb_wald[[x]]$upper
  )
}, .id = "parameter")

parameter_results_dt <- data.table(parameter_results_df)

parameter_results_dt[
  , `:=`(
    pw_wald_lower = point_est - 2 * se_wald,
    pw_wald_upper = point_est + 2 * se_wald,
    pw_boot_lower = point_est - 2 * se_boot,
    pw_boot_upper = point_est + 2 * se_boot
  )
]
levels(lme_dt$sex)
dput(unique(parameter_results_dt$parameter))
parameter_results_dt[,
  beta_label_part_1 := fcase(
    parameter == "(Intercept)", "Intercept",
    parameter == "risinjured_greater_than_2_yr", "RIS: Injured $>2$ yr.",
    parameter == "risinjured_1_to_2_yr", "RIS: Injured $1-2$ yr.",
    parameter == "risinjured_less_than_1_yr", "RIS: Injured $<1$ yr.",
    parameter == "sexfemale", "Sex",
    parameter == "speed_cent", "Speed",
    parameter == "weight_kg_cent", "Weight",
    parameter == "height_cm_cent", "Height",
    parameter == "age_cent", "Age"
  )
]
parameter_results_dt[, beta_label_part_1 := factor(
  beta_label_part_1, 
  levels = c(
    "Intercept",
    "Speed",
    "Sex",
    "RIS: Injured $>2$ yr.",
    "RIS: Injured $1-2$ yr.",
    "RIS: Injured $<1$ yr.",
    "Age",
    "Weight",
    "Height"
  ))]




hip_coef_functions_plot <- ggplot(data = parameter_results_dt[dimension == "hip"]) +
  aes(x = t, colour = beta_label_part_1) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_line(aes(y = point_est)) +
  geom_ribbon(aes(ymin = sim_boot_lower, ymax = sim_boot_upper, fill = beta_label_part_1), alpha = 0.25) +
  geom_line(aes(y = pw_boot_lower), lty = 2) +
  geom_line(aes(y = pw_boot_upper), lty = 2) +
  geom_hline(yintercept  = 0, col = "darkgrey") +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip")

knee_coef_functions_plot <- ggplot(data = parameter_results_dt[dimension == "knee"]) +
  aes(x = t, colour = beta_label_part_1) +
  geom_hline(yintercept  = 0, col = "grey") +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_line(aes(y = point_est)) +
  geom_ribbon(aes(ymin = sim_boot_lower, ymax = sim_boot_upper, fill = beta_label_part_1), alpha = 0.25) +
  geom_line(aes(y = pw_boot_lower), lty = 2) +
  geom_line(aes(y = pw_boot_upper), lty = 2) +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(knee)}_a (t)$",
       title = "Knee")

hip_coef_functions_plot
knee_coef_functions_plot

tikz(file.path(plots_path, "coefficient-functions-plot.tex"),
     width = 1 * doc_width_inches, 
     height = 1.5 * (doc_width_inches))
ggarrange(hip_coef_functions_plot, knee_coef_functions_plot, nrow =2, ncol = 1)
dev.off()




# -------------------------------------------------------------------------

ggplot(data = parameter_results_dt[dimension == "hip"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_line(aes(y = point_est)) +
  geom_ribbon(aes(ymin = sim_boot_lower, ymax = sim_boot_upper), fill = "red", alpha = 0.25) +
  geom_ribbon(aes(ymin = sim_wald_lower, ymax = sim_wald_upper), fill = "blue", alpha = 0.25)  + 
  geom_line(aes(y = pw_boot_lower), lty = 2) +
  geom_line(aes(y = pw_boot_upper), lty = 2) +
  geom_line(aes(y = pw_wald_lower), lty = 3) +
  geom_line(aes(y = pw_wald_upper), lty = 3) +
  geom_hline(yintercept  = 0, col = "darkgrey") +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip")

ggplot(data = parameter_results_dt[dimension == "knee"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_line(aes(y = point_est)) +
  geom_ribbon(aes(ymin = sim_boot_lower, ymax = sim_boot_upper), fill = "grey", alpha = 0.25, col = 1) +
  # geom_line(aes(y = pw_boot_lower), lty = 2) +
  # geom_line(aes(y = pw_boot_upper), lty = 2) +
  # geom_line(aes(y = pw_wald_lower), lty = 3) +
  # geom_line(aes(y = pw_wald_upper), lty = 3) +
  geom_hline(yintercept  = 0, col = "darkgrey") +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Knee")



# -------------------------------------------------------------------------

# SELF-SELECTED SPEEDS PLOT -----------------------------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
mean_speed <- mean(basis_transformation_results$coef_dts[[1]]$self_selected_speed_kmph)

intercept_dt <- parameter_results_dt[beta_label_part_1 == "Intercept",
                                     .(t, dimension, intercept = point_est)
]
speed_effect_dt <- parameter_results_dt[beta_label_part_1 == "Speed",
                                        .(t, dimension, speed_effect = point_est)]

speed_predictions_dt <- merge.data.table(x = intercept_dt, y = speed_effect_dt,by = c(
  "t", "dimension"
))

speed_predictions_dt[,
                     kmph_11 := intercept + (11 - mean_speed) * speed_effect
]

speed_predictions_dt[,
                     kmph_11 := intercept + (11 - mean_speed) * speed_effect
]
speed_predictions_dt[, `:=`(
  kmph_12 = kmph_11 + 1 * speed_effect,
  kmph_10 = kmph_11 - 1 * speed_effect,
  kmph_13 = kmph_11 + 2 * speed_effect,
  kmph_9 = kmph_11 - 2 * speed_effect
)]
speed_predictions_dt[, `:=`(
  intercept = NULL, speed_effect = NULL
)]
speed_predictions_dt_long <- melt.data.table(data = speed_predictions_dt,
                                             id.vars = c("t", "dimension"),
                                             variable.name = "speed",
                                             value.name = "angle")
speed_predictions_dt_long[, speed := as.numeric(stringr::str_remove(speed, "kmph_"))]
speed_predictions_dt_long[, speed := factor(ordered(speed), ordered = FALSE)]

# Create plots: -----------------------------------------------------------
p1_hip <- ggplot(data = speed_predictions_dt_long[dimension == "hip"]) +
  aes(x=t, y = angle, group = speed, colour = speed) +
  geom_line() +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Angle ($^{\\circ}$)",
       colour = "Speed (kmph)") +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  labs(title = "(a) Hip Angle")

p1_knee <- ggplot(data = speed_predictions_dt_long[dimension == "knee"]) +
  aes(x=t, y = angle, group = speed, colour = speed) +
  geom_line() +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Angle ($^{\\circ}$)",
       colour = "Speed (kmph)") +
  guides(colour = guide_legend(override.aes = list(size = 1))) +
  labs(title = "(b) Knee Angle")


speed_dt_wide <- dcast(speed_predictions_dt_long, formula = ... ~ dimension, value.var = "angle")
p2<-ggplot(data = speed_dt_wide) +
  aes(x=hip, y = knee, group = speed, colour = speed) +
  geom_path() +
  geom_label(data = speed_dt_wide[speed == 11 & t %in% c(0, 25, 50, 75, 100)],
             mapping = aes(label = paste0("$", t, "\\%$")), size = 2.25, col = 1) +
  labs(x = "Hip Angle ($^{\\circ}$)",
       y = "Knee Angle ($^{\\circ}$)",
       colour = "Speed (kmph)",
       title = "(c) Angle-Angle Diagram") +
  guides(colour = guide_legend(override.aes = list(size = 1)))
p2

barplot_dt <- basis_transformation_results$coef_dts[[1]][, .(self_selected_speed_kmph)]
barplot_dt[, fill_var := fifelse(self_selected_speed_kmph %in% 9:13, self_selected_speed_kmph, 0)]
barplot_dt[, fill_var := factor(fill_var)]
p3 <- ggplot(data = barplot_dt) +
  aes(x = self_selected_speed_kmph, fill = fill_var) +
  geom_bar(col = "black", linewidth = 0.25) +
  scale_fill_manual(values = c("grey", gg_color_hue(5))) +
  labs(x = "Speed (kmph)", y = "Count", title = "Barplot of Observed Speeds") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 145)) +
  scale_x_continuous(breaks = c(6:16)) + 
  theme(legend.position = "none") +
  ggtitle("(d) Observed Speeds")



# Save for Publishing: ---------------------------------------------------
ggarrange(p1_hip, p1_knee, p2, p3, nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom")
tikz(file.path(plots_path, "speed-predictions.tex"),
     width = 1.5 * doc_width_inches, 
     height = 0.45 * (doc_width_inches))
ggarrange(p1_hip, p1_knee, p2, p3, nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom")
dev.off()

