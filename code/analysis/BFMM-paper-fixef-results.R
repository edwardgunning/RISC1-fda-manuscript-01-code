# -------------------------------------------------------------------------
# Plotting Fixed-Effects Results + Confidence Intervals
# and
# Plot of the effect of self-selected speed!
# -------------------------------------------------------------------------

# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0
library(fda)        # CRAN v5.5.1
library(ggplot2)    # CRAN v3.4.0
library(lme4)       # CRAN v1.1-30


# -------------------------------------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning() # set theme
# + some customisations for plots
theme_update(panel.grid.major = element_blank(),
             legend.key.size = unit(0.95,"line"))

# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Path to save the outputs of analysis: -----------------------------------
plots_path <- here::here("outputs", "figures")
results_path <- here::here("outputs", "results")


# -------------------------------------------------------------------------

basis_transformation_results <- readRDS(
  file.path(results_path, "basis-transform-results.rds"))
bootstrap_results <- readRDS(file.path(results_path, "bootstrap-results.rds"))



# Format data for plotting: ----------------------------------------------


parameter_results_dt <- bootstrap_results$parameter_results_dt
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
  geom_line(aes(y = pw_boot_lower),col = "red4") +
  geom_line(aes(y = pw_boot_upper),  col = "red4") +
  geom_line(aes(y = pw_wald_lower),  col = "blue4") +
  geom_line(aes(y = pw_wald_upper), col = "blue4") +
  geom_hline(yintercept  = 0, col = "darkgrey") +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip")

ggplot(data = parameter_results_dt[dimension == "knee"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_line(aes(y = point_est)) +
  geom_line(aes(y = pw_boot_lower),col = "red4") +
  geom_line(aes(y = pw_boot_upper),  col = "red4") +
  geom_line(aes(y = pw_wald_lower),  col = "blue4") +
  geom_line(aes(y = pw_wald_upper), col = "blue4") +
  geom_hline(yintercept  = 0, col = "darkgrey") +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Knee")

ggplot(data = parameter_results_dt[dimension == "hip"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_line(aes(y = point_est)) +
  geom_line(aes(y = sim_boot_lower),col = "red4") +
  geom_line(aes(y = sim_boot_upper),  col = "red4") +
  geom_line(aes(y = sim_wald_lower),  col = "blue4") +
  geom_line(aes(y = sim_wald_upper), col = "blue4") +
  geom_hline(yintercept  = 0, col = "darkgrey") +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip")

ggplot(data = parameter_results_dt[dimension == "knee"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_line(aes(y = point_est)) +
  geom_line(aes(y = sim_boot_lower),col = "red4") +
  geom_line(aes(y = sim_boot_upper),  col = "red4") +
  geom_line(aes(y = sim_wald_lower),  col = "blue4") +
  geom_line(aes(y = sim_wald_upper), col = "blue4") +
  geom_hline(yintercept  = 0, col = "darkgrey") +
  theme(legend.position = "none") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Knee")




# Plot of the speed effect: -----------------------------------------------
# get colours:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# What is the mean speed, i.e., the value of speed at baseline (intercept)
mean_speed <- mean(basis_transformation_results$coef_dts[[1]]$self_selected_speed_kmph)

# Extract intercept function and speed fumnction
intercept_dt <- parameter_results_dt[beta_label_part_1 == "Intercept",.(t, dimension, intercept = point_est)]
speed_effect_dt <- parameter_results_dt[beta_label_part_1 == "Speed", .(t, dimension, speed_effect = point_est)]

# create dataset to score speed predictions
speed_predictions_dt <- merge.data.table(
  x = intercept_dt,
  y = speed_effect_dt,
  by = c( "t", "dimension")
  )

speed_predictions_dt[,kmph_11 := intercept + (11 - mean_speed) * speed_effect]

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

(p_speed <- ggarrange(p1_hip, p1_knee, p2, p3, 
                      nrow = 1, ncol = 4,
                      common.legend = TRUE,
                      legend = "bottom"))

# Save for Publishing: ---------------------------------------------------

tikz(file.path(plots_path, "speed-predictions.tex"),
     width = 1.5 * doc_width_inches, 
     height = 0.45 * (doc_width_inches))
ggarrange(p1_hip, p1_knee, p2, p3, nrow = 1, ncol = 4, common.legend = TRUE, legend = "bottom")
dev.off()

