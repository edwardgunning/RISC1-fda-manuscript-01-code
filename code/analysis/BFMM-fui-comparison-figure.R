# -------------------------------------------------------------------------
# Figures to compare point estimates and standard errors to 
# fast univatiate inference methods by Cui (2022).
# -------------------------------------------------------------------------

# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0
library(ggplot2)    # CRAN v3.4.0


# -------------------------------------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning() # set theme
# + some customisations for plots
theme_update(legend.key.size = unit(0.95,"line"))

# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Path to save the outputs of analysis: -----------------------------------
plots_path <- here::here("outputs", "figures")
results_path <- here::here("outputs", "results")


# -------------------------------------------------------------------------
bootstrap_results <- readRDS(file.path(results_path, "bootstrap-results.rds"))
fui_results_df <- readRDS(here::here("outputs", "results", "fui_comparison.rds"))
fui_results_dt <- data.table(fui_results_df)


# Format data for plotting: ----------------------------------------------


parameter_results_dt <- bootstrap_results$parameter_results_dt
stopifnot(nrow(parameter_results_dt) == nrow(fui_results_dt))

combined_parameter_results_dt <- merge.data.table(x = parameter_results_dt,
                                                  y = fui_results_dt,
                                                  by = c("parameter", "dimension", "t"),
                                                  all = TRUE)
stopifnot(nrow(combined_parameter_results_dt) == nrow(parameter_results_dt))

combined_parameter_results_dt[,
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
combined_parameter_results_dt[, beta_label_part_1 := factor(
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




ggplot(data = combined_parameter_results_dt[dimension == "hip"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_line(aes(y = point_est, colour = "Current", linetype = "Point Estimate")) +
  geom_line(aes(y = fui_point_est, colour = "FUI", linetype = "Point Estimate")) +
  geom_ribbon(aes(ymin = pw_fui_lower_analytic, ymax =  pw_fui_upper_analytic, fill = "FUI"), alpha = 0.25) +
  geom_ribbon(aes(ymin = pw_wald_lower, ymax = pw_wald_upper, fill = "current"), alpha = 0.25)

?geom_ribbon
# -------------------------------------------------------------------------
head(combined_parameter_results_dt)
(a1 <- ggplot(data = combined_parameter_results_dt[dimension == "hip"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_line(aes(y = point_est, colour = "Current")) +
  geom_line(aes(y = fui_point_est, colour = "FUI")) +
  geom_ribbon(aes(ymin = pw_fui_lower_analytic, ymax =  pw_fui_upper_analytic, fill = "FUI"), alpha = 0.25) +
  geom_ribbon(aes(ymin = pw_wald_lower, ymax = pw_wald_upper, fill = "Current"), alpha = 0.25) +
  # scale_colour_manual(values = c("red", "blue")) +
  # scale_fill_manual(values = c("red", "blue")) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_linetype_manual(values = c(1, 2)) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1)),
         colour = guide_legend(override.aes = list(linewidth = 1))))

(a2 <- ggplot(data = combined_parameter_results_dt[dimension == "knee"]) +
    aes(x = t) +
    facet_wrap(~ beta_label_part_1, scales = "free_y") +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_line(aes(y = point_est, colour = "Current")) +
    geom_line(aes(y = fui_point_est, colour = "FUI")) +
    geom_ribbon(aes(ymin = pw_fui_lower_analytic, ymax =  pw_fui_upper_analytic, fill = "FUI"), alpha = 0.25) +
    geom_ribbon(aes(ymin = pw_wald_lower, ymax = pw_wald_upper, fill = "Current"), alpha = 0.25) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_linetype_manual(values = c(1, 2)) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Knee") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1)),
         colour = guide_legend(override.aes = list(linewidth = 1))))


tikz(file.path(plots_path, "fui-comparison-analytic.tex"),
     width = 1 * doc_width_inches, 
     height = 1.5 * (doc_width_inches),
     standAlone = TRUE)
ggarrange(a1, a2, common.legend = TRUE, nrow = 2, legend = "bottom")
dev.off()

tinytex::lualatex(file.path(plots_path, "fui-comparison-analytic.tex"))


(b1 <- ggplot(data = combined_parameter_results_dt[dimension == "hip"]) +
  aes(x = t) +
  facet_wrap(~ beta_label_part_1, scales = "free_y") +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_line(aes(y = point_est, colour = "Current")) +
  geom_line(aes(y = fui_point_est, colour = "FUI")) +
  geom_ribbon(aes(ymin = pw_fui_lower_boot, ymax =  pw_fui_upper_boot, fill = "FUI"), alpha = 0.25) +
  geom_ribbon(aes(ymin = pw_boot_lower, ymax = pw_boot_upper, fill = "Current"), alpha = 0.25) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_linetype_manual(values = c(1, 2)) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1)),
         colour = guide_legend(override.aes = list(linewidth = 1))))

(b2 <- ggplot(data = combined_parameter_results_dt[dimension == "knee"]) +
    aes(x = t) +
    facet_wrap(~ beta_label_part_1, scales = "free_y") +
    geom_hline(yintercept = 0, colour = "grey") +
    geom_line(aes(y = point_est, colour = "Current")) +
    geom_line(aes(y = fui_point_est, colour = "FUI")) +
    geom_ribbon(aes(ymin = pw_fui_lower_boot, ymax =  pw_fui_upper_boot, fill = "FUI"), alpha = 0.25) +
    geom_ribbon(aes(ymin = pw_boot_lower, ymax = pw_boot_upper, fill = "Current"), alpha = 0.25) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  scale_linetype_manual(values = c(1, 2)) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Knee") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1)),
         colour = guide_legend(override.aes = list(linewidth = 1))))

tikz(file.path(plots_path, "fui-comparison-bootstrap.tex"),
     width = 1 * doc_width_inches, 
     height = 1.5 * (doc_width_inches),
     standAlone = TRUE)
ggarrange(b1, b2, common.legend = TRUE, nrow = 2, legend = "bottom")
dev.off()

tinytex::lualatex(file.path(plots_path, "fui-comparison-bootstrap.tex"))



