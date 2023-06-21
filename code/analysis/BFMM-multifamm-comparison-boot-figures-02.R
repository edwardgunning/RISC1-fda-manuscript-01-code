library(ggplot2)    # CRAN v3.4.0
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
# -------------------------------------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning() # set theme
# + some customisations for plots
theme_update(legend.key.size = unit(0.95,"line"))

# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937


plots_path <- here::here("outputs", "figures")
results_path <- here::here("outputs", "results")
mfamm_results <- readRDS(here::here("outputs/results/mfamm-comparison-bootstrap-02.rds"))
bootstrap_results <- readRDS(file.path(results_path, "bootstrap-results.rds"))
parameter_results_dt <- bootstrap_results$parameter_results_dt
results_list <- mfamm_results$results_list
time_vec <- mfamm_results$time_vec

sum(time_vec)/60^2

boot_index <- seq_len(10)
names(boot_index) <- paste0(boot_index)
plot_df <- purrr::map_dfr(.x = boot_index, .f = function(x) {
  df_i <- results_list[[x]]
  df_i
}, .id = "b")
plot_dt <- as.data.table(plot_df)



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

plot_dt[, parameter := param]
plot_dt[,
        beta_label_part_1 := fcase(
          parameter == "intercept", "Intercept",
          parameter == "ris_effect[2]", "RIS: Injured $>2$ yr.",
          parameter == "ris_effect[3]", "RIS: Injured $1-2$ yr.",
          parameter == "ris_effect[4]", "RIS: Injured $<1$ yr.",
          parameter == "sex_effect[2]", "Sex",
          parameter == "speed_cent_coefficient", "Speed",
          parameter == "weight_cent", "Weight",
          parameter == "height_cent", "Height",
          parameter == "age_cent", "Age"
        )
]
plot_dt[, beta_label_part_1 := factor(
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

plot_dt[, b := factor(b, levels = paste(1:10))]
hip <-ggplot(data = plot_dt[dim == "hip"]) +
  aes(x = t, y = estimate, colour = b) +
  # geom_line(alpha = 1, aes(linetype = "Point Estimate")) +
  geom_line() +
  geom_hline(yintercept = 0, col = "darkgrey") +
  facet_wrap(~ beta_label_part_1, scales = "free_y")  +
  geom_line(data = parameter_results_dt[dimension=="hip"],
            aes(x = t, y = point_est),
            inherit.aes = FALSE,
            linewidth = 0.8) +
  # scale_linetype_manual(values = c(3, 1)) +
  # geom_line(aes(y = lower, linetype = "CI")) +
  # geom_line(aes(y = upper, linetype = "CI")) +
  theme(legend.position = "bottom", 
        legend.box = "vertical") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip",
       colour = "Bootstrap Replicate:") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1)),
         colour = guide_legend(override.aes = list(linewidth = 1)))

knee <-ggplot(data = plot_dt[dim == "knee"]) +
  aes(x = t, y = estimate, colour = b) +
  # geom_line(alpha = 1, aes(linetype = "Point Estimate")) +
  geom_line() +
  geom_hline(yintercept = 0, col = "darkgrey") +
  facet_wrap(~ beta_label_part_1, scales = "free_y")  +
  geom_line(data = parameter_results_dt[dimension=="knee"],
            aes(x = t, y = point_est),
            inherit.aes = FALSE,
            linewidth = 0.8) +
  # scale_linetype_manual(values = c(3, 1)) +
  # geom_line(aes(y = lower, linetype = "CI")) +
  # geom_line(aes(y = upper, linetype = "CI")) +
  theme(legend.position = "bottom", 
        legend.box = "vertical") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Knee",
       colour = "Bootstrap Replicate:") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 1)),
         colour = guide_legend(override.aes = list(linewidth = 1)))

combined_plot <- ggpubr::ggarrange(hip, knee, common.legend = TRUE, nrow = 2, legend = "bottom")
combined_plot

tikz(file.path(plots_path, "mfamm-comparison-bootstrap-02.tex"),
     width = 1 * doc_width_inches, 
     height = 1.65 * (doc_width_inches),
     standAlone = TRUE)
combined_plot
dev.off()

tinytex::lualatex(file.path(plots_path, "mfamm-comparison-bootstrap-02.tex"))
