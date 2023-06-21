library(ggplot2)    # CRAN v3.4.0
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(xtable)     # CRAN v1.8-4
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
mfamm_results <- readRDS(here::here(results_path, "mfamm-comparison.rds"))
source(here::here("code/functions/theme_gunning.R"))

results_list <- mfamm_results$results_list
time_vec <- mfamm_results$time_vec
settings <- mfamm_results$settings

plot(interaction(settings$bf_covs, settings$mfpc_cutoff), time_vec/60)

                    # Create table of computation times: --------------------------------------
time_df <- data.frame(K = integer(settings$bf_covs), 
           PVE = settings$mfpc_cutoff,
           ct = round(time_vec/60,2))
names(time_df) <- c("$K$ Marginal", "PVE",  "Time (mins)")
time_df
bold <- function(x) {
  paste0("{\\bfseries ", x, "}") 
}
time_table <- xtable(time_df, 
                     digits =  c(0,0, 2, 2), 
                     label = "tab:mfamm-comp-time",
                     caption = "Computation time for the multiFAMM model with different settings.")
align(time_table)[1] <- "l"
print(time_table, 
      file = here::here("outputs", "tables", "mfamm-comp-time.tex"),
      sanitize.text.function = function(x){x},
      sanitize.colnames.function = bold,
      include.rownames = FALSE,
      booktabs = TRUE)

                    # Create plot of results: -------------------------------------------------
settings_index <- seq_len(8)
names(settings_index) <- paste0("setting_", settings_index)
plot_df <- purrr::map_dfr(.x = settings_index, .f = function(x) {
  df_i <- results_list[[x]]
  df_i$bf_covs <- settings[x, "bf_covs"]
  df_i$cutoff <- settings[x, "mfpc_cutoff"]
  df_i
})
plot_dt <- as.data.table(plot_df)

# -------------------------------------------------------------------------
bootstrap_results <- readRDS(file.path(results_path, "bootstrap-results.rds"))
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

plot_dt[, combo := paste0("$K =", bf_covs, "$, pve $=", cutoff, "$")]
plot_dt[, combo := factor(combo,
                          levels = c(paste0("$K =", 5, "$, pve $=", c(0.9, 0.95), "$"), 
                                     paste0("$K =", 8, "$, pve $=", c(0.9, 0.95), "$"),
                                     paste0("$K =", 10, "$, pve $=", c(0.9, 0.95), "$"),
                                     paste0("$K =", 15, "$, pve $=", c(0.9, 0.95), "$")))]




hip <-ggplot(data = plot_dt[dim == "hip"]) +
  aes(x = t, y = estimate, colour = combo) +
  geom_line(alpha = 1, aes(linetype = "Point Estimate")) +
  geom_hline(yintercept = 0, col = "grey") +
  facet_wrap(~ beta_label_part_1, scales = "free_y")  +
  geom_line(data = parameter_results_dt[dimension=="hip"],
            aes(x = t, y = point_est),
            inherit.aes = FALSE,
            linewidth = 0.8) +
  scale_linetype_manual(values = c(1, 3)) +
  geom_line(aes(y = lower, linetype = "Pointwise CI")) +
  geom_line(aes(y = upper, linetype = "Pointwise CI")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box = "vertical") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(hip)}_a (t)$",
       title = "Hip") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.8, colour = "darkslategrey")),
         colour = guide_legend(override.aes = list(linewidth = 0.8)))

knee <-ggplot(data = plot_dt[dim == "knee"]) +
  aes(x = t, y = estimate, colour = combo) +
  geom_line(alpha = 1, aes(linetype = "Point Estimate")) +
  geom_hline(yintercept = 0, col = "grey") +
  facet_wrap(~ beta_label_part_1, scales = "free_y")  +
  geom_line(data = parameter_results_dt[dimension=="knee"],
            aes(x = t, y = point_est),
            inherit.aes = FALSE,
            linewidth = 0.8) +
  scale_linetype_manual(values = c(1, 3)) +
  geom_line(aes(y = lower, linetype = "Pointwise CI")) +
  geom_line(aes(y = upper, linetype = "Pointwise CI")) +
  theme(legend.title = element_blank(),
        legend.position = "bottom", 
        legend.box = "vertical") +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Coefficient Function $\\beta^{(knee)}_a (t)$",
       title = "Knee") +
  guides(linetype = guide_legend(override.aes = list(linewidth = 0.8, colour = "darkslategrey")),
         colour = guide_legend(override.aes = list(linewidth = 0.8))) 
  
 
combined_plot <- ggpubr::ggarrange(hip, knee, common.legend = TRUE, nrow = 2, legend = "bottom")
combined_plot   

tikz(file.path(plots_path, "mfamm-comparison.tex"),
     width = 1 * doc_width_inches, 
     height = 1.65 * (doc_width_inches),
     standAlone = TRUE)
combined_plot
dev.off()

tinytex::lualatex(file.path(plots_path, "mfamm-comparison.tex"))


