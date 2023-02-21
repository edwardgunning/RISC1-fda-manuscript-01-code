# ------------------------------------------------------------------------#
# Extra plot , plotting the pointwise confidence interval coverage along 
# the functional domain.
# ------------------------------------------------------------------------#


# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.3.5
library(ggpubr)     # CRAN v0.4.0
library(tikzDevice) # CRAN v0.12.3.1


# -------------------------------------------------------------------------
# File paths: -------------------------------------------------------------
functions_path <- here::here("code", "functions")
results_path <- here::here("outputs", "results")
plots_path <- here::here("outputs", "figures")

# Source scripts: ---------------------------------------------------------
source(file.path(functions_path, "binomial_se.R"))
source(file.path(functions_path, "theme_gunning.R"))

# Read in Simulation Results ----------------------------------------------
results_list <- readRDS(
  file.path(results_path, "BFMM-tidied-simulation-results.rds")
)


# Extract Coverage Results: -----------------------------------------------
coverage_results <- results_list$coverage
cover_boot_pw_array <- coverage_results$cover_boot_pw_array
cover_wald_pw_array <- coverage_results$cover_wald_pw_array


# Rough plots: ------------------------------------------------------------
par(mfrow = c(1, 2))
matplot(apply(cover_wald_pw_array[,, 1:500],c(1,2), mean), type = "l")
matlines(apply(cover_boot_pw_array[,, 1:500],c(1,2), mean), type = "l", col = 3:4)
matplot(apply(cover_wald_pw_array[,, 501:1000],c(1,2), mean), type = "l")
matlines(apply(cover_boot_pw_array[,, 501:1000],c(1,2), mean), type = "l", col = 3:4)


# And look at range of values coverage takes on: --------------------------
range(apply(cover_wald_pw_array[,, 501:1000],c(1,2), mean))
range(apply(cover_boot_pw_array[,, 501:1000],c(1,2), mean))
range(apply(cover_wald_pw_array[,, 1:500],c(1,2), mean))
range(apply(cover_boot_pw_array[,, 1:500],c(1,2), mean))



# EGt data into format for publishable plot: ------------------------------
coverage_plot_dt <- data.table(
  scenario = factor(rep(c(1, 2), each = 404)),
  method = rep(rep(c("Wald", "Bootstrap"), each = 202), times = 2),
  dim = rep(rep(paste0("dim_", 1:2), each = 101), times = 4),
  t = rep(0:100, times = 8),
  rbind(apply(cover_wald_pw_array[,, 1:500],c(1,2), mean),
        apply(cover_boot_pw_array[,, 1:500],c(1,2), mean),
        apply(cover_wald_pw_array[,, 501:1000],c(1,2), mean),
        apply(cover_boot_pw_array[,, 501:1000],c(1,2), mean)))
coverage_plot_dt_lng<- melt.data.table(coverage_plot_dt, 
                                       id.vars = c("scenario", "method", "dim", "t"),
                                       measure.vars = c("(Intercept)", "sexfemale", "speed"),
                                       variable.name = "parameter",
                                       value.name = "coverage",
                                       variable.factor = FALSE,
                                       value.factor = FALSE
                                        )

coverage_plot_dt_lng[, se := binomial_se(coverage, n = 500)]
coverage_plot_dt_lng[, `:=`(
  lower = coverage - 2 * se,
  upper = coverage + 2 * se
)]


# want to fix y axis limits so they're same on both plots:
ylim_range <- range(coverage_plot_dt_lng[, .(upper, lower)])
coverage_plot_dt_lng[, param_name := 
                       fcase(
                         parameter == "(Intercept)", "$\\beta_0 (t)$",
                         parameter == "sexfemale", "$\\beta_1 (t)$",
                         parameter == "speed", "$\\beta_2 (t)$"
                       )]
coverage_plot_dt_lng[, dim := stringr::str_remove(dim, "dim_")]
coverage_plot_dt_lng[, dim := paste("Dimension", dim)]


# Set theme for plot: -----------------------------------------------------
theme_set(theme_bw())
theme_update(strip.text = element_text(size = 9),
             text = element_text(size = 9),
             panel.grid.minor = element_blank(),
             axis.title = element_text(size = 9),
             legend.text = element_text(size = 9, hjust = 0.5),
             legend.title = element_text(size = 9, hjust = 0.5, face = "bold"),
             plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))


# Plot --------------------------------------------------------------------
s1 <- ggplot(data = coverage_plot_dt_lng[scenario == 1]) +
  aes(x = t,
      y = coverage,
      colour = method,
      group = interaction(method, scenario)) +
  facet_wrap(dim ~ param_name, nrow = 2, ncol = 3) +
  geom_line() +
  ylim(ylim_range) +
  labs(y = "Coverage", x = "$t$",
       title = "Scenario 1",
       colour = "Method:",
       fill = "Method:") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.1)+
  geom_hline(yintercept = 0.95)


s2 <- ggplot(data = coverage_plot_dt_lng[scenario == 2]) +
  aes(x = t,
      y = coverage,
      colour = method,
      group = interaction(method, scenario)) +
  facet_wrap(dim ~ param_name, nrow = 2, ncol = 3) +
  geom_line() +
  ylim(ylim_range) +
  labs(y = "Coverage", x = "$t$",
       title = "Scenario 2",
       colour = "Method:",
       fill = "Method:") +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = method), alpha = 0.1)+
  geom_hline(yintercept = 0.95)



(coverage_plot <- ggarrange(plotlist = list(s1, s2),
                            ncol = 2,
                            nrow = 1,
                            legend = "top",
                            common.legend = TRUE))


# Save: ------------------------------------------------------------------

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

tikz(file.path(plots_path, "simulation-coverage-plot.tex"),
     width = doc_width_inches, 
     height = (0.65 * doc_width_inches))
print(coverage_plot)
dev.off()


# And print range of coverage: --------------------------------------------
coverage_range_summary <- coverage_plot_dt_lng[, .(min = round(min(coverage), 2),
                         max = round(max(coverage), 2)),
                     by = .(method, scenario)]
fwrite(x = coverage_range_summary, file = file.path(results_path, "BFMM-paper-coverage-pw-ranges.csv"))
