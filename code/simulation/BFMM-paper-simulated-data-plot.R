# -------------------------------------------------------------------------
# plot of sample of simulated data for paper
# -------------------------------------------------------------------------
# Load packages: ----------------------------------------------------------
library(data.table) # CRAN v1.14.2 # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.0 # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1
# -------------------------------------------------------------------------
functions_path <- here::here("code", "functions")
results_path <- here::here("outputs", "results")
plots_path <- here::here("outputs", "figures")

# Load all helper functions: ----------------------------------------------
source(here::here("code", "simulation", "BFMM-paper-generate-simulated-data.R"))
source(file.path(functions_path, "theme_gunning.R"))


# GGplot theme: -----------------------------------------------------------
theme_gunning()
theme_update(legend.title = element_text(size = 9, hjust = 0.5))
# rough sizing for figures:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# -------------------------------------------------------------------------
set.seed(1996)
sample_simdata_1 <- generate_data_scenario_1(N = 10, J = 2)
sample_simdata_2 <- generate_data_scenario_2(N = 10, J = 2)

df_1 <- sample_simdata_1$df
Y_df_1 <- data.frame(sample_simdata_1$Y)
names(Y_df_1) <- c(paste0("Hip_", 0:100), paste0("Knee_", 0:100))
plotting_df_1 <- cbind(df_1, Y_df_1)

df_2 <- sample_simdata_2$df
Y_df_2 <- data.frame(sample_simdata_2$Y)
names(Y_df_2) <- c(paste0("Hip_", 0:100), paste0("Knee_", 0:100))
plotting_df_2 <- cbind(df_2, Y_df_2)

# use data.table to reshape 
plotting_dt_1 <- data.table(plotting_df_1)
plotting_dt_1_lng <- melt.data.table(data = plotting_dt_1,
                                   id.vars = names(df_1),
                                   variable.name = "dim_time",
                                   value.name = "y",
                                   variable.factor = FALSE, 
                                   measure.vars = names(Y_df_1),
                                   value.factor = FALSE,
                                   verbose = TRUE)

plotting_dt_1_lng[, dim := stringr::str_remove(dim_time, "_\\d{1,3}")]
plotting_dt_1_lng[, t := as.numeric(stringr::str_extract(dim_time, "\\d{1,3}"))]

plotting_dt_2 <- data.table(plotting_df_2)
plotting_dt_2_lng <- melt.data.table(data = plotting_dt_2,
                                     id.vars = names(df_2),
                                     variable.name = "dim_time",
                                     value.name = "y",
                                     variable.factor = FALSE, 
                                     measure.vars = names(Y_df_2),
                                     value.factor = FALSE,
                                     verbose = TRUE)

plotting_dt_2_lng[, dim := stringr::str_remove(dim_time, "_\\d{1,3}")]
plotting_dt_2_lng[, t := as.numeric(stringr::str_extract(dim_time, "\\d{1,3}"))]




p1 <- ggplot(plotting_dt_1_lng[subject_id %in% 1:5]) +
  aes(x = t, y = y, group = interaction(subject_id, side),
      colour = subject_id) +
  facet_wrap( ~  dim, scales = "free_y") +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_line() + 
  labs(colour = "Subject ID:", y = "$y(t)$", x = "$t$", title = "Scenario 1") 
p1

p2 <- ggplot(plotting_dt_2_lng[subject_id %in% 1:5]) +
  aes(x = t, y = y, group = interaction(subject_id, side),
      colour = subject_id) +
  facet_wrap( ~  dim, scales = "free_y") +
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_line() + 
  labs(colour = "Subject ID:", y = "$y(t)$", x = "$t$", title = "Scenario 2") 
p2

p <- ggpubr::ggarrange(plotlist = list(p1, p2), ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")

p


# -------------------------------------------------------------------------
tikz(file.path(plots_path, "simulated-dataset.tex"),
     width =  doc_width_inches, 
     height = 1.5 * (doc_width_inches/4))
print(p)
dev.off()

