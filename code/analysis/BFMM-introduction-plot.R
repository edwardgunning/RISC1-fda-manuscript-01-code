# ------------------------------------------------------------------------#
# Check data with some graphical summaries and include them in
# summary write up of the data.
# ------------------------------------------------------------------------#

# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(ggplot2)    # CRAN v3.3.5
library(fda)        # CRAN v5.5.1
library(tikzDevice) # CRAN v0.12.3.1
library(ggrepel)

# Settings for ggplot() figures: ------------------------------------------
source(here::here("code", "functions", "theme_gunning.R"))
theme_gunning()
theme_update(legend.title = element_text(size = 9, hjust = 0.5))

# rough guide for sizing of plot outputs:
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Load some helper functions for functional data manipulation: ------------
source(here::here("code", "functions", "functions-helper-smoothing.R"))


# -------------------------------------------------------------------------
# Path to save the ouputs of analysis: ------------------------------------
plots_path <- here::here("outputs", "figures")
results_path <- here::here("outputs", "results")
data_path <- here::here("data")

# Read in data ------------------------------------------------------------
# Simplified dataset stored as an rds object.
# functional data is represented in terms of 80 basis coefficients
subject_side_coef_reg <- readRDS(file.path(data_path, "subject_side_coef_reg.rds"))


# Define basis to combine with stored coefficients
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)


# Evaluate functional data on grid: ---------------------------------------
subject_side_coef_reg[, {
  coef_mat <- coef_to_mat(dt = .SD)
  fd_obj <- fd(coef = coef_mat, basisobj = bspl80)
  plot(fd_obj)
  title(paste0("Registered ", location))},
  by = location,
  .SDcols = paste0("lm_coef_", 1:80)]

subject_side_eval_reg <-copy(subject_side_coef_reg)
subject_side_eval_reg[, paste(0:100) :=
                        {
                          coef_matrix <- as.matrix(t(.SD))
                          fd_obj <- fda::fd(coef = coef_matrix,
                                            basisobj = bspl80)
                          eval_matrix <- fda::eval.fd(evalarg = 0:100,
                                                      fdobj = fd_obj)
                          get_list_of_rows(eval_matrix)
                        },
                      .SDcols = paste("lm_coef", 1:80, sep = "_")]
subject_side_eval_reg[, paste("lm_coef", 1:80, sep = "_") := NULL]

subject_side_eval_reg <- subject_side_eval_reg[, c("subject_id", "side", "location", paste(0:100))]
# -------------------------------------------------------------------------
eval_reg_lng <- melt.data.table(subject_side_eval_reg,
                                measure.vars = paste0(0:100),
                                variable.name = "t",
                                variable.factor = FALSE,
                                value.name = "angle",
                                verbose = TRUE)
eval_reg_lng[, t := as.numeric(t)]
eval_reg_wider <- dcast.data.table(data =  eval_reg_lng,
                                   formula =  subject_id + side + t  ~ location,
                                   value.var = "angle")
eval_reg_wider_mean <- eval_reg_wider[, .(Hip = mean(Hip), Knee = mean(Knee)), by = .(t)]

knee_plot <- ggplot(data = eval_reg_wider) + 
  aes(x = t, y = Knee, group = interaction(subject_id, side)) +
  geom_line(colour = "grey15", alpha = 0.5) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Knee Angle ($^{\\circ}$)") +
  geom_line(data = eval_reg_wider[subject_id== "P_4001"],
            aes(colour = side), lwd = 1) +
  theme(legend.position = "none")

hip_plot <- ggplot(data = eval_reg_wider) + 
  aes(x = t, y = Hip, group = interaction(subject_id, side)) +
  geom_line(colour = "grey15", alpha = 0.5) +
  labs(x = "Normalised Time ($\\%$ of Stride)",
       y = "Hip Angle ($^{\\circ}$)") +
  geom_line(data = eval_reg_wider[subject_id== "P_4001"],
            aes(colour = side), lwd = 1) +
  theme(legend.position = "none")

hip_knee_plot <- ggplot(data = eval_reg_wider) + 
  aes(x = Hip, y = Knee, group = interaction(subject_id, side)) +
  geom_path(colour = "grey15", alpha = 0.5) +
  labs(x = "Hip Angle ($^{\\circ}$)",
       y = "Knee Angle ($^{\\circ}$)") +
  geom_path(data = eval_reg_wider[subject_id== "P_4001"],
            aes(colour = side), lwd = 1) +
  geom_label_repel(data = eval_reg_wider_mean[t %in% c(0, 100)], inherit.aes = FALSE,
             mapping = aes(x = Hip, y = Knee, label = paste0("$", t, "\\%$")),
             size = 2, 
             col = 1, 
             segment.color = "white",
             min.segment.length = 0.1) +
  geom_label(data = eval_reg_wider_mean[t %in% c(25, 50, 75)], inherit.aes = FALSE,
                   mapping = aes(x = Hip, y = Knee, label = paste0("$", t, "\\%$")),
                   size = 2) +
  geom_point(data = eval_reg_wider_mean[t %in% c(0, 100)], inherit.aes = FALSE, size = 0.1,
             mapping = aes(x = Hip, y = Knee), col = "white") +
  theme(legend.position = "none")



full_plot <- ggpubr::ggarrange(plotlist = list(hip_plot, knee_plot, hip_knee_plot),
                               nrow = 1, ncol = 3, labels = c("(a)", "(b)", "(c)"),
                               font.label = list(size = 11, color = "black", face = "plain", family = NULL),
                               label.x = 0.2, label.y = 0.95)

full_plot

tikz(file.path(plots_path, "intro-plot.tex"),
     width = doc_width_inches, 
     height = doc_width_inches/3)
print(full_plot)
dev.off()




# Extra - American spelling for JASA --------------------------------------

knee_plot_Am <- ggplot(data = eval_reg_wider) + 
  aes(x = t, y = Knee, group = interaction(subject_id, side)) +
  geom_line(colour = "grey15", alpha = 0.5) +
  labs(x = "Normalized Time ($\\%$ of Stride)",
       y = "Knee Angle ($^{\\circ}$)") +
  geom_line(data = eval_reg_wider[subject_id== "P_4001"],
            aes(colour = side), lwd = 1) +
  theme(legend.position = "none")

hip_plot_Am <- ggplot(data = eval_reg_wider) + 
  aes(x = t, y = Hip, group = interaction(subject_id, side)) +
  geom_line(colour = "grey15", alpha = 0.5) +
  labs(x = "Normalized Time ($\\%$ of Stride)",
       y = "Hip Angle ($^{\\circ}$)") +
  geom_line(data = eval_reg_wider[subject_id== "P_4001"],
            aes(colour = side), lwd = 1) +
  theme(legend.position = "none")


full_plot_Am <- ggpubr::ggarrange(plotlist = list(hip_plot_Am, knee_plot_Am, hip_knee_plot),
                               nrow = 1, ncol = 3, labels = c("(a)", "(b)", "(c)"),
                               font.label = list(size = 11, color = "black", face = "plain", family = NULL),
                               label.x = 0.2, label.y = 0.95)


tikz(file.path(plots_path, "intro-plot-Am.tex"),
     width = doc_width_inches, 
     height = doc_width_inches/3)
print(full_plot_Am)
dev.off()
