# -------------------------------------------------------------------------
# Assess simulation coverage and save coverage tables:
# -------------------------------------------------------------------------
library(data.table) # CRAN v1.14.2

# File paths: -------------------------------------------------------------
functions_path <- here::here("code", "functions")
results_path <- here::here("outputs", "results")
plots_path <- here::here("outputs", "figures")


# Load results: -----------------------------------------------------------
results_list <- readRDS(
  file = file.path(results_path, "BFMM-tidied-simulation-results.rds")
)

# load helper function ----------------------------------------------------
source(file.path(functions_path, "binomial_se.R"))

# Extract coverage: -------------------------------------------------------
coverage_results <- results_list$coverage
cover_boot_pw_array <- coverage_results$cover_boot_pw_array
cover_wald_pw_array <- coverage_results$cover_wald_pw_array
# and settings:
settings <- results_list$settings$settings




# Pointwise coverage: -----------------------------------------------------
# create data frame:
pw_dt <- data.table(
  rep = settings$rep,
  scenario = settings$scenario,
  t(apply(cover_wald_pw_array, 2:3, mean)),
  t(apply(cover_boot_pw_array, 2:3, mean))
)

names(pw_dt)[3:4] <- paste0("wald_", names(pw_dt)[3:4])
names(pw_dt)[5:6] <- paste0("boot_", names(pw_dt)[5:6])

# Reshape
pw_dt_lng <- melt.data.table(data = pw_dt,
                             id.vars = c("rep", "scenario"), 
                             measure.vars = c("wald_sexfemale", 
                                              "wald_speed",
                                              "boot_sexfemale", 
                                              "boot_speed"),
                             variable.name = "method", 
                             value.name = "covered",
                             variable.factor = FALSE, 
                             value.factor = FALSE,
                             verbose = TRUE
                             )


pw_dt_lng[, parameter := fifelse(stringr::str_detect(string = method, pattern =  "sexfemale"), yes = "sexfemale", no = "speed")]
pw_dt_lng[, method := fifelse(stringr::str_detect(string = method, pattern =  "wald"), yes = "wald", no = "boot")]

# Calculate average coverage:
pw_summary <- pw_dt_lng[, .(coverage = mean(covered)), by = .(method, parameter, scenario)]
# Calculate Monte-Carlo standard error:
pw_summary[, se := binomial_se(phat = coverage, n = 500)]

pw_summary[, method := factor(method, levels = c("wald", "boot"))]

pw_summary_round <- pw_summary[, .(method, parameter, scenario, coverage = round(coverage, 2), se = round(se, 2))]
setorderv(pw_summary_round, c("method", "scenario", "parameter"))
pw_summary_round

# save:
fwrite(x = pw_summary_round, 
       file = file.path(results_path, "BFMM-paper-coverage_pw.csv"))




# Simultaneous ------------------------------------------------------------
cover_boot_sim_mat <- coverage_results$cover_boot_sim_array
cover_wald_sim_mat <- coverage_results$cover_wald_sim_array

apply(apply(cover_boot_sim_mat, c(2, 3), all)[,501:1000], 1, mean)
apply(apply(cover_wald_sim_mat, c(2, 3), all)[,501:1000], 1, mean)

apply(apply(cover_boot_sim_mat, c(2, 3), all)[,1:500], 1, mean)
apply(apply(cover_wald_sim_mat, c(2, 3), all)[,1:500], 1, mean)

sim_dt <- data.table(
  rep = settings$rep,
  scenario = settings$scenario,
  t(apply(cover_wald_sim_mat, c(2:3), all)),
  t(apply(cover_boot_sim_mat, c(2:3), all))
)

names(sim_dt)[3:4] <- paste0("wald_", names(sim_dt)[3:4])
names(sim_dt)[5:6] <- paste0("boot_", names(sim_dt)[5:6])
sim_dt_lng <- melt.data.table(data = sim_dt,
                             id.vars = c("rep", "scenario"), 
                             measure.vars = c("wald_sexfemale",
                                              "wald_speed", 
                                              "boot_sexfemale",
                                              "boot_speed"),
                             variable.name = "method", 
                             value.name = "covered",
                             variable.factor = FALSE, 
                             value.factor = FALSE,
                             verbose = TRUE
)
sim_dt_lng[, parameter := fifelse(stringr::str_detect(string = method, pattern =  "sexfemale"), yes = "sexfemale", no = "speed")]
sim_dt_lng[, method := fifelse(stringr::str_detect(string = method, pattern =  "wald"), yes = "wald", no = "boot")]

# calculate average coverage
sim_summary <- sim_dt_lng[, .(coverage = mean(covered)), by = .(method, parameter, scenario)]
# Calculate Monte-Carlo standard error:
sim_summary[, se := binomial_se(phat = coverage, n = 500)]

sim_summary[, method := factor(method, levels = c("wald", "boot"))]
sim_summary_round <- sim_summary[, .(method, 
                                     parameter, 
                                     scenario, 
                                     coverage = round(coverage, 2),
                                     se = round(se, 2))]
setorderv(sim_summary_round, c("method", "scenario", "parameter"))

sim_summary_round
# save:
fwrite(x = sim_summary_round,
       file = file.path(results_path, "BFMM-paper-coverage_sim.csv"))

