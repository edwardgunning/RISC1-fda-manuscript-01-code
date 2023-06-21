# Packages ----------------------------------------------------------------
library(data.table) # CRAN v1.14.0
library(ggplot2)    # CRAN v3.3.5
library(fda)        # CRAN v5.5.1
library(multifamm)  # CRAN v0.1.1
library(sparseFLMM) # CRAN v0.4.1


source(file = here::here(
  "code",
  "functions",
  "functions-helper-smoothing.R"
))
source(here::here(
  "code/functions/rough_fit_mfamm_model.R"
))

data_path <- here::here("data")
subject_side_coef_reg <- readRDS(file.path(data_path, "subject_side_coef_reg.rds"))

# set up basis to represent the data.
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)

subject_side_eval_reg <- copy(subject_side_coef_reg)
subject_side_eval_reg[, paste0(0:100) := {
  coef_mat <- coef_to_mat(.SD)
  fd_object <- fd(coef = coef_mat, basisobj = bspl80)
  fd_eval <- eval.fd(evalarg = 0:100, fdobj = fd_object)
  return_list <- get_list_of_rows(fd_eval)
  names(return_list) <- 0:100
  return_list},
  .SDcols = paste0("lm_coef_", 1:80)]
subject_side_eval_reg[, paste0("lm_coef_", 1:80) := NULL]

rm(subject_side_coef_reg)


# Prepare the data for fitting the multiFAMM ------------------------------
# reshape to long form to fit the multiFAMM
subject_side_eval_reg_lng <- melt.data.table(subject_side_eval_reg,
                                             measure.vars = paste0(0:100),
                                             variable.name = "t",
                                             variable.factor = FALSE,
                                             value.name = "y_vec")

subject_side_eval_reg_lng[, t := as.numeric(t)]

names(subject_side_eval_reg_lng)
mfamm_data <- subject_side_eval_reg_lng[location %in% c("Hip", "Knee"),
                                        c("subject_id", 
                                          "side",
                                          "location",
                                          "t", 
                                          "y_vec",
                                          "retrospective_injury_status",
                                          "sex",
                                          "speed_cent",
                                          "age_cent",
                                          "weight_kg_cent",
                                          "height_cm_cent")]

mfamm_data[, sex := fifelse(sex == "female", 1, 0)]
mfamm_data[, ris2 := fifelse(retrospective_injury_status == "injured_greater_than_2_yr", 1, 0)]
mfamm_data[, ris3 := fifelse(retrospective_injury_status == "injured_1_to_2_yr", 1, 0)]
mfamm_data[, ris4 := fifelse(retrospective_injury_status == "injured_less_than_1_yr", 1, 0)]
mfamm_data[, retrospective_injury_status := NULL]
mfamm_data[, subject_id := as.integer(as.factor(subject_id))]
mfamm_data[, n_long := as.integer(factor(interaction(side, subject_id)))]

setnames(mfamm_data,
         old = c("subject_id", "location", "sex", "speed_cent", 
                 "ris2", "ris3", "ris4", "age_cent",
                 "weight_kg_cent", "height_cm_cent"),
         new = c("subject_long", "dim",
                 paste0("covariate.", 1:8)))

mfamm_data[, dim := factor(dim)]
stopifnot(is.numeric(mfamm_data$y_vec))
stopifnot(is.numeric(mfamm_data$t))
stopifnot(is.factor(mfamm_data$dim))
stopifnot(is.integer(mfamm_data$subject_long))
stopifnot(is.integer(mfamm_data$n_long))
stopifnot(purrr::map_lgl(mfamm_data[, paste0("covariate.", 1:8)],
                         is.numeric))
setorderv(mfamm_data, c("dim", "subject_long", "n_long", "t"))
B <- 10
results_list <- vector(mode = "list", length = B)
time_vec <- vector("numeric", length = B)
seeds_list <- vector(mode = "list", length = B)
set.seed(1)

for(mind in seq_len(B)) {
  print(paste0("iteration ", mind, " of ", B))
  # Do bootstrap of subjects:
  seeds_list[[mind]] <- .Random.seed
  mfamm_data_b <- copy(mfamm_data)
  subjects_sample <- sample(unique(mfamm_data_b$subject_long), size = length(unique(mfamm_data_b$subject_long)), replace = TRUE)
  mfamm_data_b_new <- data.table()
  for(sub_ind in seq_along(subjects_sample)) {
    mfamm_data_b_new <- rbind(mfamm_data_b_new, cbind(mfamm_data_b[subject_long == subjects_sample[sub_ind]], subject_long_new = sub_ind))
  }
  mfamm_data_b_new[, subject_long := subject_long_new]
  mfamm_data_b_new[, n_long := as.integer(as.factor(interaction(subject_long, side)))]
  # Fit model:
  time_vec[mind] <- system.time(results_list[[mind]] <- fit_multifamm(data = mfamm_data_b_new,
                                                                      bf_covs = c(5, 5),
                                                                      mfpc_cutoff = 0.95))["elapsed"]
  gc()
}

saveRDS(list(time_vec=time_vec, results_list=results_list),
        file = here::here("outputs/results/mfamm-comparison-bootstrap-02.rds"))

