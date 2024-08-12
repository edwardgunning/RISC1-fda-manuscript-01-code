# ------------------------------------------------------------------------#
# Add table with participant characteristics for review.
# ------------------------------------------------------------------------#

# Packages and functions needed -------------------------------------------
library(data.table)   # CRAN v1.14.0
library(fda)          # CRAN v5.5.1
library(modelsummary) # CRAN v1.4.1

# Helper functions:
# Custom functions to turn columns/ rows of matrix into elements of a list
source(file = here::here(
  "code",
  "functions",
  "functions-helper-smoothing.R"
))

outputs_path <- here::here("outputs")
save_path <- here::here("data", "subject-side")

# Import Data -------------------------------------------------------------
# Path to data:
# contains coefficients for registered and unregistered functions
data_path <- here::here(
  "data",
  "subject_side_coef_reg.rds"
)

# Read it in the stored basis-coefficient data:
risc1_basis_coef <- readRDS(file = data_path)

# Define basis to combine with stored coefficients
bspl80 <- fda::create.bspline.basis(
  rangeval = c(0, 100),
  nbasis = 80,
  norder = 4)


# Create subset of the data set to work with ------------------------------

# Start by looking at flexion angles:
risc1_basis_coef_subset <- risc1_basis_coef[
  subject_id != "P_4092" & # This subject only has strides for one side
    plane_of_motion == "fle" & # Let's focus on flexion for now
    !is.na(retrospective_injury_status)]


subset_info_dt <- risc1_basis_coef_subset
table1_dt <- subset_info_dt[, .(
  ris = retrospective_injury_status[1],
  speed = self_selected_speed_kmph[1],
  sex = sex[1],
  age = age[1],
  weight = weight_kg[1],
  height = height_cm[1]),
  by = subject_id]
table1_dt[, ris := factor(ris, levels = c("never_injured",
                                          "injured_greater_than_2_yr",
                                          "injured_1_to_2_yr", 
                                          "injured_less_than_1_yr"),
                          labels = c("Never Injured",
                                     "Injured $>2$ yr. ago",
                                     "Injured $1-2$ yr. ago",
                                     "Injured $<1$ yr. ago"))]
table1_dt[, sex := factor(sex, levels = c("male", "female"), labels = c("Male", "Female"))]
setnames(table1_dt,
         old = c("ris", "speed", "sex", "age", "weight", "height"),
         new = c("Retrospective Injury Status", "Speed (\\si{\\km \\per \\hour})", "Sex", "Age (years)", "Weight (kg)", "Height (cm)"))
table1 <- datasummary_balance(~ 1, data = table1_dt[, - c("subject_id")],
                              output = "latex", escape = FALSE)
table1 <- stringr::str_replace(table1,
                               pattern = "& N & Pct.\\\\\\\\\n",
                               replacement = "& N & Pct.\\\\\\\\\n\\\\midrule\n")
table1 <- stringr::str_replace(table1, pattern = "N & Pct.", 
                               replacement = "\\\\textbf{N} & $\\\\mathbf{\\\\mathbf{(\\\\%)}}$")
table1 <- stringr::str_replace(table1, pattern = "Mean & Std. Dev.", 
                               replacement = "\\\\textbf{Mean} & \\\\textbf{Std. Dev.}")
table1 <- stringr::str_replace(table1, pattern = "\\\\end\\{tabular\\}\n\\\\end\\{table\\}",
                               replacement = "\\\\end\\{tabular\\}\\\n\\\\caption\\{Summary characteristics of the participants in the RISC dataset included in this analysis.\\}\n\\\\label\\{tab:tab1chpt3.\\}\n\\\\end\\{table\\}")
writeLines(text = table1, con = file.path(outputs_path, "tables", "table_1_characteristics.tex"))

