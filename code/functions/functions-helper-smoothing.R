# Count non-NA values -----------------------------------------------------
# Function to count the number of non-NA values in an object
# Will be used to count number of non-NA values for our purpose.

## Function ---------------------------------------------------------------
count_non_na_values <- function(x) sum(!is.na(x))

## Test -------------------------------------------------------------------
# Test the function on some examples:
testthat::expect_equal(count_non_na_values(rep(c(1, NA), each = 5)),
                       5)
testthat::expect_equal(count_non_na_values(rep(c(1, NA), times = 5)),
                       5)
testthat::expect_equal(count_non_na_values(NA), 0)
testthat::expect_equal(count_non_na_values(1), 1)
testthat::expect_equal(count_non_na_values(
  data.frame(x = 1:5, y = rep(NA, 5))
  ), 5)


# Change matrix to list where elements are its columns --------------------

## Function ---------------------------------------------------------------
get_list_of_columns <- function(x) {
  stopifnot(is.matrix(x))
  lapply(seq_len(ncol(x)), function(i) x[, i])
  }


## Test -------------------------------------------------------------------
# test 1
test_vector <- 1:16
test1_matrix_col <- matrix(test_vector, ncol = 4, nrow = 4, byrow = FALSE)
test1_result <- list(1:4, 5:8, 9:12, 13:16)
testthat::expect_equal(get_list_of_columns(test1_matrix_col),
                       test1_result)

# test 2
# Should fail if not a data frame our matrix that is supplied:
testthat::expect_error(get_list_of_columns(test_vector))

# test 3
# should give error if data frame is supplied
testthat::expect_error((get_list_of_columns(data.frame(a = test_vector))))



# Change matrix to list where elements are its rows -----------------------


## Function ---------------------------------------------------------------
get_list_of_rows <-  function(x) {
  stopifnot(is.matrix(x))
  lapply(seq_len(nrow(x)), function(i) x[i,])
  }


## Test --------------------------------------------------------------------

# test 1
# re-use some of last test 1
# fill matrix by rows this time, but re-use result (should be same.)
test1_matrix_row <- matrix(test_vector, ncol = 4, nrow = 4, byrow = TRUE)
testthat::expect_equal(get_list_of_rows(test1_matrix_row),
                       test1_result)
# test 2
# Should give error if vector is supplied:
testthat::expect_error(get_list_of_rows(test_vector))

# test 3
# should give error if data frame is supplied
testthat::expect_error((get_list_of_rows(data.frame(a = test_vector))))


# Function to convert coefficients stored in dt to matrix -----------------

## Function ---------------------------------------------------------------

coef_to_mat <- function(dt) { 
  
  stopifnot(data.table::is.data.table(dt))
  
  # This is probably the most important part of the function
  
  # as we will be using this function to transpose coefficients
  # stored in a data table they must all be numeric
  # this will 
  if(!all(purrr::map_lgl(.x = dt, .f = is.numeric))) {
    stop("Not all data.table columns are numeric, cannot convert to matrix.")
  }
  
  # use data table functions to convert it to a matrix
  # and use arguments to turn columns of 
  return_matrix <- data.table:::as.matrix.data.table(
    data.table::transpose(l = dt, keep.names = "rm"),
    rownames = "rm")
  
  colnames(return_matrix) <- paste0("fun_rep", seq_len(ncol(return_matrix)))
  
  stopifnot(is.numeric(return_matrix))
  
  return_matrix
  # actually, just t(dt) would work but this is being VERY careful
  
}



## Test -------------------------------------------------------------------

# Test by 1) manually creating an result to compare to
# and compare to the conventional solution with 
test_dt <- data.table::data.table(a = 1:5, b = 5:1)
test_result <- matrix(c(1:5, 5:1), byrow = TRUE, ncol = 5, nrow = 2)
test_result_2 <- t(test_dt)

rownames(test_result) <- c("a", "b")
colnames(test_result) <- colnames(test_result_2) <- paste0("fun_rep", 1:5)

testthat::expect_equal(object = coef_to_mat(dt = test_dt), 
                       expected = test_result)
testthat::expect_equal(object = coef_to_mat(dt = test_dt), 
                       expected = test_result_2)

# Now, add a character column, it should a throw error
test_dt$char_col <- "c"
testthat::expect_error(object = coef_to_mat(dt = test_dt))



# Remove objects created for tests ----------------------------------------
# remove objects from R environment:
rm(list = c("test_vector", "test1_matrix_col", "test1_matrix_row",
            "test1_result", "test_dt", "test_result", "test_result_2"))
