library(testthat)

test_that(".get_variable extracts columns accurately and handles all inputs", {
  df <- data.frame(
    id = 1:5,
    val = c(10, 20, 30, 40, 50),
    stringsAsFactors = FALSE
  )

  get_var <- if (exists(".get_variable", mode = "function")) .get_variable else fastsae:::.get_variable
  # 1. By character column name
  expect_equal(get_var(df, "val"), df$val)
  expect_error(get_var(df, "not_existing"), 'variable "not_existing" is not found in the data')

  # 2. By formula
  expect_equal(get_var(df, ~val), df$val)
  expect_error(get_var(df, ~not_existing), "formula does not reference a valid single column in data")
  expect_error(get_var(df, ~id + val), "formula does not reference a valid single column in data")

  # 3. Vector input directly
  vec_matching <- c(1, 2, 3, 4, 5)
  expect_equal(get_var(df, vec_matching), vec_matching)

  # 4. Invalid length or type
  expect_error(get_var(df, c(1, 2)), "variable is not valid or length does not match data")
  expect_error(get_var(df, list(a = 1)), "variable is not valid or length does not match data")
})
