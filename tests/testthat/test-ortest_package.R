test_that("multiplication works", {
  n_val <- 1000
  hypo <- "alternative"
  s_val <- "e"
  
  test_data <- data_generation_cont_bin(n = n_val, hypothesis = hypo, S = s_val)
  dat_matrix <- as.matrix(test_data)
  
  # 2. Run the function
  result <- ortest(
    2, 1, c(3, 4, 5, 6, 7),
    list(
      dat = dat_matrix, 
      method = "linear", 
      sl = NULL, 
      cross_fitting = FALSE, 
      kfolds = 5, 
      two_way = FALSE, 
      three_way = FALSE
    )
  )
  
  # 3. The Actual Assertions (The "Tests")
  expect_type(result, "double")
  expect_true(all(result >= 0 & result <= 1), 
              info = "The output of ortest should be a probability between 0 and 1")
})
