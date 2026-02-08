test_that("ortest returns a valid p-value", {
  set.seed(1)
  n_val <- 200
  
  # Simple mixed-type data: y numeric, x binary, S numeric
  S1 <- runif(n_val)
  S2 <- runif(n_val)
  x  <- rbinom(n_val, 1, plogis(S1 - 0.5))
  y  <- 1 + 2 * x + S2 + rnorm(n_val)
  dat_matrix <- cbind(y, x, S1, S2)
  
  # 2. Run the function
  result <- ortest(
    x = 2, y = 1, S = c(3, 4),
    list(
      dat = dat_matrix, 
      method = "linear", 
      sl = NULL, 
      cross_fitting = FALSE, 
      error_catch = FALSE,
      kfolds = 5, 
      two_way = FALSE, 
      three_way = FALSE
    )
  )
  
  # 3. The Actual Assertions (The "Tests")
  expect_type(result, "double")
  expect_true(is.finite(result))
  expect_true(result >= 0 && result <= 1,
              info = "ortest should return a p-value in [0, 1]")
})
