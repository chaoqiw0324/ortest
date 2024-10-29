generate_dataset <- function(n = 1000) {
  # Set seed for reproducibility
  set.seed(123)
  
  # Generate numeric columns
  col1 <- rnorm(n)
  col4 <- rnorm(n)
  col7 <- rnorm(n)
  col9 <- rnorm(n)
  
  # Generate binary columns (0 or 1)
  col2 <- sample(c(0, 1), n, replace = TRUE)
  col3 <- sample(c(0, 1), n, replace = TRUE)
  col8 <- sample(c(0, 1), n, replace = TRUE)
  
  # Generate multi-level categorical columns (1, 2, 3, 4)
  col5 <- sample(1:4, n, replace = TRUE)
  col6 <- sample(1:4, n, replace = TRUE)
  col10 <- sample(1:4, n, replace = TRUE)
  
  # Combine all columns into a data frame
  data <- data.frame(
    col1 = col1,
    col2 = col2,
    col3 = col3,
    col4 = col4,
    col5 = col5,
    col6 = col6,
    col7 = col7,
    col8 = col8,
    col9 = col9,
    col10 = col10
  )
  
  return(data)
}

# Generate the dataset
set.seed(2024)
n=1000
data <- generate_dataset(n)

test_that("ortest return the correct number", {
  # Basic addition
  expect_equal(round(ortest(1,5,c(2:4,6:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE)),3), 0.334)
  expect_equal(round(ortest(5,2,c(1,3:4,6:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE)),3), 0.790)
  expect_equal(round(ortest(5,6,c(1:4,7:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE)),3), 0.555)
})

# ## c to m
# ortest(1,5,c(2:4,6:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE))
# 
# ## m to b 
# ortest(5,2,c(1,3:4,6:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE))
# 
# ## m to m
# ortest(5,6,c(1:4,7:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE))
