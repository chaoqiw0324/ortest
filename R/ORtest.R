all_numeric <- function(df) {
  all_numeric <- all(sapply(df, is.numeric))
  return(all_numeric)
}

basic_function <- function(x, y, S=NULL,method = "linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE) {
  valid_methods <- c("sl", "linear", "linear_int")
  
  # Check if the provided method is valid
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose one of: 'sl', 'linear', or 'linear_int'.")
  }
  # Check if x and y are binary or continuous
  if(! is.numeric(x) | ! is.numeric(y)){
    stop("Error: The x and y must be numeric data")
  }
  if(! all_numeric(S)){
    stop("Error: S must be empty or contain only numeric values .")
  }
  is_binary_x <- FALSE
  
  is_binary_y <- FALSE
  # Check if the variable is a factor with exactly two levels
  if (is.factor(x) && nlevels(x) == 2) {
    # Convert factor to numeric 0 and 1
    x <- as.numeric(x) - 1
    is_binary_x <- TRUE
  }
  
  # Check if the variable is numeric with exactly two unique values
  if (is.numeric(x) && length(unique(x)) == 2) {
    # Identify the two unique values
    unique_values <- unique(x)
    # Convert to 0 and 1
    x <- ifelse(x == unique_values[1], 0, 1)
    is_binary_x <- TRUE
  }
  
  # Check if the variable is a factor with exactly two levels
  if (is.factor(y) && nlevels(y) == 2) {
    # Convert factor to numeric 0 and 1
    y <- as.numeric(y) - 1
    is_binary_y <- TRUE
  }
  
  # Check if the variable is numeric with exactly two unique values
  if (is.numeric(y) && length(unique(y)) == 2) {
    # Identify the two unique values
    unique_values <- unique(y)
    # Convert to 0 and 1
    y <- ifelse(y == unique_values[1], 0, 1)
    is_binary_y <- TRUE
  }
 
  # Check if the dimension of S is big (> 4)
  
  
  if (is_binary_x) {
    exp_bin = TRUE
  } else {
    exp_bin = FALSE
  }
  
  if (is_binary_y) {
    out_bin = TRUE
  } else {
    out_bin = FALSE
  }
  
  
  if(method=="sl"){
    res <- psi_hat_sl(y=y,x=x,S=S,subset = NULL,out_bin = out_bin,exp_bin=exp_bin,sl=sl,kfolds = kfolds,cross_fitting = cross_fitting)      
  }else if (method=="linear") {
    if (!cross_fitting) {
      res <- psi_hat_linear(y=y,x=x,S=S,subset = NULL,out_bin = out_bin,exp_bin=exp_bin)
    }else{
      stop("The cross-fitting is not available for linear method")  
    }    
  }else{
    if (!cross_fitting) {
      res <- psi_hat_linear_int(y=y,x=x,S=S,subset = NULL,out_bin = out_bin,exp_bin=exp_bin,two_way = two_way,three_way = three_way)
    }else{
      stop("The cross-fitting is not available for linear method")  
    }   
    
    
  }
  return(res)
}

multi_level <- function(x, y, S,method = "linear",sl=c(), cross_fitting=FALSE, kfolds=5,two_way=FALSE,three_way=FALSE) {
  valid_methods <- c("sl", "linear", "linear_int")
  
  # Check if the provided method is valid
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose one of: 'sl', 'linear', or 'linear_int'.")
  }
  # Check if x and y are binary or continuous
  if(! is.numeric(x) | ! is.numeric(y)){
    stop("Error: The x and y must be numeric data")
  }
  if(! all_numeric(S)){
    stop("Error: S must be empty or contain only numeric values .")
  }
  dat_y <- NULL
  dat_x <- NULL
  # If the variable is a factor with more than one level
  if (is.factor(y) && nlevels(y) > 2) {
    # Use model.matrix to create dummy variables
    dummies <- model.matrix(~ y - 1)
    colnames(dummies) <- paste0("y", "_", levels(y))
    dat_y <- as.data.frame(dummies)
  }
  
  # If the variable is numeric with less than 5 unique values
  if (is.numeric(y) && length(unique(y)) > 2 && length(unique(y)) < 5) {
    # Convert numeric to factor first and then create dummy variables
    y <- factor(y)
    dummies <- model.matrix(~ y - 1)
    colnames(dummies) <- paste0("y", "_", levels(y))
    dat_y <- as.data.frame(dummies)
  } 
  
  # If the variable is a factor with more than one level
  if (is.factor(x) && nlevels(x) > 2) {
    # Use model.matrix to create dummy variables
    dummies <- model.matrix(~ x - 1)
    colnames(dummies) <- paste0("x", "_", levels(x))
    dat_x <- as.data.frame(dummies)
  }
  
  # If the variable is numeric with less than 5 unique values
  if (is.numeric(x) && length(unique(x)) > 2 && length(unique(x)) < 5) {
    # Convert numeric to factor first and then create dummy variables
    x <- factor(x)
    dummies <- model.matrix(~ x - 1)
    colnames(dummies) <- paste0("x", "_", levels(x))
    dat_x <- as.data.frame(dummies)
  } 
  ########################################  
  if(is.null(dat_y)){
    if (is.null(dat_x)) {
      res <- basic_function(x, y, S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds,two_way=two_way,three_way=three_way)
    }else{
      num_x <- ncol(dat_x)
      result <- matrix(nrow = num_x,ncol=2)
      for (i in 1:num_x) {
        res <- basic_function(dat_x[,i], y, S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds,two_way=two_way,three_way=three_way)
        result[i,1] <- res[1]
        result[i,2] <- res[2]
      }
      id <- which.max(result[, 1])
      res <- result[id,]
    }
    
    
  }else{
    if (is.null(dat_x)) {
      num_y <- ncol(dat_y)
      result <- matrix(nrow = num_y,ncol=2)
      for (i in 1:num_y) {
        res <- basic_function(x, dat_y[,i], S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds,two_way=two_way,three_way=three_way)
        result[i,1] <- res[1]
        result[i,2] <- res[2]
      }
      id <- which.max(result[, 1])
      res <- result[id,]
    }else{
      num_x <- ncol(dat_x)
      num_y <- ncol(dat_y)
      result <- matrix(nrow = num_x*num_y,ncol=2)
      for (i in 1:num_x) {
        for (j in 1:num_y) {
          res <- basic_function(dat_x[,i], dat_y[,j], S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds,two_way=two_way,three_way=three_way)
          result[(i-1)*num_y+j,1] <- res[1]
          result[(i-1)*num_y+j,2] <- res[2]
        }
      }
      id <- which.max(result[, 1])
      res <- result[id,]
    }
    
  }
  return(res)
  
}



#' Super-Learning Based Conditional Independence Test
#' 
#' We test whether \code{x} and \code{y} are conditionally associated, given
#' \code{S}, using specific methods.
#' 
#' @param x The position of the exposure variable 
#' @param y The position of the exposure variable 
#' @param S The position of the conditional variable set
#' @param suffStat A list with elements needed for the test:
#'   \describe{
#'     \item{"dat"}{Data matrix where columns represent variables.}
#'     \item{"method"}{The method used for the test: "sl", "linear", "linear_int"}
#'     \item{"sl"}{Specifies the super learning model. If \code{y} is numeric, the default method is linear regression model, mean response,MARS, random forest and XGBoost. If \code{y} is binary, the default method is LDA, mean response,MARS, random forest and XGBoost.}
#'     \item{"cross_fitting"}{Logical; whether cross-fitting is used in super learning.}
#'     \item{"kfolds"}{Number of folds for cross-fitting.}
#'     \item{"two_way"}{Logical; whether to include two-way interactions between variables in \code{S}.}
#'     \item{"three_way"}{Logical; whether to include three-way interactions between variables in \code{S}.}
#'   }
#' @details This function test whether \code{x} and \code{y} is independent conditional on
#' \code{S}. The final result only includes p value for model \code{y ~ x + S}. 
#' 
#' @return p-value for model \code{y ~ x + S}.
#' @export
#' @name ortest_function
ortest <- function(x,y,S,suffStat) {
  
  # Extract the positions of x, y, and S from the location vector
  x_pos <- x
  y_pos <- y
  s_pos <- S
  
  # Ensure that x and y positions are different
  if (x_pos == y_pos) {
    stop("x and y should be at different positions.")
  }
  dat <- suffStat$dat
  method <- suffStat$method
  sl <- suffStat$sl
  cross_fitting <- suffStat$cross_fitting
  kfolds <- suffStat$kfolds
  two_way <- suffStat$two_way
  three_way <- suffStat$three_way
  
  ### p value 1
  # Extract x, y, and S variables from the dataframe
  x <- dat[, x_pos]
  y <- dat[, y_pos]
  S <- dat[, s_pos]
  res1 <- multi_level(x = x,y= y,S=S,method = method,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds,two_way=two_way,three_way=three_way)
  
  
  
  # Two-sided p-value
  estimate <- res1[1]
  sd_estimate <- res1[2]
  z <- (estimate - 0) / sd_estimate
  
  # Calculate two-sided p-value
  p_value_two_sided <- 2 * pnorm(-abs(z))
  return(p_value_two_sided)
}