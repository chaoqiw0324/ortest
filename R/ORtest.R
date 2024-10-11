all_numeric <- function(df) {
  all_numeric <- all(sapply(df, is.numeric))
  return(all_numeric)
}

basic_function <- function(x, y, S=NULL,sl=NULL,cross_fitting=FALSE,kfolds=5) {
  # Check if x and y are binary or continuous
  if(! is.numeric(x) | ! is.numeric(y)){
    cat("The x and y must be numeric data")
  }
  if(! all_numeric(S)){
    cat("The S must be numeric data or empty")
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
  S2 <- as_tibble(S)
  colnames(S2) <- paste0("X", colnames(S2))
  is_big_S <- FALSE
  if(is.null(S2)){
    is_big_S <- FALSE
  }else{
    is_big_S <- ifelse( ncol(S2) > 4,TRUE, FALSE)
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
  
  roots="uni"
  
  if(is_big_S){
    res <- psi.hat_sl(y=y,A=x,S=S2,subset = NULL,out.bin = out_bin,exp.bin=exp_bin,exp.scalar = FALSE,root=roots,sl=sl,kfolds = kfolds,cross_fitting = cross_fitting)      
  }else{
    if (!cross_fitting) {
      res <- psi.hat_linear(y=y,A=x,S=S2,subset = NULL,out.bin = out_bin,exp.bin=exp_bin,exp.scalar = FALSE,root=roots)
    }else{
      cat("The cross-fitting is not available for columns of S is smaller than 4")  
    }
    
    
  }
  # Two-sided p-value
  # estimate <- res[1]
  # sd_estimate <- res[2]
  # S <- (estimate - 0) / sd_estimate
  # 
  # # Calculate two-sided p-value
  # p_value_two_sided <- 2 * pnorm(-abs(S))
  # 
  # # Output
  # return(list(p_value = p_value_two_sided, alpha = alpha, independent = p_value_two_sided > alpha))
  return(res)
}

multi_level <- function(x, y, S,sl=c(), cross_fitting=FALSE, kfolds=5) {
  # Check if x and y are binary or continuous
  if(! is.numeric(x) | ! is.numeric(y)){
    cat("The x and y must be numeric data")
  }
  if(! all_numeric(S)){
    cat("The S must be numeric data or empty")
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
      res <- basic_function(x, y, S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
    }else{
      num_x <- ncol(dat_x)
      result <- matrix(nrow = num_x,ncol=2)
      for (i in 1:num_x) {
        res <- basic_function(dat_x[,i], y, S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
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
        res <- basic_function(x, dat_y[,i], S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
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
          res <- basic_function(dat_x[,i], dat_y[,j], S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
          result[(i-1)*num_y+j,1] <- res[1]
          result[(i-1)*num_y+j,2] <- res[2]
        }
      }
      id <- which.max(result[, 1])
      res <- result[id,]
    }
    
  }
  # estimate <- res[1]
  # sd_estimate <- res[2]
  # S <- (estimate - 0) / sd_estimate
  # 
  # # Calculate two-sided p-value
  # p_value_two_sided <- 2 * pnorm(-abs(S))
  # 
  # # Output
  # return(list(p_value = p_value_two_sided, alpha = alpha, independent = p_value_two_sided > alpha))
  return(res)
  
}

#' Conditional Independence Test
#' We test whether \code{x} and \code{y} are conditionally associated, given
#' \code{S} 
#' 
#' @param x The position of the exposure variable 
#' @param y The position of the exposure variable 
#' @param S The position of the conditional variable set
#' @param suffStat A list with four elements, "dat" is the dat matrix ,"sl" specifies the super learning model,"cross_fitting" indicates whether super learning is used and "kfolds" decides the folds for cross-fitting.
#'
#' @details This function test whether \code{x} and \code{y} is independent conditional on
#' \code{S}. The final results includes p-value for model \code{y ~ x + S} 
#' and \code{x ~ y+S}.
#' @return p-value for model \code{y ~ x + S} and \code{x ~ y+S}.
#' @export
#'
ORtest <- function(x,y,S,suffStat) {
  
  # Extract the positions of X, Y, and Z from the location vector
  x_pos <- x
  y_pos <- y
  z_pos <- S
  
  # Ensure that X and Y positions are different
  if (x_pos == y_pos) {
    stop("X and Y should be at different positions.")
  }
  dat <- suffStat$dat
  sl <- suffStat$sl
  cross_fitting <- suffStat$cross_fitting
  kfolds <- suffStat$kfolds
  
  
  ### p value 1
  # Extract X, Y, and Z variables from the dataframe
  X <- dat[, x_pos]
  Y <- dat[, y_pos]
  Z <- dat[, z_pos]
  res1 <- multi_level(X = X,Y= Y,Z=Z,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
  ### p value 2
  # Extract X, Y, and Z variables from the dataframe
  Y <- dat[, x_pos]
  X <- dat[, y_pos]
  Z <- dat[, z_pos]
  res2 <- multi_level(X = X,Y= Y,Z=Z,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
  # Return the X, Y, and Z variables in a list
  return(list(res1,res2))
}

#' Conditional Independence Test
#' We test whether \code{x} and \code{y} are conditionally associated, given
#' \code{S} 
#' 
#' @param x The position of the exposure variable 
#' @param y The position of the exposure variable 
#' @param S The position of the conditional variable set
#' @param suffStat A list with four elements, "dat" is the dat matrix ,"sl" specifies the super learning model,"cross_fitting" indicates whether super learning is used and "kfolds" decides the folds for cross-fitting.
#'
#' @details This function test whether \code{x} and \code{y} is independent conditional on
#' \code{S}. The final result only includes p value for model \code{y ~ x + S}. 
#' 
#' @return p-value for model \code{y ~ x + S}.
#' @export
#'
ORtest_single <- function(x,y,S,suffStat) {
  
  # Extract the positions of x, y, and S from the location vector
  x_pos <- x
  y_pos <- y
  s_pos <- S
  
  # Ensure that x and y positions are different
  if (x_pos == y_pos) {
    stop("x and y should be at different positions.")
  }
  dat <- suffStat$dat
  sl <- suffStat$sl
  cross_fitting <- suffStat$cross_fitting
  kfolds <- suffStat$kfolds
  
  
  ### p value 1
  # Extract x, y, and S variables from the dataframe
  x <- dat[, x_pos]
  y <- dat[, y_pos]
  S <- dat[, s_pos]
  res1 <- multi_level(x = x,y= y,S=S,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
  
  
  
  # Two-sided p-value
  estimate <- res1[1]
  sd_estimate <- res1[2]
  z <- (estimate - 0) / sd_estimate
  
  # Calculate two-sided p-value
  p_value_two_sided <- 2 * pnorm(-abs(z))
  return(p_value_two_sided)
}