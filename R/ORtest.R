basic_function <- function(X, Y, Z=NULL,sl=NULL,cross_fitting=FALSE,kfolds=5) {
  # Check if X and Y are binary or continuous
  if(! is.numeric(X) | ! is.numeric(Y)){
    cat("The X and Y must be numeric data")
  }
  if(! all_numeric(Z)){
    cat("The Z must be numeric data or empty")
  }
  is_binary_X <- FALSE
  
  is_binary_Y <- FALSE
  # Check if the variable is a factor with exactly two levels
  if (is.factor(X) && nlevels(X) == 2) {
    # Convert factor to numeric 0 and 1
    X <- as.numeric(X) - 1
    is_binary_X <- TRUE
  }
  
  # Check if the variable is numeric with exactly two unique values
  if (is.numeric(X) && length(unique(X)) == 2) {
    # Identify the two unique values
    unique_values <- unique(X)
    # Convert to 0 and 1
    X <- ifelse(X == unique_values[1], 0, 1)
    is_binary_X <- TRUE
  }
  
  # Check if the variable is a factor with exactly two levels
  if (is.factor(Y) && nlevels(Y) == 2) {
    # Convert factor to numeric 0 and 1
    Y <- as.numeric(Y) - 1
    is_binary_Y <- TRUE
  }
  
  # Check if the variable is numeric with exactly two unique values
  if (is.numeric(Y) && length(unique(Y)) == 2) {
    # Identify the two unique values
    unique_values <- unique(Y)
    # Convert to 0 and 1
    Y <- ifelse(Y == unique_values[1], 0, 1)
    is_binary_Y <- TRUE
  }
  Z2 <- as_tibble(Z)
  colnames(Z2) <- paste0("X", colnames(Z2))
  is_big_Z <- FALSE
  if(is.null(Z2)){
    is_big_Z <- FALSE
  }else{
    is_big_Z <- ifelse( ncol(Z2) > 4,TRUE, FALSE)
  }
  # Check if the dimension of Z is big (> 4)
  
  
  if (is_binary_X) {
    exp_bin = TRUE
  } else {
    exp_bin = FALSE
  }
  
  if (is_binary_Y) {
    out_bin = TRUE
  } else {
    out_bin = FALSE
  }
  
  roots="uni"
  
  if(is_big_Z){
    res <- psi.hat_sl_2in1(Y=Y,A=X,L=Z2,subset = NULL,out.bin = out_bin,exp.bin=exp_bin,exp.scalar = FALSE,root=roots,sl=sl,kfolds = kfolds,cross_fitting = cross_fitting)      
  }else{
    if (!cross_fitting) {
      res <- psi.hat_linear(Y=Y,A=X,L=Z2,subset = NULL,out.bin = out_bin,exp.bin=exp_bin,exp.scalar = FALSE,root=roots)
    }else{
      cat("The cross-fitting is not available for columns of Z is smaller than 4")  
    }
    
    
  }
  # Two-sided p-value
  # estimate <- res[1]
  # sd_estimate <- res[2]
  # Z <- (estimate - 0) / sd_estimate
  # 
  # # Calculate two-sided p-value
  # p_value_two_sided <- 2 * pnorm(-abs(Z))
  # 
  # # Output
  # return(list(p_value = p_value_two_sided, alpha = alpha, independent = p_value_two_sided > alpha))
  return(res)
}

multi_level <- function(X, Y, Z,sl=c(), cross_fitting=FALSE, kfolds=5) {
  # Check if X and Y are binary or continuous
  if(! is.numeric(X) | ! is.numeric(Y)){
    cat("The X and Y must be numeric data")
  }
  if(! all_numeric(Z)){
    cat("The Z must be numeric data or empty")
  }
  dat_y <- NULL
  dat_x <- NULL
  # If the variable is a factor with more than one level
  if (is.factor(Y) && nlevels(Y) > 2) {
    # Use model.matrix to create dummy variables
    dummies <- model.matrix(~ Y - 1)
    colnames(dummies) <- paste0("y", "_", levels(Y))
    dat_y <- as.data.frame(dummies)
  }
  
  # If the variable is numeric with less than 5 unique values
  if (is.numeric(Y) && length(unique(Y)) > 2 && length(unique(Y)) < 5) {
    # Convert numeric to factor first and then create dummy variables
    Y <- factor(Y)
    dummies <- model.matrix(~ Y - 1)
    colnames(dummies) <- paste0("y", "_", levels(Y))
    dat_y <- as.data.frame(dummies)
  } 
  
  # If the variable is a factor with more than one level
  if (is.factor(X) && nlevels(X) > 2) {
    # Use model.matrix to create dummy variables
    dummies <- model.matrix(~ X - 1)
    colnames(dummies) <- paste0("x", "_", levels(X))
    dat_x <- as.data.frame(dummies)
  }
  
  # If the variable is numeric with less than 5 unique values
  if (is.numeric(X) && length(unique(X)) > 2 && length(unique(X)) < 5) {
    # Convert numeric to factor first and then create dummy variables
    X <- factor(X)
    dummies <- model.matrix(~ X - 1)
    colnames(dummies) <- paste0("x", "_", levels(X))
    dat_x <- as.data.frame(dummies)
  } 
  ########################################  
  if(is.null(dat_y)){
    if (is.null(dat_x)) {
      res <- basic_function(X, Y, Z,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
    }else{
      num_x <- ncol(dat_x)
      result <- matrix(nrow = num_x,ncol=2)
      for (i in 1:num_x) {
        res <- basic_function(dat_x[,i], Y, Z,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
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
        res <- basic_function(X, dat_y[,i], Z,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
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
          res <- basic_function(dat_x[,i], dat_y[,j], Z,sl=sl,cross_fitting = cross_fitting,kfolds = kfolds)
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
  # Z <- (estimate - 0) / sd_estimate
  # 
  # # Calculate two-sided p-value
  # p_value_two_sided <- 2 * pnorm(-abs(Z))
  # 
  # # Output
  # return(list(p_value = p_value_two_sided, alpha = alpha, independent = p_value_two_sided > alpha))
  return(res)
  
}

ORtest_single <- function(x,y,S,suffStat) {
  
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
  
  
  
  # Two-sided p-value
  estimate <- res1[1]
  sd_estimate <- res1[2]
  Z <- (estimate - 0) / sd_estimate
  
  # Calculate two-sided p-value
  p_value_two_sided <- 2 * pnorm(-abs(Z))
  return(p_value_two_sided)
}