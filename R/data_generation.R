trans <- function(z){
  result = exp(-(z^(2))/2)*sin(2*z)
  return(result)
}

expit <- function(x) {
  1 / (1 + exp(-x))
}

data_generation_cont_bin <- function(n,hypothesis="null",S="a"){
  allowed_methods <- c(
    "a","b","c","d","e"
  )
  
  # Check if all elements in 'sl' exist in allowed methods
  if (!all(S %in% allowed_methods)) {
    stop("The models for S doesn't exist")
  }
  
  S1 <- runif(n,0,1)
  S2 <- runif(n,0,1)
  S3 <- runif(n,0,1)
  S4 <- runif(n,0,1)
  S5 <- runif(n,0,1)
  if (S == "a") {
    S <- S1 + S2 + S3 + S4 + S5 - 2.5 # Adjusted
  } else if (S == "b") {
    S <- (S1 * S2 + S3 * S4 + S5) / 2 - 0.5 # Adjusted
  } else if (S == "c") {
    S <- (S1 * S2 * S3 + S4 + S5) # Adjusted
  } else if (S == "d") {
    S <- (S1 * S2 + S3 * S4 + S5 + S1^2 + S5^2 + S3^3) / 5 -0.4 # Adjusted
  } else {
    S <- S1 + S2 + trans(S3) + S4 + trans(S5) - 2 # Adjusted
  }   
  prob <- expit(S)
  # x <- rbinom(n,1,prob)
  x <- rbinom(n,1,prob)
  if (hypothesis=="null")  {
    y <- 1 + S + rnorm(n,0,5) 
  }else if(hypothesis=="alternative"){
    # y <- 1 + 5 * x + S + x * S + rnorm(n,0,1)
    # y <- 1 + 5 * x + S  + rnorm(n,0,1)
    y <- 1 + 5 * x + S  + rnorm(n,0,5)
  }
  data <- tibble::as_tibble(cbind(y, x, S1, S2, S3, S4, S5, prob))
  return(data)
}


create_multi_level_factor <- function(data, levels = 3, labels = NULL) {
  if (!is.numeric(data)) {
    stop("Input data must be numeric.")
  }
  
  if (levels < 2) {
    stop("Number of levels must be at least 2.")
  }
  
  # Calculate quantile cut points
  probs <- seq(0, 1, length.out = levels + 1)
  cut_points <- quantile(data, probs, na.rm = TRUE, type = 7)
  
  # Create labels if not provided
  if (is.null(labels)) {
    labels <- seq_len(levels)
  }
  
  # Cut data into levels
  factor_data <- cut(
    data,
    breaks = cut_points,
    labels = labels,
    include.lowest = TRUE
  )
  
  return(factor_data)
}


#' Generate Test Data
#'
#' @param N the number of the records
#'
#' @return test data that includes binary, numeric and multi-level factor variable
#' @export

generate_test_data <- function(N=2000){
  if (!requireNamespace("pcalg", quietly = TRUE)) {
    stop("Package 'pcalg' is required for generate_test_data(). Please install it.",
         call. = FALSE)
  }
  rDAG <- pcalg::randomDAG(n = 10, prob= 0.7, lB = 0.4, uB = 1)
  
  data <- pcalg::rmvDAG(N,rDAG,errDist="normal")
  data[,2] <- ifelse(data[,2]<=median(data[,2]),0,1)
  data[,3] <- ifelse(data[,3]<=median(data[,3]),0,1)
  data[,8] <- ifelse(data[,8]<=median(data[,8]),0,1)
  data[,5] <- create_multi_level_factor(data[,5], levels = 4)
  data[,6] <- create_multi_level_factor(data[,6], levels = 3)
  data[,10] <- create_multi_level_factor(data[,10], levels = 3)
  
  # for(i in 1:N){
  #   if(data[i,6] <= quantile(data[,6],.33)) data[i,6] <- 1
  #   else if(quantile(data[,6],.33) < data[i,6] & data[i,6] <= quantile(data[,6],.66)) data[i,6] <- 2
  #   else data[i,6] <- 3
  # }
  # 
  # for(i in 1:N){
  #   if(data[i,10] <= quantile(data[,10],.33)) data[i,10] <- 1
  #   else if(quantile(data[,10],.33) <= data[i,10] & data[i,10] <= quantile(data[,10],.66)) data[i,10] <- 2
  #   else data[i,10] <- 3
  # }
  # 
  # for(i in 1:N){
  #   if(data[i,5] <= quantile(data[,5],.25)) data[i,5] <- 1
  #   else if(quantile(data[,5],.25) <= data[i,5] & data[i,5] <= quantile(data[,5],.50)) data[i,5] <- 2
  #   else if(quantile(data[,5],.50) <= data[i,5] & data[i,5] <= quantile(data[,5],.75)) data[i,5] <- 3
  #   else data[i,5] <- 4
  # }
  
  return(data)
}

