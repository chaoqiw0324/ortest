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
  rDAG <- randomDAG(n = 10, prob= 0.7, lB = 0.4, uB = 1)
  
  data <- rmvDAG(N,rDAG,errDist="normal")
  data[,2] <- ifelse(data[,2]<=median(data[,2]),0,1)
  data[,3] <- ifelse(data[,3]<=median(data[,3]),0,1)
  data[,8] <- ifelse(data[,8]<=median(data[,8]),0,1)
  data[,5] <- create_multi_level_factor(data[,5], levels = 4)
  data[,6] <- create_multi_level_factor(data[,6], levels = 3)
  data[,8] <- create_multi_level_factor(data[,8], levels = 3)
  
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

