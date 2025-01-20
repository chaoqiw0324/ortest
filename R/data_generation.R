#' Generate Test Data
#'
#' @param N the number of the records
#'
#' @return test data that includes binary, numeric and multi-level factor variable
#' @export

generate_test_data <- function(N=2000){
  rDAG <- randomDAG(n = 10, prob= 0.7, lB = 0.4, uB = 1)
  
  data <- rmvDAG(N,rDAG,errDist="normal")
  data[,2] <- ifelse(data[,2]<=mean(data[,2]),0,1)
  data[,3] <- ifelse(data[,3]<=mean(data[,3]),0,1)
  data[,8] <- ifelse(data[,8]<=mean(data[,8]),0,1)
  
  for(i in 1:N){
    if(data[i,6] <= quantile(data[,6],.33)) data[i,6] <- 1
    else if(quantile(data[,6],.33) <= data[i,6] & data[i,6] <= quantile(data[,6],.66)) data[i,6] <- 2
    else data[i,6] <- 3
  }
  
  for(i in 1:N){
    if(data[i,10] <= quantile(data[,10],.33)) data[i,10] <- 1
    else if(quantile(data[,10],.33) <= data[i,10] & data[i,10] <= quantile(data[,10],.66)) data[i,10] <- 2
    else data[i,10] <- 3
  }
  
  for(i in 1:N){
    if(data[i,5] <= quantile(data[,5],.25)) data[i,5] <- 1
    else if(quantile(data[,5],.25) <= data[i,5] & data[i,5] <= quantile(data[,5],.50)) data[i,5] <- 2
    else if(quantile(data[,5],.50) <= data[i,5] & data[i,5] <= quantile(data[,5],.75)) data[i,5] <- 3
    else data[i,5] <- 4
  }
  
  return(data)
}
