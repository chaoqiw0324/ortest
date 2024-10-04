#' Regression or Super-Learning Based Conditional Independence Test
#' 
#' We test whether \code{x} and \code{y} are conditionally associated, given
#' \code{S} using a generalized linear model/ super learning method. 
#' 
#' @details All included variables should be either numeric or binary. If 
#' \code{S} includes less than 4 variables, regression model is fitted. If
#' \code{y} is numeric, a linear regression model is fitted. If \code{y} is 
#' binary, a linear regression model is fitted. \code{x} and \code{S} are 
#' included as explanatory variables. If \code{S} includes more than 4 variables,
#' super learning method is applied. If \code{y} is numeric, the default method 
#' is linear regression model, mean response,MARS, random forest and XGBoost. If 
#' \code{y} is binary, the default method is LDA, mean response,MARS, random forest 
#' and XGBoost.
#' This model is tested whether \code{x} and \code{y} is independent conditional on
#' \code{S}. The final result is the test statistics and standard error.
#'
#' @return Two numeric, which are the test statistics and standard error of the test. 
#' 
#' @export

all_numeric <- function(df) {
  all_numeric <- all(sapply(df, is.numeric))
  return(all_numeric)
}

psi.hat_linear <- function(y, x, S=c(), subset = NULL, out.bin = TRUE, exp.bin = FALSE, root = "uni"){
  ## Function: estimate the odds ratio parameter psi
  ## Input: 1. An outcome nuisance model, onm = f(y|L,x=0)
  ##        2. An exposure nuisance model, enm = g(x|y=0,L)
  ## Output: A real number (vector) psi.hat, an estimate of the conditional odds ratio parameter psi
  
  
  if(length(S)==0){
    fm_out <- "y ~ x"
    fm_out <- as.formula(fm_out)
    fm_exp <- "x ~ y"
    fm_exp <- as.formula(fm_exp) ## for cases where conditioning set is empty
    dat <- data.frame(y,x)
  }else{
    covnames <- colnames(S)
    fm_out <- paste0("y ~ x + ", paste(covnames, collapse = "+"))
    fm_out <- as.formula(fm_out)
    fm_exp <- paste0("x ~ y + ", paste(covnames, collapse = "+"))
    fm_exp <- as.formula(fm_exp)
    dat <- data.frame(y,x,S)
  }
  
  
  if(!is.null(subset)) y <- y[subset]
  if(!is.null(subset)) x <- x[subset,]
  if(!is.null(subset)) S <- S[subset,]
  #temp
  
  
  refx <- 0 
  refy <- 0
  
  if(out.bin && exp.bin){
    h.dag <- 0.25 ## probability f.dag(y|S) = g.dag(x|S) = 0.5, i.e., y ~ x ~ Bernoulli(0.5)
    dat1 <- dat
    dat2 <- dat
    
    outcome <- glm(fm_out,family=binomial,data = dat)
    dat1$x <- 0 ## setting x=0
    onm <- predict.glm(outcome, newdata=dat1,type="response") # may cause warning: prediction from a rank-deficient fit may be misleading
    onm[y==0] <- 1-onm[y==0]
    
    exposure <- glm(fm_exp, family = binomial,data = dat)
    dat2$y <- 0 ## setting y=0
    enm <- predict.glm(exposure, newdata=dat2,type="response") # may cause warning: prediction from a rank-deficient fit may be misleading
    enm[x==0] <- 1-enm[x==0]
    
    d.diff <- (-1)^(y+x) ## Eric's suggestion, for y,x binary
    
    # build estimate function
    estimating_function <- function(psi){
      estf = d.diff*h.dag / (exp(psi*y*x)*onm*enm)
      return(estf) 
    }
    
    estimating_equation <- function(psi){
      estf = estimating_function(psi)         
      este = sum(estf)                       
      return(este)
    }
    
    
    res <- tryCatch({
      if (root == "uni") {
        U <- function(psi,onm,enm,d.diff){ sum( d.diff*h.dag / (exp(psi*y*x)*onm*enm) ) }
        est <- uniroot(U, interval = c(-3.0, 3.0), extendInt = "yes", tol = 0.001,maxiter=1000, onm = onm, enm = enm, d.diff=d.diff)
        res <-  est$root
      }else if(root == "multi"){
        proc <- rootSolve::multiroot(f = estimating_equation,     
                                     start = c(-3.0))
        res <- proc$root
      }
    }, error = function(e) {
      # If an error occurs, set res to zero
      cat("Error encountered: ", conditionMessage(e), "\n")
      0
    })
    
    # Baking the bread (approximate derivative)
    deriv <- numDeriv::jacobian(func = estimating_equation,   
                                x = res)              
    bread <- -1*deriv / n
    
    # Cooking the filling (matrix algebra)
    outerprod <- sum(estimating_function(res) * estimating_function(res)) 
    # alternative code using matrix algebra
    # outerprod <- t(estimating_function(mu_root)) %*% estimating_function(mu_root) 
    filling <- outerprod/n 
    
    # Assembling the sandwich (matrix algebra)
    sandwich <- (bread^-1) %*% filling %*% (bread^-1)
    se <- as.numeric(sqrt(sandwich / n))
    
    result <- c(res,se)
    return(result)
  }
  
  if(!out.bin){ 
    outcome <- glm(fm_out, family = gaussian,data = dat)
  } else outcome <- glm(fm_out, family = binomial,data = dat)
  dat1 <- dat
  dat1$x <- 0 ## setting x=0
  onm <- predict.glm(outcome, newdata=dat1,type="response")
  
  
  if(!exp.bin){
    exposure <- glm(fm_exp, family = gaussian,data = dat)
  } else exposure <- glm(fm_exp, family = binomial,data = dat)
  dat2 <- dat
  dat2$y <- 0 ## setting y=0
  enm <- predict.glm(exposure, newdata=dat2,type="response")
  
  
  
  #U <- function(psi,onm,enm){ sum( (y - onm)*(x - enm)*exp(-psi*y*x) )}
  #est <- uniroot(U, interval = c(-3.0, 3.0), extendInt = "yes", tol = 0.001, onm = onm, enm = enm)
  #return(est$root)
  
  # build estimate function
  estimating_function <- function(psi){
    estf = (y - onm)*(x - enm)*exp(-psi*y*x) 
    return(estf) 
  }
  
  estimating_equation <- function(psi){
    estf = estimating_function(psi)         
    este = sum(estf)                       
    return(este)
  }   
  
  if (root == "uni") {
    U <- function(psi){ sum( (y - onm)*(x - enm)*exp(-psi*y*x) )}
    par.init <- c(0)
    sol <- BBsolve(par=par.init, fn=U, quiet=TRUE) 
    res <- sol$par
  }else if(root == "multi"){
    proc <- rootSolve::multiroot(f = estimating_equation,     
                                 start = c(-3.0))
    res <- proc$root
  }
  
  # Baking the bread (approximate derivative)
  deriv <- numDeriv::jacobian(func = estimating_equation,   
                              x = res)              
  bread <- -1*deriv / n
  
  # Cooking the filling (matrix algebra)
  outerprod <- sum(estimating_function(res) * estimating_function(res)) 
  # alternative code using matrix algebra
  # outerprod <- t(estimating_function(mu_root)) %*% estimating_function(mu_root) 
  filling <- outerprod/n 
  
  # Assembling the sandwich (matrix algebra)
  sandwich <- (bread^-1) %*% filling %*% (bread^-1)
  se <- as.numeric(sqrt(sandwich / n))
  
  result <- c(res,se)
  return(result)
  
}

psi.hat_sl_2in1 <- function(y, x, S=c(), subset = NULL, out.bin = TRUE, exp.bin = FALSE,  root = "uni",sl=NULL,cross_fitting = FALSE,kfolds=5){
  ## Function: estimate the odds ratio parameter psi
  ## Input: 1. An outcome nuisance model, onm = f(y|S,x=0)
  ##        2. An exposure nuisance model, enm = g(x|y=0,S)
  ## Output: A real number (vector) psi.hat_ranger, an estimate of the conditional odds ratio parameter psi
  
  ## Todo: generalize estimating eq to arbitrary dimensions
  
  X_exp <- data.frame(y, S)
  Y_exp <- data.frame(x)
  X_out <- data.frame(x, S)
  Y_out <- data.frame(y)
  
  
  if (!is.null(subset)) y <- y[subset]
  if (!is.null(subset)) x <- x[subset,]
  if (!is.null(subset)) S <- S[subset,]
  #temp
  
  
  refx <- 0 
  refy <- 0
  
  
  allowed_methods <- c(
    "SL.bartMachine", "SL.bayesglm", "SL.biglasso", "SL.caret", "SL.caret.rpart", 
    "SL.cforest", "SL.earth", "SL.gam", "SL.gbm", "SL.glm", "SL.glm.interaction", 
    "SL.glmnet", "SL.ipredbagg", "SL.kernelKnn", "SL.knn", "SL.ksvm", "SL.lda", 
    "SL.leekasso", "SL.lm", "SL.loess", "SL.logreg", "SL.mean", "SL.nnet", 
    "SL.nnls", "SL.polymars", "SL.qda", "SL.randomForest", "SL.ranger", 
    "SL.ridge", "SL.rpart", "SL.rpartPrune", "SL.speedglm", "SL.speedlm", 
    "SL.step", "SL.step.forward", "SL.step.interaction", "SL.stepAIC", "SL.svm", 
    "SL.template", "SL.xgboost"
  )
  
  # Check if all elements in 'sl' exist in allowed methods
  if (!all(sl %in% allowed_methods)) {
    stop("Some elements in 'sl' do not exist in the list of allowed methods.")
  }
  sl_bin <- c("SL.lda","SL.mean","SL.earth", "SL.ranger","SL.xgboost")
  sl_gau <- c("SL.glm","SL.mean","SL.earth", "SL.ranger","SL.xgboost")
  if(! is.null(sl)){
    # Combine 'sl' with 'sl_existing'
    sl_binomial <- sl
    sl_gaussian <- sl
  }else{
    sl_binomial <- unique(c(sl_bin, sl))
    sl_gaussian <- unique(c(sl_gau, sl))
  }
  
  if (cross_fitting==FALSE) {
    ################################ No Cross_Fitting
    
    
    if(out.bin && exp.bin){
      h.dag <- 0.25 ## probability f.dag(y|S) = g.dag(x|S) = 0.5, i.e., y ~ x ~ Bernoulli(0.5)
      # outcome 
      ## outcome tune parameter 
      
      X_out_new <- X_out
      Y_out1 <- Y_out
      outcome <- SuperLearner(Y =Y_out1$y,
                              X = X_out,
                              family = binomial(),
                              SL.library = sl_binomial)
      
      X_out_new$x <- 0 ## setting x=0
      onm <- predict(outcome, X=X_out_new)$pred
      # onm <- predict(outcome, data=dat1,type="prob")[,2]
      onm[y==0] <- 1-onm[y==0]
      
      # exposure 
      # exposure tune parameter
      
      X_exp_new <- X_exp
      Y_exp1 <- Y_exp
      
      exposure <- SuperLearner(Y = Y_exp1$x,
                               X = X_exp,
                               family = binomial(),
                               SL.library = sl_binomial)
      X_exp_new$y <- 0 ## setting y=0
      enm <- predict(exposure, X = X_exp_new)$pred
      # enm <- predict(exposure_train, data=dat2,type="prob")[,2]
      enm[x==0] <- 1-enm[x==0]
      
      
      d.diff <- (-1)^(y+x) ## Eric's suggestion, for y,x binary
      # build estimate function
      estimating_function <- function(psi){
        estf = d.diff*h.dag / (exp(psi*y*x)*onm*enm)
        return(estf) 
      }
      
      estimating_equation <- function(psi){
        estf = estimating_function(psi)         
        este = sum(estf)                       
        return(este)
      }
      
      
      res <- tryCatch({
        if (root == "uni") {
          U <- function(psi,onm,enm,d.diff){ sum( d.diff*h.dag / (exp(psi*y*x)*onm*enm) ) }
          est <- uniroot(U, interval = c(-3.0, 3.0), extendInt = "yes", tol = 0.001,maxiter=1000, onm = onm, enm = enm, d.diff=d.diff)
          res <-  est$root
        }else if(root == "multi"){
          proc <- rootSolve::multiroot(f = estimating_equation,     
                                       start = c(-3.0))
          res <- proc$root
        }
      }, error = function(e) {
        # If an error occurs, set res to zero
        cat("Error encountered: ", conditionMessage(e), "\n")
        0
      })      
      
      # Baking the bread (approximate derivative)
      deriv <- numDeriv::jacobian(func = estimating_equation,   
                                  x = res)              
      bread <- -1*deriv / n
      
      # Cooking the filling (matrix algebra)
      outerprod <- sum(estimating_function(res) * estimating_function(res)) 
      # alternative code using matrix algebra
      # outerprod <- t(estimating_function(mu_root)) %*% estimating_function(mu_root) 
      filling <- outerprod/n 
      
      # Assembling the sandwich (matrix algebra)
      sandwich <- (bread^-1) %*% filling %*% (bread^-1)
      se <- as.numeric(sqrt(sandwich / n))
      
      result <- c(res,se)
      return(result)
    }
    
    
    
    
    # at least one variable is non-bianry  
    
    if(!out.bin){
      
      X_out_new <- X_out
      Y_out$y1 <- (Y_out$y-mean(Y_out$y))/sd(Y_out$y)
      
      
      outcome <- SuperLearner(Y = Y_out$y1,
                              X = X_out,
                              family = gaussian(),
                              SL.library = sl_gaussian)
      # onm <- predict(outcome, data=dat1,type="response")$prediction
      onm <- predict(outcome, X = X_out_new)$pred
      onm <- onm*sd(Y_out$y)+mean(Y_out$y)
    }else{
      
      X_out_new <- X_out
      Y_out1 <- Y_out
      outcome <- SuperLearner(Y =Y_out1$y,
                              X = X_out,
                              family = binomial(),
                              SL.library = sl_binomial)
      
      X_out_new$x <- 0 ## setting x=0
      onm <- predict(outcome, X = X_out_new)$pred
      # onm <- predict(outcome, data=dat1,type="prob")[,2]
      onm[y == 0] <- 1 - onm[y == 0]
    }
    
    
    if(!exp.bin){
      X_exp_new <- X_exp
      Y_exp$x1 <- (Y_exp$x-mean(Y_exp$x))/sd(Y_exp$x)
      
      
      
      exposure <- SuperLearner(Y = Y_exp$x1,
                               X = X_exp,
                               family = gaussian(),
                               SL.library = sl_gaussian)
      X_exp_new$y <- 0 ## setting y=0
      # enm <- predict(exposure, data=dat2,type="response")$predictions
      enm <- predict(exposure, X = X_exp_new)$pred
      enm <- enm*sd(Y_exp$x)+mean(Y_exp$x)
    }else{
      X_exp_new <- X_exp
      Y_exp1 <- Y_exp
      
      exposure <- SuperLearner(Y = Y_exp1$x,
                               X = X_exp,
                               family = binomial(),
                               SL.library = sl_binomial)
      X_exp_new$y <- 0 ## setting Y=0
      enm <- predict(exposure, X=X_exp_new)$pred
      # enm <- predict(exposure_train, data=dat2,type="prob")[,2]
      enm[x == 0] <- 1 - enm[x == 0]
    }
    
    # build estimate function
    estimating_function <- function(psi){
      estf = (y - onm)*(x - enm)*exp(-psi*y*x) 
      return(estf) 
    }
    
    estimating_equation <- function(psi){
      estf = estimating_function(psi)         
      este = sum(estf)                       
      return(este)
    }   
    
    
    res <- tryCatch({
      if (root == "uni") {
        U <- function(psi){ sum( (y - onm)*(x - enm)*exp(-psi*y*x) )}
        par.init <- c(0)
        sol <- BBsolve(par=par.init, fn=U, quiet=TRUE) 
        res <- sol$par
      }else if(root == "multi"){
        proc <- rootSolve::multiroot(f = estimating_equation,     
                                     start = c(-3.0))
        res <- proc$root
      }
    }, error = function(e) {
      # If an error occurs, set res to zero
      cat("Error encountered: ", conditionMessage(e), "\n")
      0
    })    
    # Baking the bread (approximate derivative)
    deriv <- numDeriv::jacobian(func = estimating_equation,   
                                x = res)              
    bread <- -1*deriv / n
    
    # Cooking the filling (matrix algebra)
    outerprod <- sum(estimating_function(res) * estimating_function(res)) 
    # alternative code using matrix algebra
    # outerprod <- t(estimating_function(mu_root)) %*% estimating_function(mu_root) 
    filling <- outerprod/n 
    
    # Assembling the sandwich (matrix algebra)
    sandwich <- (bread^-1) %*% filling %*% (bread^-1)
    se <- as.numeric(sqrt(sandwich / n))
    
    result <- c(res,se)
    return(result)  
    
    
    
  }else{
    ################################ Cross_Fitting   
    set.seed(2024)
    folds <- createFolds(dat$y, k = kfolds, list = TRUE, returnTrain = FALSE)
    predictions <- vector("list", kfolds)
    
    for (i in 1:kfolds) {
      # Split the data into training and test sets
      test_indices <- folds[[i]]
      train_indices <- setdiff(1:nrow(dat), test_indices)
      
      if(out.bin && exp.bin){
        h.dag <- 0.25 ## probability f.dag(y|S) = g.dag(x|S) = 0.5, i.e., y ~ x ~ Bernoulli(0.5)
        X_out_train  <- as_tibble(X_out[train_indices,])
        colnames(X_out_train) <- colnames(X_out)
        
        X_out_test   <- as_tibble(X_out[test_indices,])
        colnames(X_out_test) <- colnames(X_out)
        
        Y_out_train  <- as_tibble(Y_out[train_indices,])
        colnames(Y_out_train) <- colnames(Y_out)
        
        Y_out_test   <- as_tibble(Y_out[test_indices,])
        colnames(Y_out_test) <- colnames(Y_out)
        
        outcome <- SuperLearner(Y =Y_out_train$y,
                                X = X_out_train,
                                family = binomial(),
                                SL.library = sl_binomial)
        
        X_out_test$x <- 0 ## setting x=0
        onm <- predict(outcome, newdata=X_out_test)$pred
        # onm <- predict(outcome, data=dat1,type="prob")[,2]
        
        # exposure 
        # exposure tune parameter
        X_exp_train  <- as_tibble(X_exp[train_indices,])
        colnames(X_exp_train) <- colnames(X_exp)
        
        X_exp_test   <- as_tibble(X_exp[test_indices,])
        colnames(X_exp_test) <- colnames(X_exp)
        
        Y_exp_train  <- as_tibble(Y_exp[train_indices,])
        colnames(Y_exp_train) <- colnames(Y_exp)
        
        Y_exp_test   <- as_tibble(Y_exp[test_indices,])
        colnames(Y_exp_test) <- colnames(Y_exp)
        
        exposure <- SuperLearner(Y = Y_exp_train$x,
                                 X = X_exp_train,
                                 family = binomial(),
                                 SL.library = sl_binomial)
        
        
        X_exp_test$y <- 0 ## setting y=0
        enm <- predict(exposure, newdata = X_exp_test)$pred
        # enm <- predict(exposure_train, data=dat2,type="prob")[,2]
        
        d.diff <- (-1)^(Y_out_test$y+Y_exp_test$x) ## Eric's suggestion, for y,x binary
        
        # build estimate function
        estimating_function <- function(psi){
          estf = d.diff*h.dag / (exp(psi*Y_out_test$y*Y_exp_test$x)*onm*enm)
          return(estf)
        }
        
        estimating_equation <- function(psi){
          estf = estimating_function(psi)
          este = sum(estf)
          return(este)
        }
        res <- tryCatch({
          if (root == "uni") {
            # U <- function(psi, onm, enm, d.diff) { 
            #   sum(d.diff * h.dag / (exp(psi * Y_out_test$y * Y_exp_test$x) * onm * enm)) 
            # }
            # est <- uniroot(U, interval = c(-3.0, 3.0), extendInt = "yes", tol = 0.001, maxiter = 1000, onm = onm, enm = enm, d.diff = d.diff)
            # est$root
            U <- function(psi) { 
              sum(d.diff * h.dag / (exp(psi * Y_out_test$y * Y_exp_test$x) * onm * enm)) 
              
            }
            par.init <- c(0)
            sol <- BBsolve(par=par.init, fn=U, quiet=TRUE) 
            sol$par
          } else if (root == "multi") {
            proc <- rootSolve::multiroot(f = estimating_equation, start = c(-3.0))
            proc$root
          }
        }, error = function(e) {
          # If an error occurs, set res to zero
          cat("Error encountered: ", conditionMessage(e), "\n")
          0
        })
        
        # Baking the bread (approximate derivative)
        deriv <- numDeriv::jacobian(func = estimating_equation,
                                    x = res)
        bread <- -1*deriv / (n/kfolds)
        
        # Cooking the filling (matrix algebra)
        outerprod <- sum(estimating_function(res) * estimating_function(res))
        # alternative code using matrix algebra
        # outerprod <- t(estimating_function(mu_root)) %*% estimating_function(mu_root)
        filling <- outerprod/(n/kfolds)
        
        # Assembling the sandwich (matrix algebra)
        sandwich <- (bread^-1) %*% filling %*% (bread^-1)
        se <- as.numeric(sqrt(sandwich / (n/kfolds))) 
        result <- c(res,se)
        
        predictions[[i]] <- data.frame(
          result = res,
          se = se
        )
        next
      }
      
      if(!out.bin){
        X_out_train  <- as_tibble(X_out[train_indices,])
        colnames(X_out_train) <- colnames(X_out)
        
        X_out_test   <- as_tibble(X_out[test_indices,])
        colnames(X_out_test) <- colnames(X_out)
        
        Y_out_train  <- as_tibble(Y_out[train_indices,])
        colnames(Y_out_train) <- colnames(Y_out)
        
        Y_out_test   <- as_tibble(Y_out[test_indices,])
        colnames(Y_out_test) <- colnames(Y_out)
        
        Y_out_train$y1 <- Y_out_train$y
        outcome <- SuperLearner(Y = Y_out_train$y1,
                                X = X_out_train,
                                family = gaussian(),
                                SL.library = sl_gaussian)
        
        
        X_out_test$x <- 0 ## setting x=0
        # onm <- predict(outcome, data=dat1,type="response")$prediction
        onm <- predict(outcome, newdata = X_out_test)$pred
      }else{
        X_out_train  <- as_tibble(X_out[train_indices,])
        colnames(X_out_train) <- colnames(X_out)
        
        X_out_test   <- as_tibble(X_out[test_indices,])
        colnames(X_out_test) <- colnames(X_out)
        
        Y_out_train  <- as_tibble(Y_out[train_indices,])
        colnames(Y_out_train) <- colnames(Y_out)
        
        Y_out_test   <- as_tibble(Y_out[test_indices,])
        colnames(Y_out_test) <- colnames(Y_out)
        
        outcome <- SuperLearner(Y =Y_out_train$y,
                                X = X_out_train,
                                family = binomial(),
                                SL.library = sl_binomial)
        
        
        X_out_test$x <- 0 ## setting x=0
        onm <- predict(outcome, newdata = X_out_test)$pred
        # onm <- predict(outcome, data=dat1,type="prob")[,2]
      }
      
      
      if(!exp.bin){
        X_exp_train  <- as_tibble(X_exp[train_indices,])
        colnames(X_exp_train) <- colnames(X_exp)
        
        X_exp_test   <- as_tibble(X_exp[test_indices,])
        colnames(X_exp_test) <- colnames(X_exp)
        
        Y_exp_train  <- as_tibble(Y_exp[train_indices,])
        colnames(Y_exp_train) <- colnames(Y_exp)
        
        Y_exp_test   <- as_tibble(Y_exp[test_indices,])
        colnames(Y_exp_test) <- colnames(Y_exp)
        
        Y_exp_train$x1 <- Y_exp_train$x
        exposure <- SuperLearner(Y = Y_exp_train$x1,
                                 X = X_exp_train,
                                 family = gaussian(),
                                 SL.library = sl_gaussian)
        
        X_exp_test$y <- 0 ## setting y=0
        # enm <- predict(exposure, data=dat2,type="response")$predictions
        enm <- predict(exposure, newdata = X_exp_test)$pred
        
      }else{
        X_exp_train  <- as_tibble(X_exp[train_indices,])
        colnames(X_exp_train) <- colnames(X_exp)
        
        X_exp_test   <- as_tibble(X_exp[test_indices,])
        colnames(X_exp_test) <- colnames(X_exp)
        
        Y_exp_train  <- as_tibble(Y_exp[train_indices,])
        colnames(Y_exp_train) <- colnames(Y_exp)
        
        Y_exp_test   <- as_tibble(Y_exp[test_indices,])
        colnames(Y_exp_test) <- colnames(Y_exp)
        
        exposure <- SuperLearner(Y = Y_exp_train$x,
                                 X = X_exp_train,
                                 family = binomial(),
                                 SL.library = sl_binomial)
        
        
        X_exp_test$y <- 0 ## setting y=0
        enm <- predict(exposure, newdata = X_exp_test)$pred
        
      }
      
      
      # build estimate function
      estimating_function <- function(psi){
        estf = (Y_out_test$y - onm)*(Y_exp_test$x - enm)*exp(-psi*Y_out_test$y*Y_exp_test$x)
        return(estf)
      }
      
      estimating_equation <- function(psi){
        estf = estimating_function(psi)
        este = sum(estf)
        return(este)
      }
      
      res <- tryCatch({
        if (root == "uni") {
          U <- function(psi){ sum( (Y_out_test$y - onm)*(Y_exp_test$x - enm)*exp(-psi*Y_out_test$y*Y_exp_test$x) )}
          par.init <- c(0)
          sol <- BBsolve(par=par.init, fn=U, quiet=TRUE) 
          sol$par
        } else if (root == "multi") {
          proc <- rootSolve::multiroot(f = estimating_equation, start = c(-3.0))
          proc$root
        }
      }, error = function(e) {
        # If an error occurs, set res to zero
        cat("Error encountered: ", conditionMessage(e), "\n")
        0
      })
      
      # Baking the bread (approximate derivative)
      deriv <- numDeriv::jacobian(func = estimating_equation,
                                  x = res)
      bread <- -1*deriv / (n/kfolds)
      
      # Cooking the filling (matrix algebra)
      outerprod <- sum(estimating_function(res) * estimating_function(res))
      # alternative code using matrix algebra
      # outerprod <- t(estimating_function(mu_root)) %*% estimating_function(mu_root)
      filling <- outerprod/(n/kfolds)
      
      # Assembling the sandwich (matrix algebra)
      sandwich <- (bread^-1) %*% filling %*% (bread^-1)
      se <- as.numeric(sqrt(sandwich /(n/kfolds)))
      
      predictions[[i]] <- data.frame(
        result = res,
        se = se
      )
    }# end of loop
    
    ## harmonize folds result and calculate the sd
    
    all_predictions <- do.call(rbind, predictions)
    res <- median(all_predictions$result)
    var <- median(all_predictions$se^2+(all_predictions$result-mean(all_predictions$result))^2)
    sd <- var^(1/2)
    return(c(res,sd))
    
  }  
  
  
}