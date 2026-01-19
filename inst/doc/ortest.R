## ----setup_options, include=FALSE---------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

## ----default, include = FALSE-------------------------------------------------
knitr::opts_chunk$set(comment = NA)

## ----setup,echo=FALSE,message=FALSE-------------------------------------------
# remotes::install_github("chaoqiw0324/ortest")
library(ortest)
library(pcalg)

## ----echo=FALSE---------------------------------------------------------------
n <- 1000

## ----eval=FALSE, echo=TRUE----------------------------------------------------
# L1 <- runif(n,0,1)
# L2 <- runif(n,0,1)
# L3 <- runif(n,0,1)
# L4 <- runif(n,0,1)
# L5 <- runif(n,0,1)
# L_vec <- tibble(
#   L1 = L1,
#   L2 = L2,
#   L3 = L3,
#   L4 = L4,
#   L5 = L5)
#   L.true <- 2*L1 + L2^2 + L1*L3 + 3*L4 + L1 * L5
# 
# Z <- 0.5 + 0.5*L.true
# pr <- 1/(1+8*exp(-Z))
# summary(pr)
# # control the probability close to 0.5
# Y.true <- 2*Z + rnorm(n,0,1)
# A.true <- rbinom(n,1,pr)
# dat <- tibble(
#   Y = Y.true,
#   A = A.true
# ) %>% cbind(L_vec)
# psi_hat_linear(y = dat$Y,x = dat$A,S = c(), subset = NULL, out_bin = FALSE, exp_bin = TRUE)
# psi_hat_linear_int(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = FALSE, exp_bin = TRUE,two_way = TRUE,three_way = TRUE)
# psi_hat_sl(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = FALSE, exp_bin = TRUE, sl=NULL,cross_fitting = FALSE,kfolds=5)
# psi_hat_sl(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = FALSE, exp_bin = TRUE,sl=NULL,cross_fitting = TRUE, kfolds=5)

## ----eval=FALSE---------------------------------------------------------------
# L1 <- runif(n,0,1)
# L2 <- runif(n,0,1)
# L3 <- runif(n,0,1)
# L4 <- runif(n,0,1)
# L5 <- runif(n,0,1)
# L_vec <- tibble(
#   L1 = L1,
#   L2 = L2,
#   L3 = L3,
#   L4 = L4,
#   L5 = L5)
# L.true <- 2*L1 + L2^2 + L1*L3 + 3*L4 + L1 * L5
# # L.true <- 2*L1 + 3*L2 + 4*L3
# Z <- 0.5 + 0.5*L.true
# Y.true <- 2*Z + rnorm(n,0,3)
# A.true <- 2*Z + rnorm(n,0,3)
# dat <- tibble(
#   Y = Y.true,
#   A = A.true
# ) %>% cbind(L_vec)
# 
# psi_hat_linear(y = dat$Y,x = dat$A,S = c(), subset = NULL, out_bin = FALSE, exp_bin = FALSE)
# psi_hat_linear_int(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = FALSE, exp_bin = FALSE,two_way = TRUE,three_way = TRUE)
# psi_hat_sl(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = FALSE, exp_bin = FALSE, sl=NULL,cross_fitting = FALSE,kfolds=5)
# psi_hat_sl(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = FALSE, exp_bin = FALSE, sl=NULL,cross_fitting = TRUE,kfolds=5)

## ----eval=FALSE---------------------------------------------------------------
#   L1 <- runif(n,0,1)
#   L2 <- runif(n,0,1)
#   L3 <- runif(n,0,1)
#   L4 <- runif(n,0,1)
#   L5 <- runif(n,0,1)
#   L_vec <- tibble(
#     L1 = L1,
#     L2 = L2,
#     L3 = L3,
#     L4 = L4,
#     L5 = L5)
#   L.true <- 2*L1 + L2^2 + L1*L3 + 3*L4 + L1 * L5
#   Z <- 0.5 + 0.5*L.true
#   pr <- 1/(1+8*exp(-Z))
#   summary(pr)
#   Y.true <- rbinom(n,1,pr)
#   A.true <- rbinom(n,1,pr)
#   dat <- tibble(
#     Y = Y.true,
#     A = A.true
#   ) %>% cbind(L_vec)
# 
# psi_hat_linear(y = dat$Y,x = dat$A,S = c(), subset = NULL, out_bin = TRUE, exp_bin = TRUE)
# psi_hat_linear_int(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = TRUE, exp_bin = TRUE,two_way = TRUE,three_way = TRUE)
# psi_hat_sl(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = TRUE, exp_bin = TRUE,sl=NULL,cross_fitting = FALSE,kfolds=5)
# psi_hat_sl(y = dat$Y,x = dat$A,S = L_vec, subset = NULL, out_bin = TRUE, exp_bin = TRUE,sl=NULL,cross_fitting = TRUE,kfolds=5)

## ----echo=FALSE---------------------------------------------------------------
# Generate the dataset
set.seed(2025)
data <- generate_test_data(N = 2000)

## ----eval=FALSE,message=FALSE-------------------------------------------------
# ## c to m
# ortest(1,5,c(2:4,6:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE))
# ortest(1,5,c(2:4,6:10),list(dat = data,method="sl", sl=NULL,cross_fitting=TRUE,kfolds=5,two_way=FALSE,three_way=FALSE))
# 
# ## m to b
# ortest(5,2,c(1,3:4,6:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE))
# ortest(5,2,c(1,3:4,6:10),list(dat = data,method="sl", sl=NULL,cross_fitting=TRUE,kfolds=5,two_way=FALSE,three_way=FALSE))
# 
# ## m to m
# ortest(5,6,c(1:4,7:10),list(dat = data,method="linear", sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE))
# ortest(5,6,c(1:4,7:10),list(dat = data,method="sl", sl=NULL,cross_fitting=TRUE,kfolds=5,two_way=FALSE,three_way=FALSE))
# 

## ----message=FALSE,warning=FALSE,error=FALSE----------------------------------
V <- colnames(data)
pc.fit <- pc(suffStat = list(dat = data,method="linear",sl=NULL,cross_fitting=FALSE,kfolds=5,two_way=FALSE,three_way=FALSE),
             indepTest = ortest, ## indep.test: partial correlations
             alpha=0.01, labels = V, verbose = FALSE)

## ----echo = F,message=FALSE,warning=FALSE,error=FALSE,fig.cap = "PC Plot showing the potential relationship"----
plot(pc.fit, main = "")

