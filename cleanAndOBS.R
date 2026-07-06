# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   cleanAndOBS.R
# PURPOSE:  Compare at the population level how the Covariance and Precision matrix differes of the clean data 
#           differes with the Covariance matix and the Precision matrix of the contaminated data
# AUTHOR:   Sudhir T.
# DATE:     6th July, 2026 
# ==============================================================================
# 

rm(list = ls())
library(Tlasso)
library(tensr)
library(glasso)
library(expm)
library(rTensor)
library(doParallel)
library(furrr)



source("simulation.summary.R")
source("Lfunctions.R")
source("Separate.fit.correct.R")
source("Model.R")



# RUNs <- 100
RUNs <- 3
n <- 50
dimen <- c(45, 54)
# dimen <- c(30, 36, 30)
nvars <- prod(dimen)
# dimen <- c(110, 4)

K <- length(dimen)
tensor.order <- length(dimen)


d <- 1

diff.f.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
diff.f.omega <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
diff.f.sigma.max <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
diff.f.omega.max <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode




for(itr in 1:RUNs){
  print(itr)
  data <- Model(n, itr * 123456, dimen)
  x <- data[[1]]$x
  vax <- data[[1]]$vax
  Sigma <- data[[2]]
  Omega <- data[[3]]
}




