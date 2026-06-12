# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   compareSigma.R
# PURPOSE:  With reasonable number of replicates find the Frobenius norm of difference in estimated Sigma and True Sigma
# AUTHOR:   Sudhir T.
# DATE:     12th June, 2026 
# ==============================================================================

rm(list = ls())
library(Tlasso)
library(tensr)
library(glasso)
library(expm)
library(rTensor)
library(doParallel)
library(furrr)


source("C:/Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//Separate.fit.R")
source("C://Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//simulation.summary.R")
source("Lfunctions.R")
source("Separate.fit.correct.R")
source("Model.R")









RUNs <- 100
n <- 40
# RUNs <- 1
dimen <- c(45, 57)
K <- length(dimen)


# plan(multisession, workers = ceiling(availableCores() * .7))
# RHOs <- seq(-1, 1, 0.01)
# Grho <- vector("double", length(RHOs))
# B <- 10000
# # func1 <- function(n, rho){
# #   tmp1 <- mvtnorm::rmvnorm(n, mean=c(0, 0), sigma=matrix(c(1, rho, rho, 1), nrow=2))
# #   mean(qnorm((n/(n+1)) * stats::ecdf(tmp1[, 1])(tmp1[, 1]))   * qnorm((n/(n+1)) * stats::ecdf(tmp1[, 2])(tmp1[, 2]))  )
# # }
# func1 <- function(n, rho){
#   tmp1 <- mvtnorm::rmvnorm(n, mean=c(0, 0), sigma=matrix(c(1, rho, rho, 1), nrow=2))
#   mean(qnorm((1/(n+1)) * rank(tmp1[, 1]))   * qnorm((1/(n+1)) * rank(tmp1[, 2]))  )
# }
# func2 <- function(rho, B, func1, n){
#   results <- future_map_dbl(
#     1:B, 
#     \(idx) func1(rho = rho, n = n), # Only 'idx' is iterated
#     .options = furrr_options(seed = 123)
#   )
#   mean(results)
# }
# Grho <- future_map_dbl(RHOs, func2, B = B, func1 = func1, n = 50, .options = furrr_options(seed = 123))
# 
# plot(Grho)
# abline(0, 1)
# plan(sequential)


Fsfit <- vector("double", RUNs)
FTfit <- vector("double", RUNs)
FTfit1 <- vector("double", RUNs)
for(run in 1:RUNs){
  
  data <- Model(n, run * 123456, dimen)
  x <- data[[1]]$x
  vax <- data[[1]]$vax
  Sigma <- data[[2]]
  Omega <- data[[3]]
  Tx <- NSEstimator2(x, dimen)
  Tvax <- NSEstimator2(vax, dimen)
  
  
  
  # proper candidates of tuning parameters
  lamseq <- seq(1.5e-8, 0.15, length.out = 300)
  # lamseq <- seq(0.00182, 0.001824, length.out = 100)
  lambda.list <- list() # a list containing candidates of tuning parameters for each mode 
  for (i in 1:K) {
    lambda.list[[i]] <- lamseq
  }
  
  sfit <- Separate.fit.correct(x, vax, lambda.list = lambda.list)
  Tfit <- Separate.fit.correct(Tx, Tvax, lambda.list = lambda.list)
  Tfit1 <- Separate.fit.correct(Tx, Tvax, lambda.list = lambda.list, Grho=Grho)
  
  Fsfit[[run]] <- norm(sfit[["Omegahat"]][[1]][[2]] - Sigma[[1]], type="F")
  FTfit[[run]] <- norm(Tfit[["Omegahat"]][[1]][[2]] - Sigma[[1]], type="F")
  FTfit1[[run]] <- norm(Tfit1[["Omegahat"]][[1]][[2]] - Sigma[[1]], type="F")
}


cat("Mean Forb. norm difference of sample Est. and True Sigma for p1:", mean(Fsfit), "With Sample size of", n)
cat("Mean Forb. norm difference of transformed sample Est. and True Sigma for p1:", mean(FTfit), "With Sample size of", n)
cat("Mean Forb. norm difference of transformed and corrected Est. and True Sigma for p1:", mean(FTfit1), "With Sample size of", n)




