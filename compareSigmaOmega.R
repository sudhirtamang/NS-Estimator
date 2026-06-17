# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   compareSigmaOmega.R
# PURPOSE:  With reasonable number of replicates find the Frobenius norm of difference in estimated Sigma/Omega and True Sigma/Omega
# AUTHOR:   Sudhir T.
# DATE:     15th June, 2026 
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


source("C:/Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//Separate.fit.R")
source("C://Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//simulation.summary.R")
source("Lfunctions.R")
source("Separate.fit.correct.R")
source("Model.R")






RUNs <- 100
n <- 50
dimen <- c(45, 200)
# dimen <- c(110, 4)

K <- length(dimen)
tensor.order <- length(dimen)


plan(multisession, workers = ceiling(availableCores() * .7))
RHOs <- seq(-1, 1, 0.01)
Grho <- vector("double", length(RHOs))
B <- 10000
func1 <- function(n, rho){
  tmp1 <- mvtnorm::rmvnorm(n, mean=c(0, 0), sigma=matrix(c(1, rho, rho, 1), nrow=2))
  mean(qnorm((1/(n+1)) * rank(tmp1[, 1]))   * qnorm((1/(n+1)) * rank(tmp1[, 2]))  )
}
func2 <- function(rho, B, func1, n){
  results <- future_map_dbl(
    1:B,
    \(idx) func1(rho = rho, n = n), # Only 'idx' is iterated
    .options = furrr_options(seed = 123)
  )
  mean(results)
}
Grho <- future_map_dbl(RHOs, func2, B = B, func1 = func1, n = 50, .options = furrr_options(seed = 123))
plan(sequential)



# RUNs <- 1
d <- 1
# d <- K
lambda.x <- array(0, dim = c(RUNs, d))
error.f <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
error.f.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
error.max <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
error.max.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
tpr <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
tnr <- array(0, dim = c(RUNs, d)) # true negative rate for each mode


d <- 1
# d <- K
lambda.Tx <- array(0, dim = c(RUNs, d))
error.f.T <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
error.f.T.sigma <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
error.max.T <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
error.max.T.sigma <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
tpr.T <- array(0, dim = c(RUNs, d)) # averaged true positive rate
tnr.T <- array(0, dim = c(RUNs, d)) # averaged true negative rate

d <- 1
# d <- K
lambda.Tx.C <- array(0, dim = c(RUNs, d))
error.f.T.C <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
error.f.T.C.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
error.max.T.C <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
error.max.T.C.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
tpr.T.C <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
tnr.T.C <- array(0, dim = c(RUNs, d)) # true negative rate for each mode


# proper candidates of tuning parameters
lamseq <- seq(1.5e-6, 1e-2, length.out = 300)
lambda.list <- list() # a list containing candidates of tuning parameters for each mode
for (i in 1:K) {
  lambda.list[[i]] <- lamseq
}

lamseq.C <- seq(1.5e-6, 1e-2, length.out = 300)
lambda.list.C <- list() # a list containing candidates of tuning parameters for each mode
for (i in 1:K) {
  lambda.list.C[[i]] <- lamseq.C
}

for(itr in 1:RUNs){
  
  print(itr)
  data <- Model(n, itr * 123456, dimen)
  x <- data[[1]]$x
  vax <- data[[1]]$vax
  Sigma <- data[[2]]
  Omega <- data[[3]]
  
  Tx <- NSEstimator2(x, dimen)
  Tvax <- NSEstimator2(vax, dimen)
  
  
  fitx <- Separate.fit.correct(x, vax, lambda.list = lambda.list)
  fitTx <- Separate.fit.correct(Tx, Tvax, lambda.list = lambda.list)
  fitTx.C <- Separate.fit.correct(Tx, Tvax, lambda.list = lambda.list.C, Grho=Grho)
  
  
  
  outx <- simulation.summary(list(fitx$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  outx.sigma <- simulation.summary(list(fitx$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  error.f[itr, 1] <- outx$error.f
  error.max[itr, 1] <- outx$error.max
  tpr[itr, 1] <- outx$tpr
  tnr[itr, 1] <- outx$tnr
  
  error.f.sigma[itr, 1] <- outx.sigma$error.f
  error.max.sigma[itr, 1] <- outx.sigma$error.max
  
  lambda.x[itr, 1] <- fitx$lambda[[1]]
  
  
  outTx <- simulation.summary(list(fitTx$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  outTx.sigma <- simulation.summary(list(fitTx$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  error.f.T[itr, 1] <- outTx$error.f
  error.max.T[itr, 1] <- outTx$error.max
  tpr.T[itr, 1] <- outTx$tpr
  tnr.T[itr, 1] <- outTx$tnr
  
  error.f.T.sigma[itr, 1] <- outTx.sigma$error.f
  error.max.T.sigma[itr, 1] <- outTx.sigma$error.max
  
  lambda.Tx[itr, 1] <- fitTx$lambda[[1]]
  
  outTx.C <- simulation.summary(list(fitTx.C$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  outTx.C.sigma <- simulation.summary(list(fitTx.C$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  error.f.T.C[itr, 1] <- outTx.C$error.f
  error.max.T.C[itr, 1] <- outTx.C$error.max
  tpr.T.C[itr, 1] <- outTx.C$tpr
  tnr.T.C[itr, 1] <- outTx.C$tnr
  
  error.f.T.C.sigma[itr, 1] <- outTx.C.sigma$error.f
  error.max.T.C.sigma[itr, 1] <- outTx.C.sigma$error.max
  
  lambda.Tx.C[itr, 1] <- fitTx.C$lambda[[1]]
}


cat("Comparison for OMEGA::just sample, no transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(error.f), "SD:", sd(error.f), "\n")
cat("Mean Max. Error:", colMeans(error.max), "SD:", sd(error.max), "\n")
cat("TPR:", colMeans(tpr), "SD:", sd(tpr), "\n")
cat("TNR:", colMeans(tnr), "SD:", sd(tnr), "\n")
cat("Comparison for SIGMA::just sample, no transformation, no correction", "\n")
cat("Range of lambda.best: ", sprintf("[%.10e, %.10e]\n", min(lambda.x), max(lambda.x)), "\n\n")

cat("Mean Forb. Diff.:", colMeans(error.f.sigma), "SD:", sd(error.f.sigma), "\n")
cat("Mean Max. Error:", colMeans(error.max.sigma), "SD:", sd(error.max.sigma), "\n")



cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")


cat("Comparison for OMEGA::with transformation, no correction", "\n")

cat("Mean Forb. Diff.:", colMeans(error.f.T), "SD:", sd(error.f.T), "\n")
cat("Mean Max. Error:", colMeans(error.max.T), "SD:", sd(error.max.T), "\n")
cat("TPR:", colMeans(tpr.T), "SD:", sd(tpr.T), "\n")
cat("TNR:", colMeans(tnr.T), "SD:", sd(tnr.T), "\n")
cat("Range of lambda.best: ", sprintf("[%.10e, %.10e]\n", min(lambda.Tx), max(lambda.Tx)), "\n\n")

cat("Comparison for SIGMA::with transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(error.f.T.sigma), "SD:", sd(error.f.T.sigma), "\n")
cat("Mean Max. Error:", colMeans(error.max.T.sigma), "SD:", sd(error.max.T.sigma), "\n")




cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")


cat("Comparison for OMEGA::with transformation with correction", "\n")
cat("Mean Forb. Diff.:", colMeans(error.f.T.C), "SD:", sd(error.f.T.C), "\n")
cat("Mean Max. Error:", colMeans(error.max.T.C), "SD:", sd(error.max.T.C), "\n")
cat("TPR:", colMeans(tpr.T.C), "SD:", sd(tpr.T.C), "\n")
cat("TNR:", colMeans(tnr.T.C), "SD:", sd(tnr.T.C), "\n")
cat("Range of lambda.best: ", sprintf("[%.10e, %.10e]\n", min(lambda.Tx.C), max(lambda.Tx.C)), "\n\n")

cat("Comparison for SIGMA::with transformation with correction", "\n")
cat("Mean Forb. Diff.:", colMeans(error.f.T.C.sigma), "SD:", sd(error.f.T.C.sigma), "\n")
cat("Mean Max. Error:", colMeans(error.max.T.C.sigma), "SD:", sd(error.max.T.C.sigma), "\n")





