# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   compareSigmaOmega.R
# PURPOSE:  With reasonable number of replicates find the Frobenius norm of difference in estimated Sigma/Omega and True Sigma/Omega
#           Here, the simulation is done without considering any Omega as Identity
# AUTHOR:   Sudhir T.
# DATE:     25th June, 2026 
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


source("Separate.fit.R")
source("simulation.summary.R")
source("Lfunctions.R")
source("Separate.fit.correct.R")
source("Model.R")





pctOut <- 0.1
RUNs <- 100
n <- 50
dimen <- c(45, 54)
nvars <- prod(dimen)
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


d <- 1
# d <- K
lambda.x <- array(0, dim = c(RUNs, d))
est.sigma <- array(0, dim = c(RUNs, d))
est.omega <- array(0, dim = c(RUNs, d))
error.f <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
error.f.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
error.max <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
error.max.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
tpr <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
tnr <- array(0, dim = c(RUNs, d)) # true negative rate for each mode

# RUNs <- 1
d <- 1
# d <- K
contamina.lambda.x <- array(0, dim = c(RUNs, d))
contamina.est.sigma <- array(0, dim = c(RUNs, d))
contamina.est.omega <- array(0, dim = c(RUNs, d))
contamina.error.f <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contamina.error.f.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contamina.error.max <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contamina.error.max.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contamina.tpr <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
contamina.tnr <- array(0, dim = c(RUNs, d)) # true negative rate for each mode


d <- 1
# d <- K
contamina.lambda.Tx <- array(0, dim = c(RUNs, d))
contamina.est.sigma.T <- array(0, dim = c(RUNs, d))
contamina.est.omega.T <- array(0, dim = c(RUNs, d))
contamina.error.f.T <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
contamina.error.f.T.sigma <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
contamina.error.max.T <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
contamina.error.max.T.sigma <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
contamina.tpr.T <- array(0, dim = c(RUNs, d)) # averaged true positive rate
contamina.tnr.T <- array(0, dim = c(RUNs, d)) # averaged true negative rate

d <- 1
# d <- K
contamina.lambda.Tx.C <- array(0, dim = c(RUNs, d))
contamina.est.sigma.T.C <- array(0, dim = c(RUNs, d))
contamina.est.omega.T.C <- array(0, dim = c(RUNs, d))
contamina.error.f.T.C <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contamina.error.f.T.C.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contamina.error.max.T.C <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contamina.error.max.T.C.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contamina.tpr.T.C <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
contamina.tnr.T.C <- array(0, dim = c(RUNs, d)) # true negative rate for each mode

# d <- 1
# av.error.f <- array(0, dim = c(RUNss, d)) # averaged estimation error in Frobenius norm
# av.error.f.sigma <- array(0, dim = c(RUNs, d)) 
# av.error.max <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
# av.error.max.sigma <- array(0, dim = c(RUNs, d)) 
# av.tpr <- array(0, dim = c(RUNs, d)) # averaged true positive rate
# av.tnr <- array(0, dim = c(RUNs, d)) # averaged true negative rate
# 
# d <- K
# error.f.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
# error.max.sigma <- array(0, dim = c(RUNs, d)) 
# tpr <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
# tnr <- array(0, dim = c(RUNs, d)) # true negative rate for each mode
# f.sigma <- array(0, dim = c(RUNs, d))
# f.omega <- array(0, dim = c(RUNs, d))
# 
# # ==========>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# d <- 1
# av.error.f.Tx <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
# av.error.f.sigma.Tx <- array(0, dim = c(RUNs, d)) 
# av.error.max.Tx <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
# av.error.max.sigma.Tx <- array(0, dim = c(RUNs, d)) 
# av.tpr.Tx <- array(0, dim = c(RUNs, d)) # averaged true positive rate
# av.tnr.Tx <- array(0, dim = c(RUNs, d)) # averaged true negative rate
# 
# 
# d <- K
# error.f.sigma.Tx <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
# error.max.sigma.Tx <- array(0, dim = c(RUNs, d)) 
# tpr.Tx <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
# tnr.Tx <- array(0, dim = c(RUNs, d)) # true negative rate for each mode
# f.sigma.Tx <- array(0, dim = c(RUNs, d))
# f.omega.Tx <- array(0, dim = c(RUNs, d))
# 
# # ==========>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# d <- 1
# av.error.f.Tx.C <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
# av.error.f.sigma.Tx.C <- array(0, dim = c(RUNs, d)) 
# av.error.max.Tx.C <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
# av.error.max.sigma.Tx.C <- array(0, dim = c(RUNs, d)) 
# av.tpr.Tx.C <- array(0, dim = c(RUNs, d)) # averaged true positive rate
# av.tnr.Tx.C <- array(0, dim = c(RUNs, d)) # averaged true negative rate
# 
# d <- K
# error.f.sigma.Tx.C <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
# error.max.sigma.Tx.C <- array(0, dim = c(RUNs, d)) 
# tpr.Tx.C <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
# tnr.Tx.C <- array(0, dim = c(RUNs, d)) # true negative rate for each mode
# f.sigma.Tx.C <- array(0, dim = c(RUNs, d))
# f.omega.Tx.C <- array(0, dim = c(RUNs, d))

# ==========>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# proper candidates of tuning parameters
lamseq <- seq(1.5e-15, 1, length.out = 300)
lambda.list <- list() # a list containing candidates of tuning parameters for each mode
for (i in 1:K) {
  lambda.list[[i]] <- lamseq
}

lamseq.C <- seq(1.5e-15, 1, length.out = 300)
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
  
  nOut <- ceiling(pctOut * n)
  idxcontami <- sample(n, nOut)
  contami <- rt(nOut*nvars, df=10)
  dim(contami) <- c(dimen, nOut)
  
  contam.x <- x[]
  for(i in seq_along(idxcontami)){
    contam.x[, , idxcontami[i]] <- contami[, , i]
  }
  
  contam.Tx <- NSEstimator2(contam.x, dimen)
  
  
  idxcontami <- sample(n, nOut)
  contami <- rt(nOut*nvars, df=10)
  dim(contami) <- c(dimen, nOut)
  
  contam.vax <- vax[]
  for(i in seq_along(idxcontami)){
    contam.vax[, , idxcontami[i]] <- contami[, , i]
  }
  
  contam.Tvax <- NSEstimator2(contam.vax, dimen)
  
  Tx <- NSEstimator2(x, dimen)
  Tvax <- NSEstimator2(vax, dimen)
  
  fitx <- Separate.fit.correct(contam.x, contam.vax, lambda.list = lambda.list)
  fitTx <- Separate.fit.correct(contam.Tx, contam.Tvax, lambda.list = lambda.list)
  fitTx.C <- Separate.fit.correct(contam.Tx, contam.Tvax, lambda.list = lambda.list.C, Grho=Grho)
  # fitx <- Separate.fit.correct(contam.x, contam.vax, lambda.list = lambda.list, scale.vec = c(1, 0.2))
  # fitTx <- Separate.fit.correct(contam.Tx, contam.Tvax, lambda.list = lambda.list, scale.vec = c(1, 0.2))
  # fitTx.C <- Separate.fit.correct(contam.Tx, contam.Tvax, lambda.list = lambda.list.C, Grho=Grho, scale.vec = c(1, 0.2))
  
  
  
  outx <- simulation.summary(list(fitx$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  outx.sigma <- simulation.summary(list(fitx$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  error.f[itr, 1] <- outx$error.f
  error.max[itr, 1] <- outx$error.max
  tpr[itr, 1] <- outx$tpr
  tnr[itr, 1] <- outx$tnr
  est.sigma[itr, 1] <- norm(fitx$Omegahat[[1]][[2]], type="F")
  est.omega[itr, 1] <- norm(fitx$Omegahat[[1]][[1]], type="F")
  
  error.f.sigma[itr, 1] <- outx.sigma$error.f
  error.max.sigma[itr, 1] <- outx.sigma$error.max
  
  lambda.x[itr, 1] <- fitx$lambda[[1]]
  
  
  outTx <- simulation.summary(list(fitTx$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  outTx.sigma <- simulation.summary(list(fitTx$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  error.f.T[itr, 1] <- outTx$error.f
  error.max.T[itr, 1] <- outTx$error.max
  tpr.T[itr, 1] <- outTx$tpr
  tnr.T[itr, 1] <- outTx$tnr
  est.sigma.T[itr, 1] <- norm(fitTx$Omegahat[[1]][[2]], type="F")
  est.omega.T[itr, 1] <- norm(fitTx$Omegahat[[1]][[1]], type="F")
  
  error.f.T.sigma[itr, 1] <- outTx.sigma$error.f
  error.max.T.sigma[itr, 1] <- outTx.sigma$error.max
  
  lambda.Tx[itr, 1] <- fitTx$lambda[[1]]
  
  outTx.C <- simulation.summary(list(fitTx.C$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  outTx.C.sigma <- simulation.summary(list(fitTx.C$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  error.f.T.C[itr, 1] <- outTx.C$error.f
  error.max.T.C[itr, 1] <- outTx.C$error.max
  tpr.T.C[itr, 1] <- outTx.C$tpr
  tnr.T.C[itr, 1] <- outTx.C$tnr
  est.sigma.T.C[itr, 1] <- norm(fitTx.C$Omegahat[[1]][[2]], type="F")
  est.omega.T.C[itr, 1] <- norm(fitTx.C$Omegahat[[1]][[1]], type="F")
  
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


cat("Comparison for SIGMA::just sample, no transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(error.f.sigma), "SD:", sd(error.f.sigma), "\n")
cat("Mean Max. Error:", colMeans(error.max.sigma), "SD:", sd(error.max.sigma), "\n")
cat("Mean Frob. Norm Estimated Sigma:", colMeans(est.sigma), "SD:", sd(est.sigma), "\n")
cat("Frob. Norm True Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Mean Frob. Norm Estimated Omega:", colMeans(est.omega), "SD:", sd(est.omega), "\n")
cat("Frob. Norm True Omega:", norm(Omega[[1]], type="F"), "\n")



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
cat("Mean Frob. Norm Estimated Sigma:", colMeans(est.sigma.T), "SD:", sd(est.sigma.T), "\n")
cat("Frob. Norm True Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Mean Frob. Norm Estimated Omega:", colMeans(est.omega.T), "SD:", sd(est.omega.T), "\n")
cat("Frob. Norm True Omega:", norm(Omega[[1]], type="F"), "\n")



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
cat("Mean Frob. Norm Estimated Sigma:", colMeans(est.sigma.T.C), "SD:", sd(est.sigma.T.C), "\n")
cat("Frob. Norm True Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Mean Frob. Norm Estimated Omega:", colMeans(est.omega.T.C), "SD:", sd(est.omega.T.C), "\n")
cat("Frob. Norm True Omega:", norm(Omega[[1]], type="F"), "\n")



Here are XXXXXXXXXXXX
