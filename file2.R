# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   file2.R
# PURPOSE:  BIAS correction on the NS-Estimator
# AUTHOR:   Sudhir T.
# DATE:     12th May, 2026 
# ==============================================================================

rm(list = ls())
library(Tlasso)
library(tensr)
library(glasso)
library(expm)
library(rTensor)
library(doParallel)

source("C:/Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//Separate.fit.R")
source("C://Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//simulation.summary.R")
source("C:/Users/sudhi/Desktop/PhD projects/NS-Estimator/Lfunctions.R")
source("Lfunctions.R")


dimen <- c(120, 120)
n <- 20
K <- length(dimen)

Model <- function(n, seed, dimen) {

  nvars <- prod(dimen) # number of variables
  K <- length(dimen) # order of X
  
  
  Omega <- purrr::map2(dimen, 1:K, \(x, y) ChainOmega(x, sd = y*100, norm.type = 2))
  Sigma <- purrr::map(Omega, \(x) solve(x))
  Sigma <- purrr::map(Sigma, \(x) x/x[1, 1]) # covariance matrix
  dSigma <- purrr::map(Sigma, \(x) t(chol(x))) # square root of covariance matrix
  
  
  set.seed(seed) 
  
  # Generate data observation
  # training set
  vec_x <- matrix(rnorm(nvars * n), ncol = n) 
  x <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    x[, , i] <- array(vec_x[, i], dimen)
    x[, , i] <- atrans(x[, , i], dSigma)
  }
  
  # validation set
  vec_vax <- matrix(rnorm(nvars * n), ncol = n) 
  vax <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    vax[, , i] <- array(vec_vax[, i], dimen)
    vax[, , i] <- atrans(vax[, , i], dSigma)
  }
  
  result <- list()
  result$x <- x
  result$vax <- vax
  
  return(list(result, Sigma, Omega))
  
}




Run <- 5
# Run <- 100

# initialize measurements
d <- 1
av.error.f <- array(0, dim = c(Run, d)) # averaged estimation error in Frobenius norm
av.error.max <- array(0, dim = c(Run, d)) # averaged estimation error in Maximum norm
av.tpr <- array(0, dim = c(Run, d)) # averaged true positive rate
av.tnr <- array(0, dim = c(Run, d)) # averaged true negative rate

d <- K
error.f <- array(0, dim = c(Run, d)) # estimation error in Frobenius norm for each mode
error.max <- array(0, dim = c(Run, d)) # estimation error in Maximum norm for each mode
tpr <- array(0, dim = c(Run, d)) # true positive rate for each mode
tnr <- array(0, dim = c(Run, d)) # true negative rate for each mode


d <- 1
av.error.f.T <- array(0, dim = c(Run, d)) # averaged estimation error in Frobenius norm
av.error.max.T <- array(0, dim = c(Run, d)) # averaged estimation error in Maximum norm
av.tpr.T <- array(0, dim = c(Run, d)) # averaged true positive rate
av.tnr.T <- array(0, dim = c(Run, d)) # averaged true negative rate

d <- K
error.f.T <- array(0, dim = c(Run, d)) # estimation error in Frobenius norm for each mode
error.max.T <- array(0, dim = c(Run, d)) # estimation error in Maximum norm for each mode
tpr.T <- array(0, dim = c(Run, d)) # true positive rate for each mode
tnr.T <- array(0, dim = c(Run, d)) # true negative rate for each mode


for (run in 1:Run) { 
  print(run)
  # Generate training set and validation set
  data <- Model(n, run * 123456, dimen)
  x <- data[[1]]$x
  vax <- data[[1]]$vax
  Sigma <- data[[2]]
  Omega <- data[[3]]
  
  Tx <- NSEstimator2(x, dimen)
  Tvax <- NSEstimator2(vax, dimen)
  # proper candidates of tuning parameters
  lamseq <- seq(1e-09, 1e-04, length.out = 200)
  lambda.list <- list() # a list containing candidates of tuning parameters for each mode 
  for (i in 1:K) {
    lambda.list[[i]] <- lamseq
  }
  
  
  xtildeOmega <- tildeOmega(x, dimen, n)
  # Model fitting
  fit <- Separate.fit(x, vax, lambda.list = lambda.list)

  
  ## If there is no validation set, we can use cv.Separate to tune lambda via cross-validation
  # fit <- cv.Separate(x, lambda.list=lambda.list)
  ## With a user-specified sequence of lambdas in lambda.vec, we can directly fit the model
  # fit <- Separate.fit(x, lambda.vec = lambda.vec)
  
  # Simulation summary of estimation errors, TPR and TNR
  out <- simulation.summary(fit$Omegahat, Omega, offdiag = FALSE)
  av.error.f[run] <- out$av.error.f
  av.error.max[run] <- out$av.error.max
  av.tpr[run] <- out$av.tpr
  av.tnr[run] <- out$av.tnr
  
  error.f[run, ] <- out$error.f
  error.max[run, ] <- out$error.max
  tpr[run, ] <- out$tpr
  tnr[run, ] <- out$tnr
  
  
  fit.T <- Separate.fit(Tx, Tvax, lambda.list = lambda.list)
  out <- simulation.summary(fit.T$Omegahat, Omega, offdiag = FALSE)
  av.error.f.T[run] <- out$av.error.f
  av.error.max.T[run] <- out$av.error.max
  av.tpr.T[run] <- out$av.tpr
  av.tnr.T[run] <- out$av.tnr
  
  error.f.T[run, ] <- out$error.f
  error.max.T[run, ] <- out$error.max
  tpr.T[run, ] <- out$tpr
  tnr.T[run, ] <- out$tnr
  
}

# estimation error
mean(av.error.f)
colMeans(error.f)
mean(av.error.max)
colMeans(error.max)

# TPR and TNR
mean(av.tpr)
colMeans(tpr)
mean(av.tnr)
colMeans(tnr)


# estimation error
mean(av.error.f.T)
colMeans(error.f.T)
mean(av.error.max.T)
colMeans(error.max.T)

# TPR and TNR
mean(av.tpr.T)
colMeans(tpr.T)
mean(av.tnr.T)
colMeans(tnr.T)
