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




Model <- function(n, seed, dimen) {

  nvars <- prod(dimen) # number of variables
  K <- length(dimen) # order of X
  
  
  Omega <- purrr::map2(dimen, 1:K, \(x, y) ChainOmega(x, sd = y*100, norm.type = 2))
  Omega[[2]] <- diag(dimen[[2]])
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
dimen <- c(60, 60)
n <- 40
K <- length(dimen)
run <- 1
# initialize measurements
# d <- 1
# av.error.f <- array(0, dim = c(Run, d)) # averaged estimation error in Frobenius norm
# av.error.max <- array(0, dim = c(Run, d)) # averaged estimation error in Maximum norm
# av.tpr <- array(0, dim = c(Run, d)) # averaged true positive rate
# av.tnr <- array(0, dim = c(Run, d)) # averaged true negative rate

d <- K
error.f <- array(0, dim = c(Run, d)) # estimation error in Frobenius norm for each mode
error.max <- array(0, dim = c(Run, d)) # estimation error in Maximum norm for each mode
tpr <- array(0, dim = c(Run, d)) # true positive rate for each mode
tnr <- array(0, dim = c(Run, d)) # true negative rate for each mode

# 
# d <- 1
# av.error.f.T <- array(0, dim = c(Run, d)) # averaged estimation error in Frobenius norm
# av.error.max.T <- array(0, dim = c(Run, d)) # averaged estimation error in Maximum norm
# av.tpr.T <- array(0, dim = c(Run, d)) # averaged true positive rate
# av.tnr.T <- array(0, dim = c(Run, d)) # averaged true negative rate

d <- K
error.f.T <- array(0, dim = c(Run, d)) # estimation error in Frobenius norm for each mode
error.max.T <- array(0, dim = c(Run, d)) # estimation error in Maximum norm for each mode
tpr.T <- array(0, dim = c(Run, d)) # true positive rate for each mode
tnr.T <- array(0, dim = c(Run, d)) # true negative rate for each mode


# Run <- 100
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
  lamseq <- seq(1e-09, 1e-04, length.out = 400)
  lambda.list <- list() # a list containing candidates of tuning parameters for each mode 
  for (i in 1:K) {
    lambda.list[[i]] <- lamseq
  }
  
  
  # Model fitting
  fit <- Separate.fit(x, vax, lambda.list = lambda.list)

  
  ## If there is no validation set, we can use cv.Separate to tune lambda via cross-validation
  # fit <- cv.Separate(x, lambda.list=lambda.list)
  ## With a user-specified sequence of lambdas in lambda.vec, we can directly fit the model
  # fit <- Separate.fit(x, lambda.vec = lambda.vec)
  
  # Simulation summary of estimation errors, TPR and TNR
  out <- simulation.summary(fit$Omegahat, Omega, offdiag = FALSE)
  # av.error.f[run] <- out$av.error.f
  # av.error.max[run] <- out$av.error.max
  # av.tpr[run] <- out$av.tpr
  # av.tnr[run] <- out$av.tnr
  
  error.f[run, ] <- out$error.f
  error.max[run, ] <- out$error.max
  tpr[run, ] <- out$tpr
  tnr[run, ] <- out$tnr

  
  TxtildeOmega <- tildeOmega(Tx, dimen, n)
  TxtildeOmega[[2]] <- diag(dimen[[2]])
  Txtilde_Sk <- purrr::map(1:K, \(k) tilde_Sk(Tx, TxtildeOmega, dimen, k, n))
  # purrr::walk(TxtildeOmega, \(x) print(dim(x)))
  # purrr::walk(Txtilde_Sk, \(x) print(dim(x)))
  # purrr::walk(TxtildeOmega, \(x) print(x))
  # purrr::walk(Txtilde_Sk, \(x) print(x))
  
  # ====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  RHOs <- seq(0, 1, 0.01)
  B <- 1000
  Grho <- vector("double", length(RHOs))
  for(j in seq_along(RHOs)){
    TOT1 <- 0
    for(i in 1:B){
      tmp1 <- mvtnorm::rmvnorm(n, mean=c(0, 0), sigma=matrix(c(1, RHOs[[j]], RHOs[[j]], 1), nrow=2))
      TOT1 <- TOT1 + mean(stats::ecdf(tmp1[, 1])(tmp1[, 1]) * stats::ecdf(tmp1[, 2])(tmp1[, 2]))
    }
    Grho[[j]] <- TOT1/B
  }
  
  corrected_Txtilde_Sk <- Txtilde_Sk[]
  for(i in 1:dimen[[1]]){
    for(j in 1:dimen[[1]]){
      tmp1 <- abs(Txtilde_Sk[[1]][i, j] - Grho)
      corrected_Txtilde_Sk[[1]][i, j] <- RHOs[[which.min(tmp1)]]
    }
  }
  
  
  Tfit <- Separate.fit(Tx, Tvax, lambda.list = lambda.list)
  
  Out1 = glasso(corrected_Txtilde_Sk[[1]], rho = Tfit$lambda[1], penalize.diagonal = FALSE, thr = 1.0e-4, maxit = 1e4)
  hat_Omega = as.matrix(Out1$wi)
  hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
  
  out <- simulation.summary(list(hat_Omega, Omega[[2]]), Omega, offdiag = FALSE)
  # av.error.f[run] <- out$av.error.f
  # av.error.max[run] <- out$av.error.max
  # av.tpr[run] <- out$av.tpr
  # av.tnr[run] <- out$av.tnr
  
  error.f.T[run, ] <- out$error.f
  error.max.T[run, ] <- out$error.max
  tpr.T[run, ] <- out$tpr
  tnr.T[run, ] <- out$tnr
  
  
  
  
  # ====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  

  
}

# estimation error
# mean(av.error.f)
colMeans(error.f)
# mean(av.error.max)
colMeans(error.max)

# TPR and TNR
# mean(av.tpr)
colMeans(tpr)
# mean(av.tnr)
colMeans(tnr)


# estimation error
# mean(av.error.f.T)
colMeans(error.f.T)
# mean(av.error.max.T)
colMeans(error.max.T)

# TPR and TNR
# mean(av.tpr.T)
colMeans(tpr.T)
# mean(av.tnr.T)
colMeans(tnr.T)




