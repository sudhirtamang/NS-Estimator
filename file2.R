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
library(furrr)


source("C:/Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//Separate.fit.R")
source("C://Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//simulation.summary.R")
source("C:/Users/sudhi/Desktop/PhD projects/NS-Estimator/Lfunctions.R")
source("Lfunctions.R")



Model <- function(n, seed, dimen) {

  nvars <- prod(dimen) # number of variables
  K <- length(dimen) # order of X
  
  
  Omega <- purrr::map2(dimen, 1:K, \(x, y) ChainOmega(x, sd = y*100, norm.type = 2))
  # Omega[[2]] <- diag(dimen[[2]])
  Sigma <- purrr::map(Omega, \(x) solve(x))
  # Sigma <- purrr::map(Sigma, \(x) x/x[1, 1]) # covariance matrix
  # Omega <- purrr::map(Sigma, \(x) solve(x))
  dSigma <- purrr::map(Sigma, \(x) t(chol(x))) # square root of covariance matrix
  
  
  set.seed(seed) 
  
  # Generate data observation
  # training set
  vec_x <- matrix(rnorm(nvars * n), ncol = n) 
  x <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    d <- array(vec_x[, i], dimen)
    x <- do.call("[<-", c(list(x), rep(list(substitute()), K), list(i, atrans(d, dSigma))))
  }
  
  # validation set
  vec_vax <- matrix(rnorm(nvars * n), ncol = n) 
  vax <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    d <- array(vec_vax[, i], dimen)
    vax <- do.call("[<-", c(list(vax), rep(list(substitute()), K), list(i, atrans(d, dSigma))))
  }
  
  result <- list()
  result$x <- x
  result$vax <- vax
  
  return(list(result, Sigma, Omega))
}





Run <- 5
# dimen <- c(60, 2)
# dimen <- c(30, 36, 30)
# dimen <- c(30, 36, 30, 60)
dimen <- c(100, 4)
n <- 30
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

pctOut <- 0.1
# Run <- 100
for (run in 1:Run) { 
  print(run)
  # Generate training set and validation set
  data <- Model(n, run * 123456, dimen)
  x <- data[[1]]$x
  vax <- data[[1]]$vax
  Sigma <- data[[2]]
  Omega <- data[[3]]
  
  
  nOut <- ceiling(pctOut * n)
  idxcontami <- sample(n, nOut)
  contami <- rt(nOut*nvars, df=10)
  dim(contami) <- c(dimen, nOut)
  
  contamiData <- x[]
  for(i in seq_along(idxcontami)){
    contamiData[, , , idxcontami[i]] <- contami[, , , i]
  }
  
  TcontamiData <- NSEstimator(contamiData, dimen)
  
  
  vec_vax <- matrix(rnorm(nvars * n), ncol = n) 
  vax <- array(0, dim = c(dimen, n))
  for(i in 1:n) {
    vax[, , , i] <- atrans(array(vec_vax[, i], dimen), dSigma)
  }
  
  nOut <- ceiling(pctOut * n)
  idxcontami <- sample(n, nOut)
  contami <- rt(nOut*nvars, df=10)
  dim(contami) <- c(dimen, nOut)
  
  contamiData <- vax[]
  for(i in seq_along(idxcontami)){
    contamiData[, , , idxcontami[i]] <- contami[, , , i]
  }
  
  vax_TcontamiData <- NSEstimator(contamiData, dimen)
  # 
  # xtildeOmega <- tildeOmega(x, dimen, n)
  # xtildeOmega[[2]] <- diag(dimen[[2]])
  # xtilde_Sk <- purrr::map(1:K, \(k) tilde_Sk(x, xtildeOmega, dimen, k, n))
  # 
  # 
  
  # 

  # proper candidates of tuning parameters
  lamseq <- seq(0.15, 9, length.out = 400)
  # lamseq <- seq(0.0015, 0.1, length.out = 30)
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
  # out <- simulation.summary(list(fit$Omegahat[[1]]), list(Omega[[1]]), offdiag = FALSE)
  # av.error.f[run] <- out$av.error.f
  # av.error.max[run] <- out$av.error.max
  # av.tpr[run] <- out$av.tpr
  # av.tnr[run] <- out$av.tnr
  simulation.summary(list(fit$Omegahat[[1]]), list(Omega[[1]]), offdiag = FALSE)
  error.f[run, ] <- out$error.f
  error.max[run, ] <- out$error.max
  tpr[run, ] <- out$tpr
  tnr[run, ] <- out$tnr
  # 
  # Tx <- NSEstimator2(x, dimen)
  # Tvax <- NSEstimator2(vax, dimen)
  TxtildeOmega <- tildeOmega(Tx, dimen, n)
  TxtildeOmega[[2]] <- diag(dimen[[2]])
  Txtilde_Sk <- purrr::map(1:K, \(k) tilde_Sk(Tx, TxtildeOmega, dimen, k, n))
  # purrr::walk(TxtildeOmega, \(x) print(dim(x)))
  # purrr::walk(Txtilde_Sk, \(x) print(dim(x)))
  # purrr::walk(TxtildeOmega, \(x) print(x))
  # purrr::walk(Txtilde_Sk, \(x) print(x))
  
  # ====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  plan(multisession, workers = ceiling(availableCores() * .7))
  RHOs <- seq(-1, 1, 0.01)
  Grho <- vector("double", length(RHOs))
  B <- 10000
  func1 <- function(n, rho){
    tmp1 <- mvtnorm::rmvnorm(n, mean=c(0, 0), sigma=matrix(c(1, rho, rho, 1), nrow=2))
    mean(stats::ecdf(tmp1[, 1])(tmp1[, 1]) * stats::ecdf(tmp1[, 2])(tmp1[, 2]))
  }
  func2 <- function(rho, B, func1){
    results <- future_map_dbl(
      1:B, 
      \(idx) func1(rho = rho, n = n), # Only 'idx' is iterated
      .options = furrr_options(seed = 123)
    )
    mean(results)
  }
  Grho <- future_map_dbl(RHOs, func2, B = B, func1 = func1, .options = furrr_options(seed = 123))
  plot(RHOs, Grho, pch=19, cex=0.75)
  plan(sequential)

  
  
  corrected_Txtilde_Sk <- Txtilde_Sk[]
  for(i in 1:dimen[[1]]){
    for(j in 1:dimen[[1]]){
      tmp1 <- abs(Txtilde_Sk[[1]][i, j] - Grho)
      corrected_Txtilde_Sk[[1]][i, j] <- RHOs[[which.min(tmp1)]]
    }
  }
  
  corrected_Txtilde_Sk_FINAL <- corrected_Txtilde_Sk[]
  corrected_Txtilde_Sk_FINAL[[1]] <- as.matrix(nearPD(corrected_Txtilde_Sk[[1]], corr = FALSE)$mat)
  fit <- Separate.fit(Tx, Tvax, lambda.list = lambda.list)
  
  thres <- 1.0e-4
  maxit <- 1e4
  
  
  rho <- 0.15
  Out1 <- glasso(corrected_Txtilde_Sk_FINAL[[1]], rho = rho, penalize.diagonal = FALSE, maxit = maxit, thr = thres)
  hat_Omega <- as.matrix(Out1$wi)
  # normalization
  hat_Omega <- hat_Omega / norm(hat_Omega, type = "F")



  simulation.summary(list(hat_Omega), list(Omega[[1]]), offdiag = FALSE)
  
  
  
  
  plot(c(Txtilde_Sk[[1]]), c(corrected_Txtilde_Sk[[1]]), col="blue", cex=0.65,
       xlab="Original Estimate of Sigma1", ylab="", asp = 1)
  
  points(c(diag(Txtilde_Sk[[1]])), c(diag(corrected_Txtilde_Sk[[1]])),
         col="red", cex=0.65)
  
  points(c(Txtilde_Sk[[1]]), c(Txtilde_Sk[[1]]),
         col="green", cex=0.65)
  
  legend("topleft", 
         legend = c("Off-Diagonal", "Diagonal", "Original Estimates"), 
         col = c("blue", "red", "green"), 
         lty = 1,             # Type of line (1 = solid)
         lwd = 2,             # Line width
         bty = "o")
  
  
  saveRDS(corrected_Txtilde_Sk, file="C:\\Users\\sudhi\\Desktop\\PhD projects\\NS-Estimator\\corrected_Txtilde_Sk.rds")
  saveRDS(Txtilde_Sk, file="C:\\Users\\sudhi\\Desktop\\PhD projects\\NS-Estimator\\Txtilde_Sk.rds")

  Tfit <- Separate.fit(Tx, Tvax, lambda.list = lambda.list)
  Out1 <- QUIC::QUIC(corrected_Txtilde_Sk[[1]], Tfit$lambda[1],
       tol = 1e-4,
       maxIter = 1000, path=NULL)
  simulation.summary(list(Out1$X, Omega[[2]]), Omega, offdiag = FALSE)
  
  
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




