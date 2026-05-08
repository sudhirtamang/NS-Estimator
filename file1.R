# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   file1
# PURPOSE:  Simulate the performance of NS Estimator with clean data. Compare it with
#           purely fast and separable estimator
# AUTHOR:   Sudhir T
# DATE:     4th May, 2026 
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
source("C://Users//sudhi//Desktop//Fast and Separable Estimation Replication//replication//Model 1//Model7.R")
# source("C:/Users/sudhi/Desktop/PhD projects/NS-Estimator/Lfunctions.R")
source("Lfunctions.R")


Model7 <- function(n, seed) {
  
  # seed: random seed
  # n: sample size
  
  # Model setting
  dimen <- c(30, 36, 30) # dimension of X
  nvars <- prod(dimen) # number of variables
  K <- 3 # order of X
  
  
  Omega <- purrr::map2(dimen, 1:3, \(x, y) ChainOmega(x, sd = y*100, norm.type = 2))
  Sigma <- purrr::map(Omega, \(x) solve(x))
  scales <- purrr::map_dbl(Sigma, \(x) x[1, 1])
  Omega <- purrr::map2(Omega, Sigma, \(x, y) x*y[1, 1]) # precision matrix
  Sigma <- purrr::map(Sigma, \(x) x/x[1, 1]) # covariance matrix
  dSigma <- purrr::map(Sigma, \(x) t(chol(x))) # square root of covariance matrix

  
  set.seed(seed) 
  
  # Generate data observation
  # training set
  vec_x <- matrix(rnorm(nvars * n), ncol = n) 
  x <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    x[, , , i] <- array(vec_x[, i], dimen)
    x[, , , i] <- atrans(x[, , , i], dSigma)
  }
  
  # validation set
  vec_vax <- matrix(rnorm(nvars * n), ncol = n) 
  vax <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    vax[, , , i] <- array(vec_vax[, i], dimen)
    vax[, , , i] <- atrans(vax[, , , i], dSigma)
  }
  
  result <- list()
  result$x <- x
  result$vax <- vax
  
  return(list(result, Sigma, Omega, scales))
  
}



Separate.fit = function(x, val = NULL, est.mode = NULL, lambda.vec = NULL, lambda.list = NULL, Omegatilde.list = NULL, scale.vec = NULL, normalize = TRUE, thres = 1.0e-4, maxit = 1e4, njobs = 4) {
  
  if (is.null(est.mode) == TRUE) {
    est.mode = c(1:K)
  }
  if (is.null(val)) {
    if (is.null(lambda.vec) | length(lambda.vec) != length(est.mode)) {
      stop("lambda.vec is missing or does not have the correct length")
    }
  } else {
    if (is.null(lambda.list) | length(lambda.list) != length(est.mode)) {
      stop("lambda.list is missing or does not have the correct length")
    }
  }
  
  dimen = dim(x) # dimension of dataset
  K = length(dimen) - 1 # order of tensor
  n = dimen[K + 1] # sample size of training set
  n_val = dim(val)[K + 1] # sample size of validation set
  nvars = prod(dimen) # number of variables
  m.vec = dimen[1:K] # dimension of each observation
  
  if (!(is.null(Omegatilde.list) | length(Omegatilde.list) == K)) {
    stop("argument Omegatilde.list should be a list of M matrices")
  }
  
  if (is.null(scale.vec)) {
    scale.vec = rep(1, length(est.mode))
  }
  
  ##### Calculate \tilde\Omega #####
  if (is.null(Omegatilde.list) == FALSE) {
    fit1 = Omegatilde.list  # user-specified value for \tilde\Omega
  }else {
    # Calculate \tilde\Omega by the definition in the paper
    c1 = makeCluster(njobs)
    registerDoParallel(c1)
    fit1 = foreach(k = 1:K, .export = c("x"), .combine = list, .multicombine = TRUE) %dopar% {
      # when sample size is small, use the identity matrix
      if (n * nvars < ((dimen[k]**2) * (dimen[k] - 1) / 2)) {
        Omega_tilde = diag(dimen[k])
      }
      else {
        # when sample size is large, calculate the sample estimator of the precision matrices
        S.array = array(0, c(m.vec[k], m.vec[k], n))
        for (i in 1:n) {
          d = 0
          eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to d
          Vi = rTensor::k_unfold(rTensor::as.tensor(d), m = k)@data  # unfold tensor
          S.array[, , i] = Vi %*% t(Vi)
        }
        S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # sample estimation of \Sigma_k
        Omega_tilde = solve(S.mat)
        # normalization
        if (normalize) {
          Omega_tilde = Omega_tilde / norm(Omega_tilde, type = "F")}
      }
      Omega_tilde
    }
    stopCluster(c1)
  }
  
  K1 = length(est.mode) # number of precision matrices to be estimated
  lam.best = rep(0, K1)
  loglik = list() 
  
  ###### Tuning process ######
  # When validation set is supplied, the lambdas with the maximum log-likelihood will be chosen
  if (!(is.null(val))) {
    Omega.list = list() # list of \tilde\Omega
    Omega.list.sqrt = list() # list of square root of \tilde\Omega
    for (k in 1:K) {
      Omega.list[[k]] = fit1[[k]]
      Omega.list.sqrt[[k]] = sqrtm(Omega.list[[k]])
    }
    Omega.sqrt.copy = Omega.list.sqrt
    
    for (mode_index in 1:K1) {
      k = est.mode[mode_index]
      
      # Calculate \tilde S_k using the training set
      S.array = array(0, c(m.vec[k], m.vec[k], n))
      Omega.list.sqrt[[k]] = diag(m.vec[k]) # set \tilde\Omega_k to identity matrix
      for (i in 1:n) {
        d = 0
        eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]"))) # assign the ith observation to d
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        S.array[, , i] = Vi %*% t(Vi)
      }
      S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      
      # Calculate \tilde S_k using the validation set
      testS.array = array(0, c(m.vec[k], m.vec[k], n_val))
      for (i in 1:n_val) {
        d = 0
        eval(parse(text = paste("d=val[", paste(rep(",", K), collapse = ""), "i]")))
        Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                    ms = 1:K
        )@data), m = k)@data
        testS.array[, , i] = Vi %*% t(Vi)
      }
      testS.mat = apply(testS.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
      Omega.list.sqrt[[k]] = Omega.sqrt.copy[[k]]
      
      # fit model with a sequence of lambdas
      lamk = lambda.list[[mode_index]] # a sequence of candidates for lambda_k
      lam.length = length(lamk)
      loglik2 = rep(0, lam.length)
      for (i in 1:lam.length) {
        Out1 = glasso(S.mat, rho = lamk[i], penalize.diagonal = FALSE, maxit = 1e4, thr = 1.0e-4)
        hat_Omega = Out1$wi
        loglik2[i] = -tr(testS.mat %*% hat_Omega) + log(det(hat_Omega * scale.vec[mode_index]))
        if (loglik2[i] == Inf) {
          stop(paste("Infinite value! Please choose a smaller scale for mode", mode_index))
        }
        if (loglik2[i] == -Inf) {
          stop(paste("Negative infinite value! Please choose a larger scale for mode", mode_index))
        }
      }
      ind = which.max(loglik2)
      lam.best[mode_index] = lamk[ind] # get the optimal lambda that maximizes the log-likelihood
      loglik[[mode_index]] = loglik2
    }
  } else {
    # if validation set is not provided, directly use lambda.vec to fit model
    lam.best = lambda.vec
  }
  
  
  ##### Model fitting using parallel computing #####
  # register cluster for parallel computing
  c1 = makeCluster(njobs)
  registerDoParallel(c1)
  K1 = length(est.mode)
  fit_result = foreach(mode_ind = 1:K1, .packages = c("glasso", "rTensor", "expm"), .export = c("x"), .combine = list, .multicombine = TRUE) %dopar% {
    k = est.mode[mode_ind]
    Omega.list.sqrt = list()
    for (i in 1:K) {
      Omega.list.sqrt[[i]] = sqrtm(fit1[[i]])
    }
    # Calculate \tilde S_k
    S.array = array(0, c(m.vec[k], m.vec[k], n))
    Omega.list.sqrt[[k]] = diag(m.vec[k])
    for (i in 1:n) {
      d = 0
      eval(parse(text = paste("d=x[", paste(rep(",", K), collapse = ""), "i]")))
      Vi = k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                  ms = 1:K
      )@data), m = k)@data
      S.array[, , i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array, c(1, 2), mean) * m.vec[k] / prod(m.vec) # \tilde S_k
    # fit model
    Out1 = glasso(S.mat, rho = lam.best[mode_ind], penalize.diagonal = FALSE, maxit = maxit, thr = thres)
    hat_Omega = as.matrix(Out1$wi)
    normalization
    if (normalize) {
      hat_Omega = (hat_Omega / norm(hat_Omega, type = "F"))
    }
    list(hat_Omega, S.mat)
  }
  stopCluster(c1)
  
  result = list()
  result$fit_result = fit_result
  result$lambda = lam.best
  result$loglik = loglik
  return(result)
}


# Model setting
n <- 20 # sample size
dimen <- c(30, 36, 30) # dimension of tensor
nvars <- prod(dimen) # number of variables
K <- 3 # order of tensor
Run <- 100
Run <- 2


d <- 1
av.error.f <- array(0, dim = c(Run, d)) # averaged estimation error in Frobenius norm
av.error.max <- array(0, dim = c(Run, d)) # averaged estimation error in Maximum norm
av.tpr <- array(0, dim = c(Run, d)) # averaged true positive rate
av.tnr <- array(0, dim = c(Run, d)) # averaged true negative rate

d <- 3
error.f <- array(0, dim = c(Run, d)) # estimation error in Frobenius norm for each mode
error.max <- array(0, dim = c(Run, d)) # estimation error in Maximum norm for each mode
tpr <- array(0, dim = c(Run, d)) # true positive rate for each mode
tnr <- array(0, dim = c(Run, d)) # true negative rate for each mode

d <- 1
av.error.f.T <- array(0, dim = c(Run, d)) # averaged estimation error in Frobenius norm
av.error.max.T <- array(0, dim = c(Run, d)) # averaged estimation error in Maximum norm
av.tpr.T <- array(0, dim = c(Run, d)) # averaged true positive rate
av.tnr.T <- array(0, dim = c(Run, d)) # averaged true negative rate

d <- 3
error.f.T <- array(0, dim = c(Run, d)) # estimation error in Frobenius norm for each mode
error.max.T <- array(0, dim = c(Run, d)) # estimation error in Maximum norm for each mode
tpr.T <- array(0, dim = c(Run, d)) # true positive rate for each mode
tnr.T <- array(0, dim = c(Run, d)) # true negative rate for each mode


for (run in 1:Run) { 
  print(run)
  # Generate training set and validation set
  data <- Model7(n, run * 123456)
  x <- data[[1]]$x
  vax <- data[[1]]$vax
  Sigma <- data[[2]]
  Omega <- data[[3]]
  scales <- data[[4]]
  
  Tx <- NSEstimator(x, dimen)
  Tvax <- NSEstimator(vax, dimen)
  
  # proper candidates of tuning parameters
  # lamseq <- seq(1.5e-09, 0.2, length.out = 100)
  lamseq <- seq(1e-09, 1e-04, length.out = 2000)
  lambda.list <- list() # a list containing candidates of tuning parameters for each mode 
  for (i in 1:K) {
    lambda.list[[i]] <- lamseq
  }
  
  # Model fitting
  fit <- Separate.fit(x, vax, lambda.list = lambda.list)
  # lamseq <- seq(1.5e-03, 0.2, length.out = 100)
  # # lamseq <- seq(1.5e-09, 2e-05, length.out = 50)
  # lambda.list <- list() # a list containing candidates of tuning parameters for each mode 
  # for (i in 1:K) {
  #   lambda.list[[i]] <- lamseq
  # }
  # fit.T <- Separate.fit(Tx, lambda.vec = fit$lambda)
  fit.T <- Separate.fit(Tx, Tvax, lambda.list = lambda.list)
  
  ## If there is no validation set, we can use cv.Separate to tune lambda via cross-validation
  # fit <- cv.Separate(x, lambda.list=lambda.list)
  ## With a user-specified sequence of lambdas in lambda.vec, we can directly fit the model
  # fit <- Separate.fit(x, lambda.vec = lambda.vec)
  
  # Simulation summary of estimation errors, TPR and TNR
  out <- simulation.summary(purrr::map2(fit$fit_result, scales, \(x, y) x[[1]]*y), Omega, offdiag = FALSE)
  # out <- simulation.summary(purrr::map(fit$fit_result, \(x) x[[1]]), purrr::map(fit.T$fit_result, \(x) x[[1]]), offdiag = FALSE)
  # out <- simulation.summary(purrr::map(fit$fit_result, \(x) x[[2]]), purrr::map(fit.T$fit_result, \(x) x[[2]]), offdiag = FALSE)
  av.error.f[run] <- out$av.error.f
  av.error.max[run] <- out$av.error.max
  av.tpr[run] <- out$av.tpr
  av.tnr[run] <- out$av.tnr
  
  error.f[run, ] <- out$error.f
  error.max[run, ] <- out$error.max
  tpr[run, ] <- out$tpr
  tnr[run, ] <- out$tnr
  
  # out2 <- simulation.summary(purrr::map(fit.T$fit_result, \(x) x[[2]]), Sigma, offdiag = FALSE)
  out <- simulation.summary(purrr::map2(fit.T$fit_result, scales, \(x, y) x[[1]]*y), Omega, offdiag = FALSE)
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

########################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

Sigma_fit <- purrr::map(fit$fit_result, \(x) x[[2]])
Sigma_fit.T <- purrr::map(fit.T$fit_result, \(x) x[[2]])
simulation.summary(Sigma_fit.T, Sigma_fit, offdiag = FALSE)
simulation.summary(purrr::map(fit$fit_result, \(x) x[[1]]), purrr::map(fit.T$fit_result, \(x) x[[1]]), offdiag = FALSE)
simulation.summary(purrr::map(fit$fit_result, \(x) x[[1]]), Omega, offdiag = FALSE)
simulation.summary(purrr::map(fit.T$fit_result, \(x) x[[1]]), Omega, offdiag = FALSE)




# 1. Plot the first vector
plot(c(Sigma[[1]]), c(Sigma_fit[[1]]), 
     type = "l", 
     col = "blue", 
     lwd = 2,
     xlab = "Index / Sigma", 
     ylab = "Values", 
     asp=1,
     # xlim = c(0, 1),
     # ylim = c(0, 1),
     main = "Comparison of Estimated Vectors")
abline(a = 0, b = 1, col = "red", lwd = 2)

# 2. Add the second vector as a line
lines(c(Sigma_fit[[1]]), c(est_Sigmax[[1]]), 
      col = "red", 
      lwd = 2)

lines(c(Sigma_fit[[1]]), c(Sigma[[1]]), 
      col = "black", 
      lwd = 2)

# 3. Add a legend to distinguish the lines
legend("bottomright", 
       legend = c("est_SigmaTx", "est_Sigmax"), 
       col = c("blue", "red",), 
       lty = 1, 
       lwd = 2)


