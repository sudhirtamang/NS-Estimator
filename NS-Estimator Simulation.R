# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   NS-Estimator
# PURPOSE:  Compare the performance of NS-Estimator in estimating precision matrix
#           in the presence of contamination
# AUTHOR:   Sudhir T
# DATE:     5th April, 2026
# ==============================================================================
rm(list = ls())
library(Tlasso)
library(tensr)
library(expm)
library(rTensor)
library(glasso)
library(mvtnorm)
library(parallel)
library(doParallel)
source("simulation.summary.R")
source("Lfunctions.R")
source("C://Users//sudhi//Desktop//Coordinatewise Gaussianization//replication//List_functions_for_cellwise_corruption.R")

dimen <- c(30, 36, 30)
nvars <- prod(dimen) 
K <- 3
seed <- 1
n <- 200
pctOut <- 0.1
normalize <- TRUE
maxit <- 1e4
thres <- 1.0e-4

Omega <- purrr::map2(dimen, 100 * 1:3, ChainOmega, norm.type = 2) # precision matrix for each mode
# Omega <- purrr::map2(dimen, dimen, generate_sparse_omega) # precision matrix for each mode
# Omega <- purrr::map(Omega, \(x) x/norm(x, type="F"))
Sigma <- purrr::map(Omega, solve) 
dSigma <- purrr::map(purrr::map(Sigma, chol), t)
Omega.list.sqrt <- purrr::map(Omega, sqrtm)
x <- array(0, dim = c(dimen, n))
vec_x <- matrix(rnorm(nvars * n), ncol = n) 
for(i in 1:n) {
  x[, , , i] <- atrans(array(vec_x[, i], dimen), dSigma)
}

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


# lambdas <- exp(seq(log(15e-09), log(0.004), length.out=50))
# lambdas <- seq(0.001, 0.002, length.out=100)
# av.tnr <- vector("double", length(lambdas))
# av.tpr <- vector("double", length(lambdas))
# fit1 = omega_from_sample(TcontamiData, dimen, n)
# for(i in seq_along(lambdas)){
#   fit_result <- fast_sep(TcontamiData, dimen, n, fit1, rep(lambdas[[i]], 3), normalize=TRUE, thres=1.0e-4, maxit=1e4)
#   
#   res <- simulation.summary(fit_result, Omega, offdiag = FALSE)
#   av.tnr[[i]] <- res$av.tnr
#   av.tpr[[i]] <- res$av.tpr
# }
# 
# plot(lambdas, av.tnr,type="l", col="red", ylab="")
# lines(lambdas, av.tpr, col="black")
# legend("bottomright", legend = c("av.tnr", "av.tpr"), col = c("red", "black"), lty = 1)





# ===================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ===================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ===================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ====================================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# k <- 2
# datax <- x
# datavax <- vax
# lambdas <- seq(0.000015, 0.06, length.out = 50)
# 
# k <- 2
datax <- TcontamiData
datavax <- vax_TcontamiData
# lambdas <- seq(1.5e-09, 2e-4, length.out = 50)
# 
# fit1 = omega_from_sample(datax, dimen, k, n)
# fit2 = omega_from_sample(datavax, dimen, k, n)
# # Mat_Sk_x <- tilde_Sk(datax, fit1, dimen, k, n)
# # Mat_Sk_vax <- tilde_Sk(datavax, fit1, dimen, k, n)
# Mat_Sk_x <- tilde_Sk(datax, Omega, dimen, k, n)
# Mat_Sk_vax <- tilde_Sk(datavax, Omega, dimen, k, n)


lam.best <- vector("double", length(dimen))
for(j in seq_along(lam.best)){
  fit1 = omega_from_sample(datax, dimen, j, n)
  fit2 = omega_from_sample(datavax, dimen, j, n)
  Mat_Sk_x <- tilde_Sk(datax, fit1, dimen, j, n)
  Mat_Sk_vax <- tilde_Sk(datavax, fit1, dimen, j, n)
  lambdas <- seq(0.0000000015, 2e-5, length.out = 50)
  Loglikelihood <- vector("double", length(lambdas))
  for(i in seq_along(lambdas)){
    out <- glasso::glasso(Mat_Sk_x, rho = lambdas[[i]], penalize.diagonal = FALSE, maxit = maxit, thr = thres, n=n)
    # Loglikelihood[[i]] <- -tr(fit2[[1]] %*% Omega[[1]]) + log(det(Omega[[1]]))
    Loglikelihood[[i]] <- -tr(Mat_Sk_vax %*% out$wi) + log(det(out$wi))
  }
  lam.best[[j]] <- lambdas[[which.max(Loglikelihood)]]
  plot(lambdas, Loglikelihood)
}
# plot(lambdas, Loglikelihood)
# print(Loglikelihood)
# ====================================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



fit_result <- fast_sep(TcontamiData, dimen, n, fit1, lam.best, normalize=TRUE, thres=1.0e-6, maxit=1e4)
# fit_result <- fast_sep(TcontamiData, dimen, n, fit1, rep(1.6e-05, 3), normalize=TRUE, thres=1.0e-4, maxit=1e4)

simulation.summary(fit_result, Omega, offdiag = FALSE)






