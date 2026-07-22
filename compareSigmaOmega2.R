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



source("simulation.summary.R")
source("Lfunctions.R")
source("Separate.fit.correct.R")
source("Separate.fit.correct.plot.R")
source("Model.R")





pctOut <- 0.35
RUNs <- 100
# RUNs <- 2
n <- 50
dimen <- c(45, 2)
# dimen <- c(30, 36, 30)
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
contam.lambda.x <- array(0, dim = c(RUNs, d))
contam.est.sigma <- array(0, dim = c(RUNs, d))
contam.est.omega <- array(0, dim = c(RUNs, d))
contam.error.f <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contam.error.f.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contam.error.max <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contam.error.max.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contam.tpr <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
contam.tnr <- array(0, dim = c(RUNs, d)) # true negative rate for each mode


d <- 1
# d <- K
contam.lambda.Tx <- array(0, dim = c(RUNs, d))
contam.est.sigma.Tx <- array(0, dim = c(RUNs, d))
contam.est.omega.Tx <- array(0, dim = c(RUNs, d))
contam.error.f.Tx <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
contam.error.f.Tx.sigma <- array(0, dim = c(RUNs, d)) # averaged estimation error in Frobenius norm
contam.error.max.Tx <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
contam.error.max.Tx.sigma <- array(0, dim = c(RUNs, d)) # averaged estimation error in Maximum norm
contam.tpr.Tx <- array(0, dim = c(RUNs, d)) # averaged true positive rate
contam.tnr.Tx <- array(0, dim = c(RUNs, d)) # averaged true negative rate

d <- 1
# d <- K
contam.lambda.Tx.C <- array(0, dim = c(RUNs, d))
contam.est.sigma.Tx.C <- array(0, dim = c(RUNs, d))
contam.est.omega.Tx.C <- array(0, dim = c(RUNs, d))
contam.error.f.Tx.C <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contam.error.f.Tx.C.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Frobenius norm for each mode
contam.error.max.Tx.C <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contam.error.max.Tx.C.sigma <- array(0, dim = c(RUNs, d)) # estimation error in Maximum norm for each mode
contam.tpr.Tx.C <- array(0, dim = c(RUNs, d)) # true positive rate for each mode
contam.tnr.Tx.C <- array(0, dim = c(RUNs, d)) # true negative rate for each mode

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
lamseq <- seq(1.5e-3, 1.5e-1, length.out = 300)
lambda.list <- list() # a list containing candidates of tuning parameters for each mode
for (i in 1:K) {
  lambda.list[[i]] <- lamseq
}

lamseq.C <- seq(1.5e-6, 1, length.out = 300)
lambda.list.C <- list() # a list containing candidates of tuning parameters for each mode
for (i in 1:K) {
  lambda.list.C[[i]] <- lamseq.C
}


set.seed(1234)
all.values <- sample(RUNs, 8)
x.values <- sort(all.values[1:2])
contam.x.values <- sort(all.values[3:4])
contam.Tx.values <- sort(all.values[5:6])
contam.Tx.values.C <- sort(all.values[7:8])
xlab <- "Tuning parameter"
ylab <- "Log-likelihood"

for(itr in 1:RUNs){
  
  print(itr)
  data <- Model(n, itr * 123456, dimen)
  x <- data[[1]]$x
  vax <- data[[1]]$vax
  Sigma <- data[[2]]
  Omega <- data[[3]]
  
  nOut <- ceiling(pctOut * n)
  idxcontami <- sample(n, nOut)
  contami <- rt(nOut*nvars, df=6)
  dim(contami) <- c(dimen, nOut)
  
  contam.x <- x[]
  for(i in seq_along(idxcontami)){
    # contam.x[, , , idxcontami[i]] <- contami[, , , i]
    contam.x[, , idxcontami[i]] <- contami[, , i]
  }
  
  contam.Tx <- NSEstimator2(contam.x, dimen)
  # contam.Tx <- NSEstimator(contam.x, dimen)
  
  
  idxcontami <- sample(n, nOut)
  contami <- rt(nOut*nvars, df=6)
  dim(contami) <- c(dimen, nOut)
  
  contam.vax <- vax[]
  for(i in seq_along(idxcontami)){
    # contam.vax[, , , idxcontami[i]] <- contami[, , , i]
    contam.vax[, , idxcontami[i]] <- contami[, , i]
  }
  
  contam.Tvax <- NSEstimator2(contam.vax, dimen)
  # contam.Tvax <- NSEstimator(contam.vax, dimen)
  
  # Tx <- NSEstimator(x, dimen)
  # Tvax <- NSEstimator(vax, dimen)
  Tx <- NSEstimator2(x, dimen)
  Tvax <- NSEstimator2(vax, dimen)
  if(itr %in% x.values){
    main <- paste0("Clean data, no transformation, no correction at ", itr, " Iteration")
    fitx <- Separate.fit.correct.plot(x, vax, lambda.list = lambda.list, main=main, xlab=xlab, ylab=ylab)
  } else {
    fitx <- Separate.fit.correct(x, vax, lambda.list = lambda.list)
  }
  if(itr %in% contam.x.values){
    main <- paste0("Contaminated data, no transformation, no correction at ", itr, " Iteration")
    contam.fitx <- Separate.fit.correct.plot(contam.x, contam.vax, lambda.list = lambda.list, main=main, xlab=xlab, ylab=ylab)
  } else {
    contam.fitx <- Separate.fit.correct(contam.x, contam.vax, lambda.list = lambda.list.C)
  }
  if(itr %in% contam.Tx.values){
    main <- paste0("Contaminated data, with transformation, no correction at ", itr, " Iteration")
    contam.fitTx <- Separate.fit.correct.plot(contam.Tx, contam.Tvax, lambda.list = lambda.list, main=main, xlab=xlab, ylab=ylab)
  } else {
    contam.fitTx <- Separate.fit.correct(contam.Tx, contam.Tvax, lambda.list = lambda.list.C)
  }
  if(itr %in% contam.Tx.values.C){
    main <- paste0("Contaminated data, with transformation, with correction at ", itr, " Iteration")
    contam.fitTx.C <- Separate.fit.correct.plot(contam.Tx, contam.Tvax, lambda.list = lambda.list, main=main, xlab=xlab, ylab=ylab)
  } else {
    contam.fitTx.C <- Separate.fit.correct(contam.Tx, contam.Tvax, lambda.list = lambda.list.C, Grho=Grho)
  }
  
  
  
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
  
  
  
  contam.outx <- simulation.summary(list(contam.fitx$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  contam.outx.sigma <- simulation.summary(list(contam.fitx$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  contam.error.f[itr, 1] <- contam.outx$error.f
  contam.error.max[itr, 1] <- contam.outx$error.max
  contam.tpr[itr, 1] <- contam.outx$tpr
  contam.tnr[itr, 1] <- contam.outx$tnr
  contam.est.sigma[itr, 1] <- norm(contam.fitx$Omegahat[[1]][[2]], type="F")
  contam.est.omega[itr, 1] <- norm(contam.fitx$Omegahat[[1]][[1]], type="F")
  
  contam.error.f.sigma[itr, 1] <- contam.outx.sigma$error.f
  contam.error.max.sigma[itr, 1] <- contam.outx.sigma$error.max
  
  contam.lambda.x[itr, 1] <- contam.fitx$lambda[[1]]
  
  
  
  
  
  contam.outTx <- simulation.summary(list(contam.fitTx$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  contam.outTx.sigma <- simulation.summary(list(contam.fitTx$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  contam.error.f.Tx[itr, 1] <- contam.outTx$error.f
  contam.error.max.Tx[itr, 1] <- contam.outTx$error.max
  contam.tpr.Tx[itr, 1] <- contam.outTx$tpr
  contam.tnr.Tx[itr, 1] <- contam.outTx$tnr
  contam.est.sigma.Tx[itr, 1] <- norm(contam.fitTx$Omegahat[[1]][[2]], type="F")
  contam.est.omega.Tx[itr, 1] <- norm(contam.fitTx$Omegahat[[1]][[1]], type="F")
  
  contam.error.f.Tx.sigma[itr, 1] <- contam.outTx.sigma$error.f
  contam.error.max.Tx.sigma[itr, 1] <- contam.outTx.sigma$error.max
  
  contam.lambda.Tx[itr, 1] <- contam.fitTx$lambda[[1]]
  
  
  
  
  contam.outTx.C <- simulation.summary(list(contam.fitTx.C$Omegahat[[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)
  contam.outTx.C.sigma <- simulation.summary(list(contam.fitTx.C$Omegahat[[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
  contam.error.f.Tx.C[itr, 1] <- contam.outTx.C$error.f
  contam.error.max.Tx.C[itr, 1] <- contam.outTx.C$error.max
  contam.tpr.Tx.C[itr, 1] <- contam.outTx.C$tpr
  contam.tnr.Tx.C[itr, 1] <- contam.outTx.C$tnr
  contam.est.sigma.Tx.C[itr, 1] <- norm(contam.fitTx.C$Omegahat[[1]][[2]], type="F")
  contam.est.omega.Tx.C[itr, 1] <- norm(contam.fitTx.C$Omegahat[[1]][[1]], type="F")
  
  contam.error.f.Tx.C.sigma[itr, 1] <- contam.outTx.C.sigma$error.f
  contam.error.max.Tx.C.sigma[itr, 1] <- contam.outTx.C.sigma$error.max
  
  contam.lambda.Tx.C[itr, 1] <- contam.fitTx.C$lambda[[1]]
}

cat("Comparison for OMEGA::clean data, no transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(error.f), "SD:", sd(error.f), "\n")
cat("Mean Max. Error:", colMeans(error.max), "SD:", sd(error.max), "\n")
cat("TPR:", colMeans(tpr), "SD:", sd(tpr), "\n")
cat("TNR:", colMeans(tnr), "SD:", sd(tnr), "\n")
cat("Comparison for SIGMA::just sample, no transformation, no correction", "\n")
cat("Range of lambda.best: ", sprintf("[%.10e, %.10e]\n", min(lambda.x), max(lambda.x)), "\n\n")


cat("Comparison for SIGMA::clean data, no transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(error.f.sigma), "SD:", sd(error.f.sigma), "\n")
cat("Mean Max. Error:", colMeans(error.max.sigma), "SD:", sd(error.max.sigma), "\n")
cat("Mean Frob. Norm Estimated Sigma:", colMeans(est.sigma), "SD:", sd(est.sigma), "\n")
cat("Frob. Norm True Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Mean Frob. Norm Estimated Omega:", colMeans(est.omega), "SD:", sd(est.omega), "\n")
cat("Frob. Norm True Omega:", norm(Omega[[1]], type="F"), "\n")



cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("Comparison for OMEGA:: contaminated data, no transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(contam.error.f), "SD:", sd(contam.error.f), "\n")
cat("Mean Max. Error:", colMeans(contam.error.max), "SD:", sd(contam.error.max), "\n")
cat("TPR:", colMeans(contam.tpr), "SD:", sd(contam.tpr), "\n")
cat("TNR:", colMeans(contam.tnr), "SD:", sd(contam.tnr), "\n")
cat("Comparison for SIGMA::just sample, no transformation, no correction", "\n")
cat("Range of lambda.best: ", sprintf("[%.10e, %.10e]\n", min(contam.lambda.x), max(contam.lambda.x)), "\n\n")


cat("Comparison for SIGMA:: contaminated data, no transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(contam.error.f.sigma), "SD:", sd(contam.error.f.sigma), "\n")
cat("Mean Max. Error:", colMeans(contam.error.max.sigma), "SD:", sd(contam.error.max.sigma), "\n")
cat("Mean Frob. Norm Estimated Sigma:", colMeans(contam.est.sigma), "SD:", sd(contam.est.sigma), "\n")
cat("Frob. Norm True Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Mean Frob. Norm Estimated Omega:", colMeans(contam.est.omega), "SD:", sd(contam.est.omega), "\n")
cat("Frob. Norm True Omega:", norm(Omega[[1]], type="F"), "\n")



cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")


cat("Comparison for OMEGA::Contaminated data with transformation, no correction", "\n")

cat("Mean Forb. Diff.:", colMeans(contam.error.f.Tx), "SD:", sd(contam.error.f.Tx), "\n")
cat("Mean Max. Error:", colMeans(contam.error.max.Tx), "SD:", sd(contam.error.max.Tx), "\n")
cat("TPR:", colMeans(contam.tpr.Tx), "SD:", sd(contam.tpr.Tx), "\n")
cat("TNR:", colMeans(contam.tnr.Tx), "SD:", sd(contam.tnr.Tx), "\n")
cat("Range of lambda.best: ", sprintf("[%.10e, %.10e]\n", min(contam.lambda.Tx), max(contam.lambda.Tx)), "\n\n")

cat("Comparison for SIGMA::Contaminated data with transformation, no correction", "\n")
cat("Mean Forb. Diff.:", colMeans(contam.error.f.Tx.sigma), "SD:", sd(contam.error.f.Tx.sigma), "\n")
cat("Mean Max. Error:", colMeans(contam.error.max.Tx.sigma), "SD:", sd(contam.error.max.Tx.sigma), "\n")
cat("Mean Frob. Norm Estimated Sigma:", colMeans(contam.est.sigma.Tx), "SD:", sd(contam.est.sigma.Tx), "\n")
cat("Frob. Norm True Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Mean Frob. Norm Estimated Omega:", colMeans(contam.est.omega.Tx), "SD:", sd(contam.est.omega.Tx), "\n")
cat("Frob. Norm True Omega:", norm(Omega[[1]], type="F"), "\n")



cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")
cat("====================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "\n")


cat("Comparison for OMEGA::Contaminated data with transformation with correction", "\n")
cat("Mean Forb. Diff.:", colMeans(contam.error.f.Tx.C), "SD:", sd(contam.error.f.Tx.C), "\n")
cat("Mean Max. Error:", colMeans(contam.error.max.Tx.C), "SD:", sd(contam.error.max.Tx.C), "\n")
cat("TPR:", colMeans(contam.tpr.Tx.C), "SD:", sd(contam.tpr.Tx.C), "\n")
cat("TNR:", colMeans(contam.tnr.Tx.C), "SD:", sd(contam.tnr.Tx.C), "\n")
cat("Range of lambda.best: ", sprintf("[%.10e, %.10e]\n", min(contam.lambda.Tx.C), max(contam.lambda.Tx.C)), "\n\n")

cat("Comparison for SIGMA::Contaminated data with transformation with correction", "\n")
cat("Mean Forb. Diff.:", colMeans(contam.error.f.Tx.C.sigma), "SD:", sd(contam.error.f.Tx.C.sigma), "\n")
cat("Mean Max. Error:", colMeans(contam.error.max.Tx.C.sigma), "SD:", sd(contam.error.max.Tx.C.sigma), "\n")
cat("Mean Frob. Norm Estimated Sigma:", colMeans(contam.est.sigma.Tx.C), "SD:", sd(contam.est.sigma.Tx.C), "\n")
cat("Frob. Norm True Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Mean Frob. Norm Estimated Omega:", colMeans(contam.est.omega.Tx.C), "SD:", sd(contam.est.omega.Tx.C), "\n")
cat("Frob. Norm True Omega:", norm(Omega[[1]], type="F"), "\n")




