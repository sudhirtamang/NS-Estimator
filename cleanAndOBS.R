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
source("Model.R")


RUNs <- 1
n <- 50
dimen <- c(45, 54)



data <- Model(50, 123456, dimen)
x <- data[[1]]$x
vax <- data[[1]]$vax
Sigma <- data[[2]]
Omega <- data[[3]]


pct <- 0.5
Sigma_obs <- (1 - pct)*Sigma[[1]] + pct*diag(dimen[[1]])
Omega_obs <- solve(Sigma_obs)


# PCTs <- seq(0, 1, 0.1)
# NORMdiff <- vector("double", length=length(contaminations))
# NORMmax <- vector("double", length=length(contaminations))
# for(i in seq_along(contaminations)){
#   Sigma_obs <- (1 - PCTs[[i]])*Sigma[[1]] + PCTs[[i]]*diag(dimen[[1]])
#   NORMdiff[[i]] <- norm(Sigma_obs - Sigma[[1]], type="F")
#   NORMmax[[i]] <- norm(Sigma_obs - Sigma[[1]], type="M")
# }
# 
# plot(NORMdiff)
# points(NORMmax)



cat("for 50% Contamination", "\n")

cat("Frob. norm of Clean Sigma:", norm(Sigma[[1]], type="F"), "\n")
cat("Norm Max. of Clean Sigma:", norm(Sigma[[1]], type="M"), "\n")
cat("Frob. norm of Contaminated Sigma:", norm(Sigma_obs, type="F"), "\n")
cat("Norm Max. of Contaminated Sigma:", norm(Sigma_obs, type="M"), "\n")

cat("\n\n")

cat("Frob. norm of Clean Precision Matrix:", norm(Omega[[1]], type="F"), "\n")
cat("Norm Max. of Clean Precision Matrix:", norm(Omega[[1]], type="M"), "\n")
cat("Frob. norm of Contaminated Precision Matrix:", norm(Omega_obs, type="F"), "\n")
cat("Norm Max. of Contaminated Precision Matrix:", norm(Omega_obs, type="M"), "\n")

cat("\n\n")

cat("Frob. norm diff between Clean Sigma and Contaminated Sigma:", norm(Sigma[[1]] - Sigma_obs, type="F"), "\n")
cat("Norn Max. of diff between Clean Sigma and Contaminated Sigma:", norm(Sigma[[1]] - Sigma_obs, type="M"), "\n")
cat("Frob. norm diff between Clean Precision Matrix and Contaminated Precision Matrix:", norm(Omega[[1]] - Omega_obs, type="F"), "\n")
cat("Norn Max. of diff between Clean Precision Matrix and Contaminated Precision Matrix:", norm(Omega[[1]] - Omega_obs, type="M"), "\n")





