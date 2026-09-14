rm(list = ls())
library(Tlasso)
library(tensr)
library(glasso)
library(expm)
library(rTensor)
library(doParallel)
library(furrr)
library(purrr)
library(mvnfast)
library(tictoc)


source("C:/Users/sudhi/Desktop/Coordinatewise Gaussianization/replication/List_functions_for_cellwise_corruption.R")
source("C:/Users/sudhi/Desktop/PhD projects/NS-Estimator/Lfunctions.R")



sampleCOVest <- function(X, hat_sigma) {
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  return(crossprod(X_centered) / (nrow(X) - 1))
}



Model <- function(n, seed, p, rmvnorm_precision_n) {
  
  Omega <- ChainOmega(p, sd = 100, norm.type = 2)
  Sigma <- solve(Omega)
  
  set.seed(seed)
  
  # Generate data observation
  # training set
  x <- rmvnorm_precision_n(n, p, Omega)
  
  # validation set
  vax <- rmvnorm_precision_n(n, p, Omega)
  
  result <- list()
  result$x <- x
  result$vax <- vax
  
  return(list(result, Sigma, Omega))
}


list.packages <- c("mvnfast", "Tlasso", "purrr",  "glasso", "QUIC", "caret", "R.utils", "tictoc")
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



doSimulation2 <- function(Grho, is.corrected, isNS.transform, lower, p, n, R, p_conta, estSIGMA_hat, generate_Omega, medDev,
                          rmvnorm_precision_n, TPrate, TNrate, FNrate, FPrate, text1, text2){
  total_cores <- parallel::detectCores()
  c1 <- parallel::makeCluster(floor(total_cores * 0.7))
  
  doParallel::registerDoParallel(c1)
  Results <- foreach(k = seq_len(R), .packages = list.packages,	.combine = list, .multicombine = TRUE) %dopar% {
    Omega <- generate_Omega(p, p)
    Solved_Omega <- solve(Omega)
    data <- rmvnorm_precision_n(n, p, Omega)
    
    
    #   Add contamination to the data
    for(i in 1:n){
      B <- diag(rbinom(p, 1, prob=p_conta))
      data[i, ] <- (diag(p) - B) %*% data[i, ] + B %*% rnorm(p, 10, sqrt(0.2))
    }
    
    if(isNS.transform){
      for(i in 1:p){
        data[, i] <- qnorm(rank(data[, i])/(n + 1))
      }
    }
    
    q75 <- qnorm(0.75)
    hat_sigma <- apply(data, 2, medDev)/q75
    
    # tic()
    # Covariance Mat est. ----------------------------------------->>>>>>>>>>>>>>
    SIGMA_hat <- estSIGMA_hat(data, hat_sigma)
    
    if(is.corrected){
      RHOs <- seq(-1, 1, 0.01)
      for(i in 1:p){
        for(j in 1:p){
          tmp1 <- abs(SIGMA_hat[i, j] - Grho)
          SIGMA_hat[i, j] <- RHOs[[which.min(tmp1)]]
        }
      }
      SIGMA_hat <- as.matrix(Matrix::nearPD(SIGMA_hat)$mat)
    }
    # Cross-validation to choose optimal penalizing parameter
    lambda_max <- max(abs(SIGMA_hat[row(SIGMA_hat) != col(SIGMA_hat)]))
    lambda_min <- lower * lambda_max
    LAMBDAs <- exp(seq(log(lambda_min), log(lambda_max), length.out = 15))
    
    fold_indices <- caret::createFolds(1:n, k = 5, list = TRUE)
    # --- Precompute fold covariance estimates ---
    fold_SIGMA_train <- vector("list", 5)
    fold_SIGMA_test  <- vector("list", 5)
    
    for (j in seq_along(fold_indices)) {
      test  <- data[fold_indices[[j]], ]
      train <- data[-fold_indices[[j]], ]
      
      hat_sigma_test  <- apply(test, 2, medDev) / q75
      hat_sigma_train <- apply(train, 2, medDev) / q75
      
      fold_SIGMA_test[[j]]  <- estSIGMA_hat(test, hat_sigma_test)
      if(is.corrected){
        RHOs <- seq(-1, 1, 0.01)
        for(i in 1:p){
          for(k in 1:p){
            tmp1 <- abs(fold_SIGMA_test[[j]][i, j] - Grho)
            fold_SIGMA_test[[j]][i, k] <- RHOs[[which.min(tmp1)]]
          }
        }
        fold_SIGMA_test[[j]] <- as.matrix(Matrix::nearPD(fold_SIGMA_test[[j]])$mat)
      }
      fold_SIGMA_train[[j]] <- estSIGMA_hat(train, hat_sigma_train)
      if(is.corrected){
        RHOs <- seq(-1, 1, 0.01)
        for(i in 1:p){
          for(k in 1:p){
            tmp1 <- abs(fold_SIGMA_train[[j]][i, j] - Grho)
            fold_SIGMA_train[[j]][i, k] <- RHOs[[which.min(tmp1)]]
          }
        }
        fold_SIGMA_train[[j]] <- as.matrix(Matrix::nearPD(fold_SIGMA_train[[j]])$mat)
      }
    }
    # --- Lambda loop ---
    NegLoglikelihood <- numeric(length(LAMBDAs))
    for (i in seq_along(LAMBDAs)) {
      tmp0 <- numeric(5)
      rho_matrix <- matrix(LAMBDAs[[i]], nrow = p, ncol = p)
      diag(rho_matrix) <- 0
      for (j in seq_along(fold_indices)) {
        fit <- QUIC(fold_SIGMA_train[[j]], rho_matrix,
                    tol = 1e-4, maxIter = 1000, msg = 0, path = NULL)
        Theta <- matrix(fit$X, p, p)
        R1 <- chol(Theta)
        log_det <- 2 * sum(log(diag(R1)))
        
        trace   <- sum(fold_SIGMA_test[[j]] * Theta)
        tmp0[j] <- -log_det + trace
      }
      NegLoglikelihood[[i]] <- mean(tmp0)
    }
    rho_matrixopt <- matrix(LAMBDAs[which.min(NegLoglikelihood)], nrow = p, ncol = p)
    diag(rho_matrixopt) <- 0
    fit <- QUIC(SIGMA_hat, rho_matrixopt,
                tol = 1e-4,
                maxIter = 1000,
                msg = 2, path=NULL)
    OMEGA_hat <- matrix(fit$X, p, p)
    # OMEGA_hat <- OMEGA_hat/norm(OMEGA_hat, type="F")
    list( FN = FNrate(OMEGA_hat, Omega)
          ,FP = FPrate(OMEGA_hat, Omega)
          ,TP = TPrate(OMEGA_hat, Omega)
          ,TN = TNrate(OMEGA_hat, Omega)
          ,MAXInfOmega = max(abs(OMEGA_hat-Omega))
          ,ForbNormOmega = norm(OMEGA_hat-Omega, typ="F")
          ,MAXInfSigma = max(abs(SIGMA_hat-Solved_Omega))
          ,ForbNormSigma = norm(SIGMA_hat-Solved_Omega, typ="F")
    )
  }
    # toc()
  
  stopCluster(c1)
  # cat( text1,
  #      "\n", paste(p_conta*100, "%", sep=""), "Max infinity norm Sigma: ", mean(map_dbl(Results, \(x) x[["MAXInfSigma"]]))
  #      ,paste0("(", formatC(sd(map_dbl(Results, \(x) x[["MAXInfSigma"]]))/sqrt(R), format="e", digits=5), ")")
  #      
  #      ,"\n", paste(p_conta*100, "%", sep=""), "Frob. Norm Sigma: ", mean(map_dbl(Results, \(x) x[["ForbNormSigma"]]))
  #      ,paste0("(", formatC(sd(map_dbl(Results, \(x) x[["ForbNormSigma"]]))/sqrt(R), format="e", digits=5), ")")
  #      
  #      ,"\n", paste(p_conta*100, "%", sep=""), "Max infinity norm Omega: ", mean(map_dbl(Results, \(x) x[["MAXInfOmega"]]))
  #      ,paste0("(", formatC(sd(map_dbl(Results, \(x) x[["MAXInfOmega"]]))/sqrt(R), format="e", digits=5), ")")
  #      
  #      ,"\n", paste(p_conta*100, "%", sep=""), "Frob. Norm Omega: ", mean(map_dbl(Results, \(x) x[["ForbNormOmega"]]))
  #      ,paste0("(" , formatC(sd(map_dbl(Results, \(x) x[["ForbNormOmega"]]))/sqrt(R), format="e", digits=5), ")")
  #      
  #      
  #      ,"\n", paste(p_conta*100, "%", sep=""), "contamination FP: ", mean(map_dbl(Results, \(x) x[["FP"]]))
  #      ,paste0("(", formatC(sd(map_dbl(Results, \(x) x[["FP"]]))/sqrt(R), format="e", digits=5), ")")
  #      
  #      
  #      ,"\n", paste(p_conta*100, "%", sep=""), "contamination FN: ", mean(map_dbl(Results, \(x) x[["FN"]]))
  #      ,paste0("(", formatC(sd(map_dbl(Results, \(x) x[["FN"]]))/sqrt(R), format="e", digits=5), ")")
  #      
  #      
  #      ,"\n", paste(p_conta*100, "%", sep=""), "contamination TP: ", mean(map_dbl(Results, \(x) x[["TP"]]))
  #      ,paste0("(", formatC(sd(map_dbl(Results, \(x) x[["TP"]]))/sqrt(R), format="e", digits=5), ")")
  #      
  #      
  #      ,"\n", paste(p_conta*100, "%", sep=""), "contamination TN: ", mean(map_dbl(Results, \(x) x[["TN"]]))
  #      ,paste0("(", formatC(sd(map_dbl(Results, \(x) x[["TN"]]))/sqrt(R), format="e", digits=5), ")"),  "\n",
  #      text2, "\n"
  # )
}

R <- 1
n <- 200
p <- 200
p_conta <- 0.05
# cat("Computation time for k-Tau\n\n")
tic("Computation time for k-Tau\n\n")
# time.ktau <- system.time(
doSimulation2(Grho, FALSE, FALSE, lower = 0.01, p = p, n=n, R=R, p_conta=p_conta, robust_KTau, generate_sparse_omega, medDev,
              rmvnorm_precision_n, TPrate, TNrate, FNrate, FPrate,
              paste0("===================================== START:::::::::", p_conta1*100,  "% cellwise n = ", n, " p = ", p,
                      "Kendall, SPARSE: ChainOmega", "With ", R, " Replications"),
              paste0("===================================== END:::::::::::", p_conta*100,  "% cellwise n = ", n, " p = ", p,
                       "Kendall, SPARSE: ChainOmega", "With ", R, " Replications")
)
# )
toc()
# print(time.ktau)

tic("Computation time for NS-Estimator\n\n")
# time.NS <- system.time(
doSimulation2(Grho, FALSE, TRUE, lower = 0.01, p = p, n=n, R=R, p_conta=p_conta, sampleCOVest, generate_sparse_omega, medDev,
              rmvnorm_precision_n, TPrate, TNrate, FNrate, FPrate,
              paste0("===================================== START:::::::::", p_conta*100,  "% cellwise n = ", n, " p = ", p,
                      "Kendall, SPARSE: ChainOmega", "With ", R, " Replications"),
              paste0("===================================== END:::::::::::", p_conta*100,  "% cellwise n = ", n, " p = ", p,
                      "Kendall, SPARSE: ChainOmega", "With ", R, " Replications")
)
# )
toc()
# cat("Computation time for NS-Estimator\n\n")
# print(time.NS)