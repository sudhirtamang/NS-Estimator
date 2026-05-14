# ==============================================================================
# PROJECT:  NS-Estimator
# SCRIPT:   Lfunctions
# PURPOSE:  Functions for the NS-Estimator project
# AUTHOR:   Sudhir T
# DATE:     15th April, 2026
# ==============================================================================

omega_from_sample <- function(x, dimen, k, n, normalize=TRUE){
  K <- length(dimen)
  nvars <- prod(dimen)
  njobs <- ceiling(0.65 * parallel::detectCores())
  c1 = parallel::makeCluster(ifelse(K <= njobs, K, njobs))
  doParallel::registerDoParallel(c1)
  fit1 = foreach::foreach(k = 1:K, .combine = list,  .multicombine = TRUE) %dopar% {
    if (n * nvars < ((dimen[k]**2) * (dimen[k] - 1) / 2)) {
      Omega_tilde = diag(dimen[k])
    }else {
      S.array = array(0, c(dimen[k], dimen[k], n))
      for (i in 1:n) {
        Vi = rTensor::k_unfold(rTensor::as.tensor(x[, , , i]), m = k)@data  
        S.array[, , i] = Vi %*% t(Vi)
      }
      S.mat = apply(S.array, c(1, 2), mean) * dimen[k] / prod(dimen) 
      Omega_tilde = solve(S.mat)
      if (normalize) {
        Omega_tilde = Omega_tilde / norm(Omega_tilde, type = "F")}
    }
    Omega_tilde
  }
  parallel::stopCluster(c1)
  fit1
}

# =============================>>>>>>>>>>>>>>>>>>>>
tilde_Sk <- function(x, Omega, dimen, k, n){
  K <- length(dimen)
  Omega.list.sqrt <- purrr::map(Omega, sqrtm)
  Omega.list.sqrt[[k]] <- diag(dimen[[k]])
  tmparr <- array(0, c(dimen[[k]], dimen[[k]], n))
  for (i in 1:n) {
    d <- do.call("[", c(list(x), rep(list(substitute()), K), i))
    Vi <- k_unfold(as.tensor(ttl(as.tensor(d), Omega.list.sqrt,
                                ms = 1:length(dimen)
    )@data), m = k)@data
    tmparr[, , i] <- Vi %*% t(Vi)
  }
  apply(tmparr, c(1, 2), mean) * (dimen[k] / prod(dimen))
}

# ============================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

fast_sep <- function(x, dimen, n, fit1, lam.best, normalize=TRUE, thres=1.0e-4, maxit=1e4){
  K <- length(dimen)
  nvars <- prod(dimen)
  njobs <- ceiling(0.65 * parallel::detectCores())
  c1 = parallel::makeCluster(ifelse(K <= njobs, K, njobs))
  doParallel::registerDoParallel(c1)
  fit_result = foreach(k = 1:K, .packages = c("glasso", "rTensor", "expm"), .combine = list, .multicombine = TRUE) %dopar% {
    Omega.list.sqrt = list()
    for (i in 1:K) {
      Omega.list.sqrt[[i]] = sqrtm(fit1[[i]])
    }
    S.array = array(0, c(dimen[k], dimen[k], n))
    Omega.list.sqrt[[k]] = diag(dimen[k])
    for (i in 1:n) {
      Vi = rTensor::k_unfold(rTensor::as.tensor(rTensor::ttl(rTensor::as.tensor(x[, , , i]), Omega.list.sqrt,
                                  ms = 1:K
      )@data), m = k)@data
      S.array[, , i] = Vi %*% t(Vi)
    }
    S.mat = apply(S.array, c(1, 2), mean) * dimen[k] / prod(dimen)
    Out1 = glasso::glasso(S.mat, rho = lam.best[k], penalize.diagonal = FALSE, maxit = maxit, thr = thres)
    hat_Omega = as.matrix(Out1$wi)
    if (normalize) {
      hat_Omega = hat_Omega / norm(hat_Omega, type = "F")
    }
    hat_Omega
  }
  stopCluster(c1)
  fit_result
}


# ===================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ===================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ===================================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# x[, , , 1] is the first tensor, dimen should be din(x[, , , 1])
NSEstimator <- function(x, dimen){
  tmp <- array(NA, dim(x))
  for(j in 1:dimen[[1]]){
    for(k in 1:dimen[[2]]){
      for(l in 1:dimen[[3]]){
        tmpCDF <- stats::ecdf(x[j, k, l, ])
        tmp[j, k, l, ] <- qnorm((n/(n+1)) * tmpCDF(x[j, k, l, ]))
      }
    }
  }
  return(tmp)
}


NSEstimator2 <- function(x, dimen){
  tmp <- array(NA, dim(x))
  for(j in 1:dimen[[1]]){
    for(k in 1:dimen[[2]]){
        tmpCDF <- stats::ecdf(x[j, k, ])
        tmp[j, k, ] <- qnorm((n/(n+1)) * tmpCDF(x[j, k, ]))
    }
  }
  return(tmp)
}


tildeOmega <- function(x, dimen, n){
  # x is an array.
  # if dim(x) = c(30, 34, 35, 200) then the first tensor observation is x[, , , 1].
  # dimen is the dimension of each tensor observation
  # n is the sample size, or the length of the last dimension of x, for the e.g. before it is 200.
  K <- length(dimen)
  nvars <- prod(dimen)
  Omega_tilde <- vector("list", K)
  for(k in 1:K){
  # fit1 = foreach(k = 1:K, .export = c("x"), .combine = list, .multicombine = TRUE) %dopar% {
    # when sample size is small, use the identity matrix
    if (n * nvars < ((dimen[k]**2) * (dimen[k] - 1) / 2)) {
      Omega_tilde[[k]] <- diag(dimen[k])
    }else {
      # when sample size is large, calculate the sample estimator of the precision matrices
      S.array <- array(dim=c(dimen[k], dimen[k], n))
      for (i in 1:n) {
        d <- do.call("[", c(list(x), rep(list(substitute()), K), i))
        Vi <- rTensor::k_unfold(rTensor::as.tensor(d), m = k)@data  # unfold tensor
        S.array[, , i] <- Vi %*% t(Vi)
      }
      S.mat <- apply(S.array, c(1, 2), mean) * dimen[k] / nvars 
      Omega_tilde[[k]] <- solve(S.mat)
      Omega_tilde[[k]] <- Omega_tilde[[k]] / norm(Omega_tilde[[k]], type = "F")}
  }
  Omega_tilde
}


# 
# fold_indices <- caret::createFolds(1:n, k = 5, list = TRUE)
# fold_SIGMA_train <- vector("list", 5)
# fold_SIGMA_test  <- vector("list", 5)
# 
# for (j in seq_along(fold_indices)) {
#   test  <- data[fold_indices[[j]], ]
#   train <- data[-fold_indices[[j]], ]
#   
#   hat_sigma_test  <- apply(test, 2, medDev) / q75
#   hat_sigma_train <- apply(train, 2, medDev) / q75
#   
#   fold_SIGMA_test[[j]]  <- estSIGMA_hat(test, hat_sigma_test)
#   fold_SIGMA_train[[j]] <- estSIGMA_hat(train, hat_sigma_train)
# }
# # --- Lambda loop ---
# NegLoglikelihood <- numeric(length(LAMBDAs))
# for (i in seq_along(LAMBDAs)) {
#   tmp0 <- numeric(5)
#   rho_matrix <- matrix(LAMBDAs[[i]], nrow = p, ncol = p)
#   diag(rho_matrix) <- 0
#   for (j in seq_along(fold_indices)) {
#     fit <- QUIC(fold_SIGMA_train[[j]], rho_matrix,
#                 tol = 1e-4, maxIter = 1000, msg = 0, path = NULL)
#     Theta <- matrix(fit$X, p, p)
#     R1 <- chol(Theta)
#     log_det <- 2 * sum(log(diag(R1)))
#     
#     trace   <- sum(fold_SIGMA_test[[j]] * Theta)
#     tmp0[j] <- -log_det + trace
#   }
#   NegLoglikelihood[[i]] <- mean(tmp0)
# }




# Setup: A 3D array (could be any dimension)
my_array <- array(1:24, dim = c(2, 3, 4))

# Goal: Subset the 2nd index of the 2nd dimension
target_dim <- 2
target_index <- 2

# 1. Create a list of 'empty' arguments for every dimension
args <- rep(list(substitute()), length(dim(my_array)))

# 2. Fill in the specific dimension you want to subset
args[[target_dim]] <- target_index

# 3. Use do.call to call the "[" function
# Equivalent to: my_array[, 2, ]
result <- do.call("[", c(list(my_array), args))

