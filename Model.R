
Model <- function(n, seed, dimen) {
  
  nvars <- prod(dimen) # number of variables
  tensor.order <- length(dimen) # order of X
  
  
  Omega <- purrr::map2(dimen, 1:tensor.order, \(p, i) ChainOmega(p, sd = i*100, norm.type = 2))
  for(i in seq_along(dimen[2:tensor.order])){
    Omega[[i+1]] <- diag(dimen[[i + 1]])
  }
  Sigma <- purrr::map(Omega, \(omega) solve(omega))
  multiples <- purrr::map_dbl(Sigma, \(sigma) sigma[1, 1])
  Sigma <- purrr::map(Sigma, \(sigma) sigma/sigma[1, 1]) # covariance matrix
  Omega <- purrr::map2(Omega, multiples, \(omega, c) omega*c)
  dSigma <- purrr::map(Sigma, \(sigma) t(chol(sigma))) # square root of covariance matrix
  
  
  set.seed(seed) 
  
  # Generate data observation
  # training set
  vec_x <- matrix(rnorm(nvars * n), ncol = n) 
  x <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    d <- array(vec_x[, i], dimen)
    x <- do.call("[<-", c(list(x), rep(list(substitute()), tensor.order), list(i, atrans(d, dSigma))))
  }
  
  # validation set
  vec_vax <- matrix(rnorm(nvars * n), ncol = n) 
  vax <- array(0, dim = c(dimen, n))
  for (i in 1:n) {
    d <- array(vec_vax[, i], dimen)
    vax <- do.call("[<-", c(list(vax), rep(list(substitute()), tensor.order), list(i, atrans(d, dSigma))))
  }
  
  result <- list()
  result$x <- x
  result$vax <- vax
  
  return(list(result, Sigma, Omega))
}



