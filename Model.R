
Model <- function(n, seed, dimen) {
  
  nvars <- prod(dimen) # number of variables
  K <- length(dimen) # order of X
  
  
  Omega <- purrr::map2(dimen, 1:K, \(x, y) ChainOmega(x, sd = y*100, norm.type = 2))
  Omega[[2]] <- diag(dimen[[2]])
  Sigma <- purrr::map(Omega, \(x) solve(x))
  multiples <- purrr::map_dbl(Sigma, \(x) x[1, 1])
  Sigma <- purrr::map(Sigma, \(x) x/x[1, 1]) # covariance matrix
  Omega <- purrr::map2(Omega, multiples, \(x, y) x*y)
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