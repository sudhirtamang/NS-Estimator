

source("C:/Users/sudhi/Desktop/Coordinatewise Gaussianization/replication/List_functions_for_cellwise_corruption.R")
Model <- function(n, seed, dimen) {
  
  nvars <- prod(dimen) # number of variables
  
  
  Omega <- generate_sparse_omega(dimen[[1]], dimen[[1]])
  Sigma <- solve(Omega)

  
  
  set.seed(seed) 
  
  # Generate data observation
  # training set
  x <- rmvnorm_precision_n(n, dimen[[1]], Omega)
  
  # validation set
  vax <- rmvnorm_precision_n(n, dimen[[1]], Omega)
  
  result <- list()
  result$x <- x
  result$vax <- vax
  
  return(list(result, Sigma, Omega))
}



