
x <- x
val <- vax
est.mode <- NULL
lambda.vec <- NULL
lambda.list <- NULL
Omegatilde.list <- NULL
scale.vec <- NULL
normalize <- TRUE
thres <- 1.0e-4
maxit <- 1e4
njobs <- 4
Grho <- NULL

# proper candidates of tuning parameters
lamseq <- seq(1.5e-2, 1.5, length.out = 100)
lambda.list <- list() # a list containing candidates of tuning parameters for each mode
for (i in 1:K) {
  lambda.list[[i]] <- lamseq
}

lamseq.C <- seq(1.5e-3, 3, length.out = 100)
lambda.list.C <- list() # a list containing candidates of tuning parameters for each mode
for (i in 1:K) {
  lambda.list.C[[i]] <- lamseq.C
}

mode_index <- 1
