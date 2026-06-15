sfit <- Separate.fit(x, vax, lambda.list = lambda.list)
sfit1 <- Separate.fit(x, vax, lambda.list = lambda.list, Grho=Grho)



Tfit <- Separate.fit(Tx, Tvax, lambda.list = lambda.list)
Tfit1 <- Separate.fit(Tx, Tvax, lambda.list = lambda.list, Grho=Grho)


norm(sfit[["Omegahat"]][[1]][[2]] - Sigma[[1]], type="F")
norm(Tfit[["Omegahat"]][[1]][[2]] - Sigma[[1]], type="F")
norm(Tfit1[["Omegahat"]][[1]][[2]] - Sigma[[1]], type="F")



plot(log(c(Sigma[[1]])), log(c(Sigma[[1]])), main="Compare Sigma", cex=0.5, xlab="True Sigma", ylab="Estimated Sigmas", pch=16)
abline(0, 1)
points(log(c(Sigma[[1]])), log(c(sfit[["Omegahat"]][[1]][[2]])), col="red", asp=1, cex=0.5, pch=16)
points(log(c(Sigma[[1]])), log(c(Tfit[["Omegahat"]][[1]][[2]])), col="yellow", cex=0.5, pch=16)
points(log(c(Sigma[[1]])), log(c(Tfit1[["Omegahat"]][[1]][[2]])), col="green", cex=0.5, pch=16)
# points(log(c(Tfit[["Omegahat"]][[1]][[2]])), log(c(Tfit1[["Omegahat"]][[1]][[2]])), col="blue", cex=0.5, pch=16)
legend(
  "topleft",
  legend = c("True", "sfit", "Tfit", "Tfit1", "asc"),
  col = c("black", "red", "yellow", "green", "blue"),
  pch = 16
  # bty = "n"
)


plot((c(Sigma[[1]])), (c(Sigma[[1]])), main="Compare Sigma", cex=0.5, xlab="True Sigma", ylab="Estimated Sigmas", pch=16)
abline(0, 1)
points((c(Sigma[[1]])), (c(sfit[["Omegahat"]][[1]][[2]])), col="red", asp=1, cex=0.5, pch=16)
points((c(Sigma[[1]])), (c(Tfit[["Omegahat"]][[1]][[2]])), col="yellow", cex=0.5, pch=16)
points((c(Sigma[[1]])), (c(Tfit1[["Omegahat"]][[1]][[2]])), col="green", cex=0.5, pch=16)
points((c(Tfit[["Omegahat"]][[1]][[2]])), (c(Tfit1[["Omegahat"]][[1]][[2]])), col="blue", cex=0.5, pch=16)
# plot((c(Tfit[["Omegahat"]][[1]][[2]])), (c(Tfit1[["Omegahat"]][[1]][[2]])), col="blue", cex=0.5, pch=16)
# abline(0, 1)
legend(
  "topleft",
  legend = c("True", "sfit", "Tfit", "Tfit1", "asc"),
  col = c("black", "red", "yellow", "green", "blue"),
  pch = 16
  # bty = "n"
)



idx.zeros <- which(c(Omega[[1]]) == 0)
y <- c(Omega[[1]])[idx.zeros]
plot(idx.zeros, y, main="Compare Omega", cex=0.5, xlab="True Omega", ylab="Est. Omega", pch=16, ylim=c(-1, 4))

idx <- which(c(sfit[["Omegahat"]][[1]][[1]])[idx.zeros] != 0)
print(length(idx))
y <- (c(sfit[["Omegahat"]][[1]][[1]])[idx.zeros])[idx]
points(idx, y + 1, col="red", cex=0.2, pch=16)

idx <- which(c(Tfit[["Omegahat"]][[1]][[1]])[idx.zeros] != 0)
print(length(idx))
y <- (c(Tfit[["Omegahat"]][[1]][[1]])[idx.zeros])[idx]
points(idx, y + 2, col="blue", cex=0.2, pch=16)

idx <- which(c(Tfit1[["Omegahat"]][[1]][[1]])[idx.zeros] != 0)
print(length(idx))
y <- (c(Tfit1[["Omegahat"]][[1]][[1]])[idx.zeros])[idx]
points(idx, y + 3, col="green", cex=0.2, pch=16)
legend(
  "topleft",
  legend = c("True", "sfit", "Tfit", "Tfit1"),
  col = c("black", "red", "blue", "green"),
  pch = 16
  # bty = "n"
)


tmp1 <- Tfit1[["Omegahat"]][[1]][[2]]

sum(c(tmp1) < 1e-6)
simulation.summary(list(tmp1), list(Omega[[1]]), offdiag = FALSE)
simulation.summary(list(sfit[["Omegahat"]][[1]][[1]]), list(Omega[[1]]), offdiag = FALSE)






