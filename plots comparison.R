sfit <- Separate.fit(x, vax, lambda.list = lambda.list)
sfit1 <- Separate.fit(x, vax, lambda.list = lambda.list, Grho=Grho)



Tfit <- Separate.fit(Tx, Tvax, lambda.list = lambda.list)
Tfit1 <- Separate.fit(Tx, Tvax, lambda.list = lambda.list, Grho=Grho)



plot(c(Sigma[[1]]), c(Sigma[[1]]), main="Compare Sigma", cex=0.5, xlab="True Sigma", ylab="Estimated Sigmas")
points(c(Sigma[[1]]), c(sfit[["Omegahat"]][[1]][[2]]), col="red", asp=1, cex=0.5)
points(c(Sigma[[1]]), c(Tfit[["Omegahat"]][[1]][[2]]), col="yellow", cex=0.5)
points(c(Sigma[[1]]), c(Tfit1[["Omegahat"]][[1]][[2]]), col="green", cex=0.5)
legend(
  "topleft",
  legend = c("True", "sfit", "Tfit", "Tfit1"),
  col = c("black", "red", "yellow", "green"),
  pch = 16
  # bty = "n"
)

plot(c(Omega[[1]]), c(Omega[[1]]), main="Compare Omega", cex=0.5, xlab="True Omega", ylab="Est. Omega")
points(c(Omega[[1]]), c(sfit[["Omegahat"]][[1]][[1]]), col="red", asp=1, cex=0.5)
points(c(Omega[[1]]), c(Tfit[["Omegahat"]][[1]][[1]]), col="yellow", cex=0.5)
points(c(Omega[[1]]), c(Tfit1[["Omegahat"]][[1]][[1]]), col="green", cex=0.5)
legend(
  "topleft",
  legend = c("True", "sfit", "Tfit", "Tfit1"),
  col = c("black", "red", "yellow", "green"),
  pch = 16
  # bty = "n"
)



simulation.summary(list(Tfit[["Omegahat"]][[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
simulation.summary(list(sfit[["Omegahat"]][[1]][[2]]), list(Sigma[[1]]), offdiag = FALSE)
