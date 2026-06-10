n<-50
p<-30

g<-function(n,rho=seq(-1,1,0.01),B=1000){
  n.rho<-length(rho)
  g.value<-rep(0,n.rho)
  for(i.rho in 1:n.rho){
    rho.tmp<-rho[i.rho]
    g.value.b<-rep(0,B)
    for(b in 1:B){
      W<-matrix(0,n,2)
      W[,1]<-rnorm(n)
      W[,2]<-rho.tmp*W[,1]+sqrt(1-rho.tmp^2)*rnorm(n)
      W.T<-matrix(0,n,2)
      W.T[,1]<-qnorm(rank(W[,1])/(n+1))
      W.T[,2]<-qnorm(rank(W[,2])/(n+1))
      g.value.b[b]<-mean(W.T[,1]*W.T[,2])
    }
    g.value[i.rho]<-mean(g.value.b)
  }
  object<-list(rho=rho,g.value=g.value)
  object
}



sigma<-matrix(0,p,p)
r<-0.9
for(i in 1:p){
  for(j in 1:p){
    sigma[i,j]<-r^abs(i-j)
  }
}

#CS
#r<-0.8
#sigma<-matrix(r,p,p)
#diag(sigma)<-1


sigma.eigen<-eigen(sigma)
sigma.sqrt<-sigma.eigen$vectors%*%diag(sqrt(sigma.eigen$values))%*%t(sigma.eigen$vectors)



obj<-g(n)

n.rep<-100
bias.correct<-rep(0,n.rep)
bias.T<-rep(0,n.rep)
bias.ave<-rep(0,n.rep)
error.correct<-rep(0,n.rep)
error.T<-rep(0,n.rep)
error.ave<-rep(0,n.rep)


for(i.rep in 1:n.rep){

  X<-matrix(rnorm(n*p),n,p)
  X<-X%*%sigma.sqrt  
  
X.T<-matrix(0,n,p)
for(j in 1:p){
  X.T[,j]<-qnorm(rank(X[,j])/(n+1))
}

sigma.hat<-t(X)%*%X/n
sigma.hat.T<-t(X.T)%*%X.T/n
g.value<-obj$g.value
rho<-obj$rho

sigma.hat.correct<-matrix(0,p,p)

for(i in 1:p){
  for(j in 1:p){
    sigma.hat.correct[i,j]<-rho[which.min(abs(g.value-sigma.hat.T[i,j]))]
  }
}

bias.correct[i.rep]<-mean(sigma.hat.correct[1,2]-sigma[1,2])
bias.T[i.rep]<-mean(sigma.hat.T[1,2]-sigma[1,2])
bias.ave[i.rep]<-mean(sigma.hat[1,2]-sigma[1,2])

#error.correct[i.rep]<-max(abs(sigma.hat.correct-sigma))
#error.T[i.rep]<-max(abs(sigma.hat.T-sigma))
#error.ave[i.rep]<-max(abs(sigma.hat-sigma))

error.correct[i.rep]<-mean((sigma.hat.correct-sigma)^2)
error.T[i.rep]<-mean(abs(sigma.hat.T-sigma)^2)
error.ave[i.rep]<-mean(abs(sigma.hat-sigma)^2)
}

par(mfrow=c(1,2))
boxplot(bias.correct,bias.T,bias.ave)
boxplot(error.correct,error.T,error.ave)
