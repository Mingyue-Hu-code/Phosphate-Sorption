data<-read.csv("projdat.txt",sep = "")

X<-data$x
Y<-data$y

mle <- function(par, x,y) {
  
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  theta <- par[4]
  sigma <-par[5]
  
  mu <- alpha*x^(beta*x^(-gamma))
  
  l <- -length(x)/2*log(2*pi*sigma)-1/(2*sigma)*((y[1]-mu[1])^2+sum((y[-1]-mu[-1]-theta*(y[-25]-mu[-25]))^2))
  
  return(-l)
}

##### maximum log likelihood

initial <- c(0.5063, 1.222, 0.1614, 0, 1.0149)


model1 <- optim(initial, mle, x=X, y=Y, hessian = T, method = "BFGS")
est1<-model1$par

#confidence band
covb<-solve(model1$hessian)[1:3,1:3]
covb1<-solve(model1$hessian)

0.7322784-(1.96)*sqrt(covb[1,1])
0.7322784+(1.96)*sqrt(covb[1,1])

1.0335912-(1.96)*sqrt(covb[2,2])
1.0335912+(1.96)*sqrt(covb[2,2])

0.1575802-(1.96)*sqrt(covb[3,3])
0.1575802+(1.96)*sqrt(covb[3,3])

0.5584688-(1.96)*sqrt(covb1[4,4])
0.5584688+(1.96)*sqrt(covb1[4,4])

0.6244067-(1.96)*sqrt(covb1[5,5])
0.6244067+(1.96)*sqrt(covb1[5,5])


#plot
fctn<-function(xmat, cps){
  cps[1]*xmat^(cps[2]*xmat^(-cps[3]))
}

ders<-function(xmat, cps){
  cbind(xmat^(cps[2]*xmat^(-cps[3])), cps[1]*xmat^(cps[2]*xmat^(-cps[3])-cps[3])*log(xmat), 
        -cps[1]*cps[2]*xmat^(-cps[3])*(log(xmat))^2*xmat^(cps[2]*xmat^(-cps[3])))
}

vhatmufun <- function(par, derfun, X, vcov){
  diag(derfun(X,par)%*%vcov%*%t(derfun(X,par)))
}

wts<-function(xmat, cps){
  diag(1,nrow = 25,ncol = 25)
}

sehatmu <- sqrt(vhatmufun(c(0.7322784, 1.0335912, 0.1575802), ders, data$x, covb))


plot(data$x, data$y, xlab = "covariate values", ylab = "response variables")
lines(data$x[order(data$x)], fctn(data$x,c(0.7322784, 1.0335912, 0.1575802))[order(data$x)])
lines(data$x[order(data$x)], fctn(data$x,c(0.7322784, 1.0335912, 0.1575802))[order(data$x)] + 1.96*sehatmu[order(data$x)], lty = 2)
lines(data$x[order(data$x)], fctn(data$x,c(0.7322784, 1.0335912, 0.1575802))[order(data$x)] - 1.96*sehatmu[order(data$x)], lty = 2)


####MLE 2

mle2 <- function(par, x,y) {
  
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  theta <- par[4]
  sigma <-par[5]
  
  mu <- alpha*x^(beta*x^(-gamma))
  
  l <- -(length(x)-1)/2*log(2*pi*sigma)-1/(2*sigma)*(sum((y[-1]-mu[-1]-theta*(y[-25]-mu[-25]))^2))
  
  return(-l)
}

##### maximum log likelihood

initial <- c(0.5063, 1.222, 0.1614, 0, 1.0149)


model2 <- optim(initial, mle2, x=X, y=Y, hessian = T, method = "BFGS")
model2
est2<-model2$par
#### confidence band
#confidence band (jackknife)
C1<-rep(0,5)
for (i in 1:25) {

model_j<- optim(initial, mle2, x=X[-i], y=Y[-i], hessian = T, method = "BFGS")
C1<-C1+(est2-model_j$par)^2

}

C1<-(24/25)*C1;C1


0.7322784-(1.96)*sqrt(covb[1,1])
0.7322784+(1.96)*sqrt(covb[1,1])


####plot
sehatmu <- sqrt(vhatmufun(c(0.7322784, 1.0335912, 0.1575802), ders, data$x, covb))


plot(data$x, data$y, xlab = "covariate values", ylab = "response variables")
lines(data$x[order(data$x)], fctn(data$x,c(0.7322784, 1.0335912, 0.1575802))[order(data$x)])
lines(data$x[order(data$x)], fctn(data$x,c(0.7322784, 1.0335912, 0.1575802))[order(data$x)] + 1.96*sehatmu[order(data$x)], lty = 2)
lines(data$x[order(data$x)], fctn(data$x,c(0.7322784, 1.0335912, 0.1575802))[order(data$x)] - 1.96*sehatmu[order(data$x)], lty = 2)

