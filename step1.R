source("gfunctionsforFishExample.R")
source("glsfunctions.R")
source("part1.R")

data<-read.csv("projdat.txt",sep = "")

plot(data$x, data$y, xlab = "covariate values", ylab = "response variables")

fctn<-function(xmat, cps){
  cps[1]*xmat^(cps[2]*xmat^(-cps[3]))
}

ders<-function(xmat, cps){
  cbind(xmat^(cps[2]*xmat^(-cps[3])), cps[1]*xmat^(cps[2]*xmat^(-cps[3])-cps[3])*log(xmat), 
        -cps[1]*cps[2]*xmat^(-cps[3])*(log(xmat))^2*xmat^(cps[2]*xmat^(-cps[3])))
}

wts<-function(xmat, cps){
  diag(1,nrow = 25,ncol = 25)
}

cps<-gauss.newton(data[,1],data[,2],c(-1.4675,0.7038,0.001),fctn, ders, wts)
a<-nonlin(data[,1],data[,2],c(-1.4675,0.7038,0.001), fctn, ders, wts)

0.506264-(1.96)*sqrt(0.18229784)
0.506264+(1.96)*sqrt(0.18229784)

1.221994-(1.96)*sqrt(0.26592934)
1.221994+(1.96)*sqrt(0.26592934)

0.1613655-(1.96)*sqrt(0.0005802841)
0.1613655+(1.96)*sqrt(0.0005802841)

ci<-c(sqrt(0.18229784),sqrt(0.26592934),sqrt(0.0005802841))

a<-lm(log(data[,2])~log(data[,1]))
summary(a)
exp((-1.4675)/0.7038)

vhatmufun <- function(par, derfun, X, vcov){
  diag(derfun(X,par)%*%vcov%*%t(derfun(X,par)))
}

sehatmu <- sqrt(vhatmufun(c(0.506264, 1.221994, 0.1613655), ders, data$x, a[["covb"]]))


#plot
plot(data$x, data$y, xlab = "covariate values", ylab = "response variables")
lines(data$x[order(data$x)], fctn(data$x,c(0.506264, 1.221994, 0.1613655))[order(data$x)])
lines(data$x[order(data$x)], fctn(data$x,c(0.506264, 1.221994, 0.1613655))[order(data$x)] + 1.96*sehatmu[order(data$x)], lty = 2)
lines(data$x[order(data$x)], fctn(data$x,c(0.506264, 1.221994, 0.1613655))[order(data$x)] - 1.96*sehatmu[order(data$x)], lty = 2)

#residual
res<-(data$y-fctn(data$x,c(0.506264, 1.221994, 0.1613655)))/sqrt(fctn(data$x,c(0.506264, 1.221994, 0.1613655)))
plot(data$x, a[["stdres"]], xlab = "covariate values", ylab = "standardized residuals")
abline(a=0, b=0)

#question 2
w<-data$y-fctn(data$x,c(0.506264, 1.221994, 0.1613655))
w1<-w[-1]
w2<-w[-25]
wi<-sum(w1*w2)
wj<-sum(w2^2)
theta1<-wi/wj
x<-w1-theta1*w2
sigma1<-(1/24)*sum(x^2)

#Question 2e
B <- 2500
b <- 0
theta <- c()
sigma <- c()

repeat{
  
  b <- b + 1
  
  ws <- arima.sim(list(order = c(1,0,0), ar = theta1), n = 25, rand.gen= rnorm, sd = sqrt(sigma1)) 
  
  
  ws1<-ws[-1]
  ws2<-ws[-25]
  wsi<-sum(ws1)
  wsj<-sum(ws2)
  thetahat<-wsi/wsj
  xs<-ws1-thetahat*ws2
  sigmahat<-(1/24)*sum(xs^2)
  
  theta <- c(theta, thetahat)
  sigma <- c(sigma, sigmahat)
  
  if(b == B){break}
  
}


biashat<-mean(theta1-theta)
bias<-theta1-theta

lowerEndpoint <- 2*theta1-quantile(theta, prob = 0.975)
upperEndpoint <- 2*theta1-quantile(theta, prob = 0.025)
c(lowerEndpoint, upperEndpoint)


#####  Percentile Bootstrap confidence interval:

c(quantile(bias, prob = 0.025), quantile(bias, prob = 0.975))


#####3

repeat{
  
  b <- b + 1
  
  ws<-rep(0,25)
  ws[1]<-w[1]
  for (i in 2:25) {
    ws[i]<-theta1*ws[i-1]+rnorm(1,mean = 0, sd = sqrt(sigma1))
  } 
  
  
  ws1<-ws[-1]
  ws2<-ws[-25]
  wsi<-sum(ws1)
  wsj<-sum(ws2)
  thetahat<-wsi/wsj
  xs<-ws1-thetahat*ws2
  sigmahat<-(1/24)*sum(xs^2)
  
  theta <- c(theta, thetahat)
  sigma <- c(sigma, sigmahat)
  
  if(b == B){break}
  
}

