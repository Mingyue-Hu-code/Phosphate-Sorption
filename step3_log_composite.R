
data<-read.csv("projdat.txt",sep = "")


fctn<-function(xmat, cps){
  cps[1]*xmat^(cps[2]*xmat^(-cps[3]))
}


w<-data$y-fctn(data$x,c(0.506264, 1.221994, 0.1613655))
w1<-w[-1]
w2<-w[-25]
wi<-sum(w1*w2)
wj<-sum(w2^2)
theta1<-wi/wj
x<-w1-theta1*w2
sigma1<-(1/24)*sum(x^2)

B <- 2500
b <- 0
theta <- c()
sigma <- c()

repeat{
  
  b <- b + 1
  
  ws<-rep(0,25)
  ws[1]<-w[1]
  for (i in 2:25) {
    ws[i]<-theta1*ws[i-1]+rnorm(1,mean = 0, sd = sqrt(sigma1))
  } 
  
  
  ws1<-ws[-1]
  ws2<-ws[-25]
  wsi<-sum(ws1*ws2)
  wsj<-sum(ws2^2)
  thetahat<-wsi/wsj
  xs<-ws1-thetahat*ws2
  sigmahat<-(1/24)*sum(xs^2)
  
  theta <- c(theta, thetahat)
  sigma <- c(sigma, sigmahat)
  
  if(b == B){break}
  
}

biashat<-mean(theta1-theta)
#bias<-theta1-theta

lowerEndpoint <- 2*theta1-quantile(theta, prob = 0.975)
upperEndpoint <- 2*theta1-quantile(theta, prob = 0.025)
c(lowerEndpoint, upperEndpoint)


#####  Percentile Bootstrap confidence interval:

c(quantile(theta, prob = 0.025), quantile(theta, prob = 0.975))

######for sigma2
biashat<-mean(sigma1-sigma)

lowerEndpoint <- 2*sigma1-quantile(sigma, prob = 0.975)
upperEndpoint <- 2*sigma1-quantile(sigma, prob = 0.025)
c(lowerEndpoint, upperEndpoint)


c(quantile(sigma, prob = 0.025), quantile(sigma, prob = 0.975))



#for naive
w<-data$y-fctn(data$x,c(0.506264, 1.221994, 0.1613655))
w1<-w[-1]
w2<-w[-25]
wi<-sum(w1)
wj<-sum(w2)
theta1<-wi/wj
x<-w1-theta1*w2
sigma1<-(1/24)*sum(x^2)

B <- 2500
b <- 0
theta <- c()
sigma <- c()

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

biashat<-mean(theta1-theta)
#bias<-theta1-theta

lowerEndpoint <- 2*theta1-quantile(theta, prob = 0.975)
upperEndpoint <- 2*theta1-quantile(theta, prob = 0.025)
c(lowerEndpoint, upperEndpoint)


#####  Percentile Bootstrap confidence interval:

c(quantile(theta, prob = 0.025), quantile(theta, prob = 0.975))

######for sigma2
biashat<-mean(sigma1-sigma)

lowerEndpoint <- 2*sigma1-quantile(sigma, prob = 0.975)
upperEndpoint <- 2*sigma1-quantile(sigma, prob = 0.025)
c(lowerEndpoint, upperEndpoint)


c(quantile(sigma, prob = 0.025), quantile(sigma, prob = 0.975))
