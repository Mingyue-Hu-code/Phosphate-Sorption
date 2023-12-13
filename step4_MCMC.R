freundpost<-function(pars,dat,priorpars,proposalpars,B,M){
  #overall Gibbs algorithm for extended freundlich model with AR1 errors
  #three Gibbs steps, one is Metropolis for alpha, beta
  #one is Metropolis for  gamma
  #other is Metropolis for theta, sigma
  #pars is vector of starting values alpha, beta, gamma, theta, sigma
  #dat has names $x and $y
  #priorpars is vector M,V,A,B,L,P,S
  #proposalpars is vector valp,vbet,vgam,vthet,vsig
  #B is burn-in
  #M is number of MC iterations to collect after burn-in
  ind1<-NULL; ind2<-NULL; ind3<-NULL
  alps<-NULL; bets<-NULL; gams<-NULL; thets<-NULL; sigs<-NULL
  cpars<-pars
  cnt<-0
  repeat{
    cnt<-cnt+1
    #metropolis for alpha, beta
    abstar<-abproposals(cpars,proposalpars)
    abres<-abmetrop(cpars,abstar,priorpars,proposalpars,dat)
    if(cnt>B) ind1<-c(ind1,abres[1])
    newpars1<-abres[-1]
    #metropolis for gamma
    gamstar<-gproposal(newpars1,proposalpars)
    gres<-gmetrop(newpars1,gamstar,priorpars,proposalpars,dat)
    if(cnt>B) ind2<-c(ind2,gres[1])
    newpars2<-gres[-1]
    #metropolis for theta, sigma
    tsstar<-tsproposals(newpars2,proposalpars)
    tsres<-tsmetrop(newpars2,tsstar,priorpars,proposalpars,dat)
    if(cnt>B) ind3<-c(ind3,tsres[1])
    newpars<-tsres[-1]
    if(cnt>B){
      alps<-c(alps,newpars[1])
      bets<-c(bets,newpars[2])
      gams<-c(gams,newpars[3])
      thets<-c(thets,newpars[4])
      sigs<-c(sigs,newpars[5])}
    cpars<-newpars
    if(cnt==B+M) break
  }
  res<-data.frame(ind1=ind1,ind2=ind2,ind3=ind3,alp=alps,bet=bets,gam=gams,thet=thets,sig=sigs)
  return(res)
}
#-------------------------------------------------------------------------
prioralp<-function(alp,M,V){
  #compute prior pdf for alpha
  #truncated normal (at zero)
  #
  num<-dnorm(alp,M,sqrt(V))
  den<-pnorm(alp,M,sqrt(V))
  palp<-num/den
  return(palp)
}
#------------------------------------------------
priorbet<-function(bet,A,B){
  #compute prior pdf for beta
  #gamma
  pbet<-dgamma(bet,A,B)
  return(pbet)
}
#------------------------------------------------
priorgam<-function(gam,L){
  #compute prior for gamma
  #exponential
  pgam<-dgamma(gam,1,L)
  return(pgam)
}
#---------------------------------------------
priorthet<-function(thet,P){
  #compute prior for theta
  #Laplace mean 0 and variance 2 phi^2
  #truncated at -1 and 1
  pthet<-(1/(2*P))*exp(-abs(thet)/P)
  pthet2<-function(thet,phi){(1/(2*phi))*exp(-abs(thet)/phi)}
  den<-integrate(pthet2,lower=-1,upper=1,phi=P)$value
  pthet3<-pthet/den
  return(pthet3)
}
#----------------------------------------------
priorsig<-function(sig,S){
  #compute prior for sigma
  #uniform (0,S)
  psig<-(1/S)*((sig>0) & (sig<S))
  return(psig)
}
#---------------------------------------------
datmod<-function(dat,pars){
  #compute joint density for observed responses
  #pars is vector alpha, beta, gamma, theta, sigma (not sigma2)
  #dat has names $x and $y
  #
  ys<-dat$y
  xs<-dat$x
  n<-length(ys)
  thet<-pars[4]; sig<-pars[5]
  mus<-freundlich(xs,pars)
  zs<-ys-mus
  zsm1<-c(0,zs[1:(n-1)])
  k<-(n/2)*log(2*pi*sig^2)
  eterm<-(0.5/sig^2)*sum((zs-thet*zsm1)^2)
  joint<-exp(-1*k-eterm)
  return(joint)
}
#----------------------------------------------------
freundlich<-function(xs,pars){
  #compute extended Freundlich response curve
  #pars is alpha, beta, gamma
  #xs are covariates (e.g., phosphorus concentration in soil)
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]
  epart<-bet*xs^(-gam)
  fs<-alp*xs^epart
  #res<-data.frame(x=xs,f=fs)
  res<-fs
  return(res)
}
#--------------------------------------------------------
abproposals<-function(pars,proposalpars){
  #proposal values for alpha, beta
  #proposalpars is vector valp, vbet, vgam, vthet, vsig
  #all from random walks
  #all restricted to be positive
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  valp<-proposalpars[1]; vbet<-proposalpars[2]
  zalp<-rnorm(1,0,sqrt(valp))
  zbet<-rnorm(1,0,sqrt(vbet))
  alpstar<-alp+zalp
  betstar<-bet+zbet
  if(alpstar<0) alpstar<-alp
  if(betstar<0) betstar<-bet
  res<-c(alpstar,betstar,gam,thet,sig)
  return(res)
}
#-------------------------------------------------------
gproposal<-function(pars,proposalpars){
  #proposal value for gamma
  #proposalpars is vector valp, vbet, vgam, vthet, vsig
  #from random walk
  #restricted to be positive
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  vgam<-proposalpars[3]
  zgam<-rnorm(1,0,sqrt(vgam))
  gamstar<-gam+zgam
  if(gamstar<0) gamstar<-gam
  res<-c(alp,bet,gamstar,thet,sig)
  return(res)
}
#-------------------------------------------------------
tsproposals<-function(pars,proposalpars){
  #proposal values for theta and sigma
  #both random walks
  #theta restricted to be on the interval (-1,1)
  #and sigma restricted to be positive 
  #pars is vector alpha, beta, gamma, theta, sigma
  #proposalpars is vector valp, vbet, vgam, vthet, vsig
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  vthet<-proposalpars[4]; vsig<-proposalpars[5]
  zthet<-rnorm(1,0,sqrt(vthet))
  zsig<-rnorm(1,0,sqrt(vsig))
  thetstar<-thet+zthet
  sigstar<-sig+zsig
  if((thetstar< -1) | (thetstar>1)) thetstar<-thet
  if(sigstar<0) sigstar<-sig
  res<-c(alp,bet,gam,thetstar,sigstar)
  return(res)
}
#------------------------------------------------------------
abmetrop<-function(pars,stars,priorpars,vars,dat){
  #conduct metropolis  for alpha, beta
  #pars is vector alpha, beta, gamma, theta, sigma -- current values
  #stars are jump proposals (alphastar, betastar, gamma, theta, sigma)
  #as produced by function abproposals
  #priorpars is vector M,V (for alpha), A,B (for beta), L (for gamma), PHI (for theta), S (for sigma)
  #vars are variances for proposal distributions
  #all computations are in log scale
  #then compbine into Metropolis acceptance probability for alpha, beta
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  alpstar<-stars[1]; betstar<-stars[2]
  valp<-vars[1]; vbet<-vars[2]
  M<-priorpars[1]; V<-priorpars[2]; A<-priorpars[3]; B<-priorpars[4]
  pialp<-log(prioralp(alp,M,V))
  pialpstar<-log(prioralp(alpstar,M,V))
  pibet<-log(priorbet(bet,A,B))
  pibetstar<-log(priorbet(betstar,A,B))
  fypar<-log(datmod(dat,pars))
  fystar<-log(datmod(dat,stars))
  metprob<-exp(pialpstar+pibetstar+fystar-pialp-pibet-fypar)
  accept<-min(c(1,metprob))
  ustar<-runif(1,0,1)
  if(ustar<=accept){ newpars<-stars; ind1<-1}
  if(ustar>accept){ newpars<-pars; ind1<-0}
  res<-c(ind1,newpars)
  return(res)
}
#-------------------------------------------------------
gmetrop<-function(pars,stars,priorpars,vars,dat){
  #conduct metropolis  for gamma
  #pars is vector alpha, beta, gamma, theta, sigma -- current values
  #stars are jump proposals (alpha, beta, gammastar, theta, sigma)
  #as produced by function gproposal
  #priorpars is vector M,V (for alpha), A,B (for beta), L (for gamma), PHI (for theta), S (for sigma)
  #vars are variances for proposal distributions
  #all computations are in log scale
  #then compbine into Metropolis acceptance probability for gamma
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  gamstar<-stars[3]
  vgam<-vars[3]
  L<-priorpars[5]
  pigam<-log(priorgam(gam,L))
  pigamstar<-log(priorgam(gamstar,L))
  fypar<-log(datmod(dat,pars))
  fystar<-log(datmod(dat,stars))
  qgamstar<-(-0.5/vgam)*(gamstar-gam)^2
  qgam<-(-0.5/vgam)*(gam-gamstar)^2
  metprob<-exp(pigamstar+fystar-pigam-fypar)
  accept<-min(c(1,metprob))
  ustar<-runif(1,0,1)
  if(ustar<=accept){ newpars<-stars; ind1<-1}
  if(ustar>accept){ newpars<-pars; ind1<-0}
  res<-c(ind1,newpars)
  return(res)
}
#--------------------------------------------------------------
tsmetrop<-function(pars,stars,priorpars,vars,dat){
  #compute proposal densities for theta and sigma
  #pars is vector alpha, beta, gamma, theta, sigma -- current values
  #stars are jump proposals (alpha, beta, gamma, thetastar, sigmastar)
  #with the star values produced by tsproposals
  #priorpars is vector M V (for alpha), A B (for beta), L (for gamma),  PHI (for theta), S (for sigma)
  #vars are variances for proposal distributions
  #all computations are in log scale
  #then compbine into Metropolis acceptance probability for theta and sigma
  #
  alp<-pars[1]; bet<-pars[2]; gam<-pars[3]; thet<-pars[4]; sig<-pars[5]
  thetstar<-stars[1]; sigstar<-stars[2]
  vthet<-vars[4]; vsig<-vars[5]
  PHI<-priorpars[6]; S<-priorpars[7]
  pithet<-log(priorthet(thet,PHI))
  pisig<-log(priorsig(sig,S))
  pithetstar<-log(priorthet(thetstar,PHI))
  pisigstar<-log(priorsig(sigstar,S))
  fypar<-log(datmod(dat,pars))
  fystar<-log(datmod(dat,stars))
  metprob<-exp(pithetstar+pisigstar+fystar-pithet-pisig-fypar)
  accept<-min(c(1,metprob))
  ustar<-runif(1,0,1)
  if(ustar<=accept){ newpars<-stars; ind2<-1}
  if(ustar>accept){ newpars<-pars; ind2<-0}
  res<-c(ind2,newpars)
  return(res)
}
#--------------------------------------------------------------
setwd("~/Desktop")
data <- read.table('projdat.txt',header = T)
b = 5000
b = 0
m <- 50000

pars_list <- list(
  c(0.7323, 1.0336, 0.1576, 0.5585, 0.6244),
  c(1.0793, 0.8042, 0.1444, 0.5762, 0.6437),
  c(0.85, 0.9, 0.15, 0.56, 0.6)
)

result_burning_list <- lapply(pars_list, function(pars) {
  freundpost(pars = pars,
             dat = data, priorpars = c(0.5,10,0.01,0.01,15,0.35,20),
             proposalpars = c(0.01,0.01,0.1,0.01,0.01), B = b, M = m)
})

library(patchwork)
library(tidyverse)

result_burning_list[[1]] %>% head
burning_plot_list <- lapply(result_burning_list, function(df) {
  df$index <- 1:nrow(df)
  p <- df %>% ggplot() +
    geom_line(aes(x = index, y = alp))
  p
})
alpha_plot <- burning_plot_list[[1]] / burning_plot_list[[2]] / burning_plot_list[[3]]
alpha_plot

# re-run with b <- 5000
burning_plot_list2 <- lapply(result_burning_list, function(df) {
  df$index <- 1:nrow(df)
  p <- df %>% ggplot() +
    geom_line(aes(x = index, y = alp))
  p
})
alpha_plot2 <- burning_plot_list2[[1]] / burning_plot_list2[[2]] / burning_plot_list2[[3]]
alpha_plot2

# question 2
b <- 5000

# c(0.7,1.02,0.15,0.6,0.8)

result <- freundpost(pars = pars_list[[1]],
                     dat = data, priorpars = c(0.5,10,0.01,0.01,15,0.35,20),
                     proposalpars = c(0.01,0.01,0.1,0.01,0.01),B = b, M = m)
cat(sum(result$ind1)/m,sum(result$ind2)/m,sum(result$ind3)/m)
plot(result$alp,type = "l", ylab = "alpha")
result1 <- result[result$ind1==1,]
# summary(result4$alp)
# summary(result4$bet)
result2 <- result[result$ind2==1,]
# summary(result2$gam)
result3 <- result[result$ind3==1,]
# summary(result3$thet)
# summary(result3$sig)
#credible interval

five_ns <- rbind(summary(result1$alp),
                 summary(result1$bet),
                 summary(result2$gam),
                 summary(result3$thet),
                 summary(result3$sig))
rownames(five_ns) <- c('alpha', 'beta', 'gamma', 'theta', 'sigma')
five_ns %>% round(4)

est_list <- list(result1, result2, result3)
lower_upper_list <- lapply(est_list, function(est) {
  lower <- apply(est, 2, quantile, 0.025)
  upper <- apply(est, 2, quantile, 0.975)
  
  return(data.frame(lower = lower, upper = upper))
})
lower_upper_result <- rbind(
  lower_upper_list[[1]][4:5, ],
  lower_upper_list[[2]][6, ],
  lower_upper_list[[3]][7:8, ]
)


# est<-result2[,4:8]
# lower <- apply(est, 2, quantile, 0.025);lower
# upper <- apply(est, 2, quantile, 0.975);upper
result4 <- result[result$ind1==1 &result$ind2==1 &result$ind3==1,]
est1<-result4[,4:8]
cor(est1)
xtable(cor(est1))

#########################
# get tables
xtable(five_ns, digits = 4)
xtable(lower_upper_result, digits = 4)
xtable(cor(est1), digits = 4)

########################

median(est)
library(tidyverse)

fctn<-function(xmat, cps){
  cps[1]*xmat^(cps[2]*xmat^(-cps[3]))
}

result_tb <- as_tibble(result)
#result_tb <- result_tb[(b+1):nrow(result_tb), ]
result_tb <- result_tb[result_tb$ind1==1 &result_tb$ind2==1 &result_tb$ind3==1,]

result_tb <- result_tb %>% mutate(
  response = purrr::pmap(list(alp, bet, gam), function(aa, bb, gg, xmat) {
    fctn(xmat, c(aa, bb, gg)  )
  }, xmat = data$x)
)

response_mat <- do.call(cbind, result_tb$response) %>% t()
response_mat_after_b <- response_mat
#response_mat_after_b <- response_mat[(b+1):nrow(response_mat), ]

response_median <- apply(response_mat_after_b, 2, median)
response_lower <- apply(response_mat_after_b, 2, quantile, 0.025)
response_upper <- apply(response_mat_after_b, 2, quantile, 0.975)


plot(data$x, response_median, type = 'l', ylim = c(min(response_lower), max(response_upper)),
     xlab = 'x', ylab = 'y')
lines(data$x, response_lower)
lines(data$x, response_upper)
points(data$x, data$y)

result <- freundpost(pars = c(0.7,1.02,0.15,0.6,0.8),
                     dat = data, priorpars = c(0.5,10,0.01,0.01,15,0.35,20),
                     proposalpars = c(0.01,0.01,0.1,0.01,0.01),B = b, M = m)
result <- freundpost(pars = pars_list[[1]],
                     dat = data, priorpars = c(0.5,10,0.01,0.01,15,0.35,20),
                     proposalpars = c(0.01,0.01,0.1,0.01,0.01),B = b, M = m)

###################################################


#compare
data <- read.table('projdat.txt',header = T)
b = 5000
m <- 50000
result <- freundpost(pars = c(0.7,1.02,0.15,0.6,0.8),
                     dat = data, priorpars = c(0.5,10,0.01,0.01,15,0.35,20),
                     proposalpars = c(0.01,0.01,0.1,0.01,0.01),B = b, M = m)
library(xtable)

result_frame <- data.frame(lower = lower, 
                           upper = upper)
xtable(result_frame)






