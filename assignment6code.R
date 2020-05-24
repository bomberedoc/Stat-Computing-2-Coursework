#install.packages(c("LaplacesDemon","invgamma"))
library(LaplacesDemon)
library(invgamma)

TIME_=Sys.time()

generate_data<-function(n,mu1=10,mu2=25,sig=5,p=0.35)
{
  y<-rnorm(n,ifelse(runif(n)<p,mu1,mu2),sig)
  return(y)
}
count_k<-function(k,z)
{
  count<-numeric(k)
  for(i in 1:k)
    count[i]<-sum(z==i)
  return(count)
}
get_prob<-function(y,p,sig,mu)
{
  spl<-numeric(length(y))
  for(i in 1:length(y))
  {
    prob<-p*(dnorm(y[i],mu,sig))
    prob<-prob/sum(prob)
    spl[i]<-sample(1:length(p),1,prob = prob)
  }
  return(spl)
}
simulate_mu<-function(y,z,sig,k)
{
  mu<-numeric(k)
  for(i in 1:k)
  {
    mu[i]<-rnorm(1,sum(y[z==i])*100/(sig^2+100*sum(z==i)),sqrt(100*sig^2/(sig^2+100*sum(z==i))))
  }
  return(mu)
}
find_max<-function(z)
{
  return(which(table(z)==max(table(z)))[1])
}
ret_bic<-function(n,k,y)
{
  ### VECTOR/MATRIX INITIALISATIONS
  mu<-matrix(, nrow = 20000, ncol = k)
  sig<-numeric(20000)
  p<-matrix(,nrow = 20000,ncol = k)
  z<-matrix(,nrow = 20000,ncol = n)
  ### PRIORS
  mu_init<-rnorm(k,0,10)
  sig_init<-sqrt(1/rgamma(1,0.01,0.01))
  p_init<-rdirichlet(1,rep(1,k))
  z_init<-sample(1:k,n,replace = T,prob = p_init)
  ### ITERATIONS
  for(i in 1:20000)
  {
    mu[i,]<-simulate_mu(y=y,z=z_init,sig=sig_init,k=k)
    sig[i]<-sqrt(1/rgamma(1,n/2+0.01,0.01+sum((y-mu[i,z_init])^2)/2))
    p[i,]<-rdirichlet(1,count_k(k=k,z=z_init)+1)
    z[i,]<-get_prob(y=y,p=p[i,],sig=sig[i],mu=mu[i,])
    
    mu_init<-mu[i,]
    z_init<-z[i,]
    p_init<-p[i,]
    sig_init<-sig[i]
  }
  ### BURN-IN AND THINNING
  index<-seq(from=2001,to=20000,by=10)
  mu1<-apply(mu[index,], 2, mean)
  p1<-apply(p[index,], 2, mean)
  sig1<-mean(sig[index])
  z1<-apply(z[index,],2,find_max)
  ### LOG-LIKELIHOOD AND BIC CALCULATION
  lik<-sum(dnorm(y,mu1[z1],sig1,log = T))+sum(dnorm(mu1,0,10,log = T))+
    dinvgamma(sig1^2,0.01,0.01,log = T)+sum(count_k(k=k,z=z1)*log(p1))
  BIC<-k*log(n)-2*lik
  
  return(list(mu=mu[index,],p=p[index,],sigma=sig[index],z=z,BIC=BIC))
}
### DATA GENERATION
n<-200
k<-2:5
y<-generate_data(n=n)
### FINDING OPTIMUM VALUE OF K
a<-numeric(length(k))
for(i in 1:length(k))
  a[i]<-ret_bic(n=n,k=k[i],y=y)$BIC
k_opt<-k[which(a==min(a))]
### USING THE OPTIMUM VALUE OF K TO FIND ESTIMATES OF PARAMETERS OF OPTIMUM MODEL
output<-ret_bic(n=n,k=k_opt,y=y)

apply(output$p,2,mean)
apply(output$mu, 2, mean)
mean(output$sigma)
apply(output$p,2,quantile,c(0.025,0.975))
apply(output$mu,2,quantile,c(0.025,0.975))
quantile(output$sigma,c(0.025,0.975))

Sys.time()-TIME_