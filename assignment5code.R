### DATA GENERATION
mu<-c(10,25)
sig<-c(5,5)
p<-c(.35,.65)
z<-sample(1:length(p),prob = p,replace = T,size = 200)
x<-rnorm(200,mu[z],sig[z])

### PRIOR INITIALISATIONS AND VECTOR/MATRIX INITIALISATIONS
p_init<-runif(1)
K_init<-rbinom(200,1,p_init)
mu_init<-rnorm(2)
sig_init<-sqrt(1/rgamma(1,0.01,0.01))



### ITERATIONS
for(i in 1:20001)
{
  K<-rbinom(200,1,p_init*dnorm(x,mu_init[1],sig_init)/
              (p_init*dnorm(x,mu_init[1],sig_init)+(1-p_init)*dnorm(x,mu_init[2],sig_init)))
  p<-rbeta(sum(K==0)+1,sum(K==1)+1)
  sig<-sqrt(1/rgamma())
  mu[1]<-rnorm()
  mu[2]<-rnorm()
  p_init<-p
  sig_init<-sig
  mu_init<-mu
}
index<-seq(from=2001,to=20001,by=10)
mean(p[index])
apply(mu[index,],2,mean)
mean(sig[index])