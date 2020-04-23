### DATA GENERATION
mu<-c(10,25)
sig<-c(5,5)
p<-c(.35,.65)
z<-sample(1:length(p),prob = p,replace = T,size = 200)
x<-rnorm(200,mu[z],sig[z])

### PRIOR INITIALISATIONS AND VECTOR/MATRIX INITIALISATIONS
p_init<-runif(1)
K_init<-rbinom(200,1,p_init)
mu_init<-rnorm(2,0,10)
sig_init<-sqrt(1/rgamma(1,0.01,0.01))

sig<-numeric(20000)
p<-numeric(20000)
mu<-matrix(,20000,2)

### ITERATIONS
for(i in 1:20001)
{
  K<-rbinom(200,1,p_init*dnorm(x,mu_init[1],sig_init)/
              (p_init*dnorm(x,mu_init[1],sig_init)+(1-p_init)*dnorm(x,mu_init[2],sig_init)))
  p[i]<-rbeta(sum(K==0)+1,sum(K==1)+1)
  sig[i]<-sqrt(1/rgamma(1,,))
  mu[i,1]<-rnorm(1,sum(K==0)*mean(x[K==0])*100/(100*sum(K==0)+sig[i]^2),100*sig[i]^2/(100*sum(K==0)+sig[i]^2))
  mu[i,2]<-rnorm(1,sum(K==1)*mean(x[K==1])*100/(100*sum(K==1)+sig[i]^2),100*sig[i]^2/(100*sum(K==1)+sig[i]^2))
  p_init<-p[i]
  sig_init<-sig[i]
  mu_init<-mu[i,]
}
### BURN-IN AND THINNING
index<-seq(from=2001,to=20001,by=10)
mean(p[index])
apply(mu[index,],2,mean)
mean(sig[index])