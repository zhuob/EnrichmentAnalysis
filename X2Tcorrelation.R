#  General conclusion:  
#  1. This simulation shows that, generally, correlation of T statistics
#     are not the same as sample correlation after accounting for treatment effects.


library(MASS)

mu1 <- c(0, 8); mu2 <- c(8,0)
sig1 <- 2; sig2 <- 1.5; rho <- 0.5
sigma <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
size <- 100
nsim <- 100
t <- matrix(NA, ncol=2 ,nrow=nsim)
rho <- c()
for ( i in 1: nsim)
{
  z1 <- mvrnorm(size, mu1, sigma);  z2 <- mvrnorm(size, mu2, sigma)
  z <- rbind(z1, z2)
  t[i, 1] <- t.test(z1[, 1], z2[, 1])$stat; 
  t[i, 2] <- t.test(z1[, 2], z2[, 2])$stat
  
  rm1 <- z1[, 1] - mean(z1[, 1])
  rm2 <- z1[, 2] - mean(z1[, 2])
  
  rm3 <- z2[, 1] - mean(z2[, 1])
  rm4 <- z2[, 2] - mean(z2[, 2])
  rho[i] <- cor(c(rm1, rm3), c(rm2, rm4))

}



## t to z
# the conclusion here is that, if one of the test is non-null, then
# the correlation of t statistic is not consistent

nsim <- 100
d <- 20
df <- 20
x <- rt(nsim, df)
z <- qnorm(pt(x + d, df= df))
cor(x, z)

N <- 1000
nsim <- 20
stat.t <- matrix(NA, N, 3)
rho <- matrix(NA, N, 2)

for ( i in 1:N)
{ 
  sig <- 10
  sigma <- sig *matrix(c(1, 0.5, 0.5, 1), 2)
  set.seed(i)
  x1 <- mvrnorm(nsim, mu = c(20, 30), sigma)
  set.seed(1e4 + i+5)
  x2 <- mvrnorm(nsim, mu= c(0, 32), sigma)
  s1 <- apply(x1, 2, var)
  s2 <- apply(x2, 2, var)
  mu1 <- apply(x1, 2, mean)
  mu2 <- apply(x2, 2, mean)
  denom1 <-  (sqrt(s1[1]/nsim + s2[1]/nsim))
  denom2 <- (sqrt(s1[2]/nsim + s2[2]/nsim))
  
  # denom1 <- 1
  #    denom2 <- 10
  
  t1 <- (mu1[1]-mu2[1]-18)/denom1
  t2 <- (mu1[2]-mu2[2])/denom2
  t3 <- (mu1[1]-mu2[1])/denom1
  
  #   t1 <- t.test(x1[, 1]-18, x2[, 1])$stat
  #   t2 <- t.test(x1[, 2], x2[, 2])$stat
  
  stat.t[i,] <- c(t1, t2, t3)
  rho[i, ] <- c(cor(x1)[1, 2], cor(x2)[1, 2])
}

apply(rho,2, mean)
cor(stat.t)  #col 1-3: Null t, null t, non-null t





## plot z test correlation, sample.cor and true cor

plot.cor.z.test <- function(mu1, mu2, sig1, sig2, sig3, sig4, sigqua){
  
  
  true.cor  <- seq(-0.9, 0.9, 0.01) 
  z.cor <- samp.cor <- c()
  
  for ( j in 1: length(true.cor))
  {
    permu.B <- 100
    rho <- true.cor[j]
    
    # mean 
    #mu1 <- c(20, 30)
    #mu2 <- c(0, 32)
    
    # covariance 
    sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
    sigma2 <- matrix(c(sig3^2, rho*sig3*sig4, rho*sig3*sig4, sig4^2), 2, 2)
    nsim1 <- 20
    nsim2 <- 30
    
    z.stat <- matrix(NA, permu.B, 2)
    s.cor <- c()
    for ( i in 1:permu.B)
    {
      
      x1 <- mvrnorm(nsim1, mu = mu1, sigma1)
      x2 <- mvrnorm(nsim2, mu = mu2, sigma2)
      
      s1 <- sqrt(sigma1[1,1]/nsim1 + sigma2[1, 1]/nsim2)
      s2 <- sqrt(sigma1[2,2]/nsim1 + sigma2[2, 2]/nsim2)
      
      t1 <- (mean(x1[, 1]) - mean(x2[, 1]))/s1
      t2 <- (mean(x1[, 2]) - mean(x2[, 2]))/s2
      
      
      
      rm1 <- x1[, 1] - mean(x1[, 1])
      rm2 <- x1[, 2] - mean(x1[, 2])
      
      rm3 <- x2[, 1] - mean(x2[, 1])
      rm4 <- x2[, 2] - mean(x2[, 2])
      s.cor[i] <- cor(c(rm1, rm3), c(rm2, rm4)) # sample corr

      z.stat[i, ] <- c(t1, t2)
      
    }
    
    z.cor[j] <- cor(z.stat)[1, 2]
    samp.cor[j] <- mean(s.cor)
    
    
  }
  
  plot(true.cor, samp.cor, pch=3, cex=0.4, main=sigqua)
  abline(0, 1)
  points(true.cor, z.cor, pch=20, cex=0.5, col='red')
  legend("topleft", legend=c("sample", "test stat"), pch=c(3, 20), col=c("black", "red"))
}

mu1 <- c(20, 60)
mu2 <- c(5, 30)

plot.cor.z.test(mu1, mu2, sig1=11, sig2= 3,sig3=0.2, sig4=0.2, sigqua="DE")

mu1 <- c(20, 60)
mu2 <- c(200, 120)
plot.cor.z.test(mu1, mu2, sig1=100, sig2= 5,sig3=100, sig4=5, sigqua="rho constant")
plot.cor.z.test(mu1, mu2, sig1=10, sig2= 10, sig3=10, sig4=10, sigqua="Null") 

#####################################################################
# for ,  sig1*sig4 must equal to sig2*sig3, which is often true because 
# we assume for two genes, sig1=sig3 and sig2=sig4
#####################################################################
plot.cor.z.test(mu1, mu2, sig1=0.2, sig2= 0.3,sig3=0.2, sig4=0.3, sigqua="DE")
plot.cor.z.test(mu1, mu2, sig1=50, sig2= 0.3,sig3=0.03, sig4=10, sigqua="DE")





### use pooled correlation? ??
## be cautious 

## THIS PART SHOWS THAT YOU CANNOT COMBINE TWO DATA SETS IF THEY DON'T HAVE 
## HOMOGENEOUS VARIANCE
#
## THIS IS THE CASE WHERE TWO LAB EXPERIMENTS ARE NOT COMPARABLE ( i.e., in this
# case, we cannot combine two data sets to calculate sample correlation)

rho <- 0.3
sig1=50
sig2= 0.1*sig1
mu1 <- mu2 <-  c(0, 0)
sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
rand1 <- mvrnorm(1000, mu1, sigma1)
cor(rand1)

sig3=50
sig4=.1*sig3

sigma2 <- matrix(c(sig3^2, rho*sig3*sig4, rho*sig3*sig4, sig4^2), 2, 2)
rand2 <- mvrnorm(1000, mu2, sigma2)
cor(rand2)
rand <- rbind(rand1, rand2)
cor(rand)
tt1 <- c(rand1[, 1], rand2[, 1])
tt2 <- c(rand1[, 2], rand2[, 2])
cor(tt1, tt2)


#####################################################################################
# READ FROM HERE
#####################################################################################



## explore the correlations of two sample t test 

#   sig: sigma^2
# coeff: true correlation
#   mu1: mean vector for lab 1
#   mu2: mean vector for lab 2

preseve.cor <- function(mu1, mu2, sig, rho,  nsim) 
  # coeff is the correlation between two genes
{
  # covariance matrix
  # sigma <- sig*matrix(c(1, coeff, coeff, 1), 2, 2)
    sigma <- sig

  t.stat1 <- matrix(NA, nsim, 2) # store the test stats
  # rho1 <- t.stat1 # store the sample corr
  rho1 <- c()
  num1 <- 20  # number of sample size 
  
  # study correlation 
  for ( i in 1:nsim)
  {
    z1 <- mvrnorm(num1, mu=mu1, sigma)
    #  z1 contains obs from two correlated  genes,  rho=0.6
    z2 <- mvrnorm(num1, mu=mu2, sigma)
    # z2  is just another replicate of the same genes
    #  gene1 Not DE, gene2 DE
    
    # two sample t.test
    t1 <- t.test(z1[, 1], z2[, 1])$stat
    t2 <- t.test(z1[, 2], z2[, 2])$stat
    
    # rho1[i, ] <- c(cor(z1)[1, 2], cor(z2)[1, 2])
    rm1 <- z1[, 1] - mean(z1[, 1])
    rm2 <- z1[, 2] - mean(z1[, 2])
    
    rm3 <- z2[, 1] - mean(z2[, 1])
    rm4 <- z2[, 2] - mean(z2[, 2])
    rho1[i] <- cor(c(rm1, rm3), c(rm2, rm4))
    
    t.stat1[i, ] <- c(t1, t2)
    
  }
  
  t.cor <- cor(t.stat1)[1, 2]# correlation between t test stats
  #p.t <- pt(t.stat1, 2*num1-2) # p value does not preserve the correlation
  #p.cor <- cor(p.t)[1, 2] # correlation between p values
  # mean.coeff <- apply(rho1, 2, mean)
  mean.coeff <- mean(rho1)
  # z <- qnorm(p.t) # z value is similar to test stat and thus sample correlation
  # z.cor <- cor(z)[1, 2] # correlation between z test stats
  
  return(c(t.cor, mean.coeff, rho))
}



## the correlation cannot be preserved if the null are not true
## even if we consider the test statistics

#################################################################
# HOWEVER, IN USUAL ANALYSIS, WE WILL REMOVE THE GROUP MEANS FIRST
# AND IN THAT WAY, THE T STAT CORRELATION IS ABOUT THE SAME AS SAMPLE CORRELATION
#################################################################
coeff <- seq(-0.9, 0.9, by = 0.01)
mu1 <- c(0, 0)
mu2 <- c(0, 0)
sig <- 10  # if the variance is too large, the correlation preserves 


corr.matrix <- data.frame(matrix(NA, length(coeff), 3))
colnames(corr.matrix) <- c("t","rho", "rho.true")

for (i in 1:length(coeff))
{
  sig1 <- 2; sig2 <- 2; rho <- coeff[i]
  # print(rho)
    sig <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
  corr.matrix[i, ] <- preseve.cor(mu1, mu2, sig,rho, nsim = 100)
  
}

plot(corr.matrix$rho.true, corr.matrix$rho, pch=20, cex=0.3, 
     xlab="rho.true", ylab="rho",
     main = "NON-DE, equal variance")
abline(0, 1)
points(corr.matrix$rho.true, corr.matrix$t, cex= 0.5, pch =3, col="magenta")




coeff <- seq(-0.9, 0.9, by = 0.01)
mu1 <- c(0, 0)
mu2 <- c(0, 0)
sig <- 10  # if the variance is too large, the correlation preserves 


corr.matrix <- data.frame(matrix(NA, length(coeff), 3))
colnames(corr.matrix) <- c("t","rho", "rho.true")

for (i in 1:length(coeff))
{
  sig1 <- 2; sig2 <- 10; rho <- coeff[i]
  # print(rho)
  sig <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
  corr.matrix[i, ] <- preseve.cor(mu1, mu2, sig,rho, nsim = 100)
  
}

plot(corr.matrix$rho.true, corr.matrix$rho, pch=20, cex=0.3, 
     xlab="rho.true", ylab="rho",
     main = "NON-DE, UNequal variance")
abline(0, 1)
points(corr.matrix$rho.true, corr.matrix$t, cex= 0.5, pch =3, col="magenta")




coeff <- seq(-0.9, 0.9, by = 0.01)
mu1 <- c(8, 0)
mu2 <- c(0, 0)
sig <- 10  # if the variance is too large, the correlation preserves 


corr.matrix <- data.frame(matrix(NA, length(coeff), 3))
colnames(corr.matrix) <- c("t","rho", "rho.true")

for (i in 1:length(coeff))
{
  sig1 <- 2; sig2 <- 2; rho <- coeff[i]
  # print(rho)
  sig <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
  corr.matrix[i, ] <- preseve.cor(mu1, mu2, sig,rho, nsim = 100)
  
}

plot(corr.matrix$rho.true, corr.matrix$rho, pch=20, cex=0.3, 
     xlab="rho.true", ylab="rho",
     main = "DE, equal variance")
abline(0, 1)
points(corr.matrix$rho.true, corr.matrix$t, cex= 0.5, pch =3, col="magenta")


coeff <- seq(-0.9, 0.9, by = 0.01)
mu1 <- c(8, 0)
mu2 <- c(0, 0)
sig <- 10  # if the variance is too large, the correlation preserves 


corr.matrix <- data.frame(matrix(NA, length(coeff), 3))
colnames(corr.matrix) <- c("t","rho", "rho.true")

for (i in 1:length(coeff))
{
  sig1 <- 2; sig2 <- 10; rho <- coeff[i]
  # print(rho)
  sig <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
  corr.matrix[i, ] <- preseve.cor(mu1, mu2, sig,rho, nsim = 100)
  
}

plot(corr.matrix$rho.true, corr.matrix$rho, pch=20, cex=0.3, 
     xlab="rho.true", ylab="rho",
     main = "DE, UNequal variance")
abline(0, 1)
points(corr.matrix$rho.true, corr.matrix$t, cex= 0.5, pch =3, col="magenta")



## The correlation of test stat decreases as DE becomes more significant.



rho <- c()


for ( j in 1:100)
{
mu1 <- c(0, 0)
mu2 <- c(0, 50)  ### change the magnitude of DE values
sigma <- matrix(c(4, 3, 3, 9), 2, 2)

trt1 <- mvrnorm(n=1000, mu1, Sigma =sigma )
trt2 <- mvrnorm(n=1000, mu2, Sigma =sigma )

t1 <- t2 <- c()
for ( i in 1:100)
{
  ids <- sample(1:1000, 100)
  trt1.samp <- trt1[ids, ] # get 100 samples from 
  trt2.samp <- trt2[ids, ]
  
  t1[i] <- t.test(trt1.samp[, 1], trt2.samp[, 1])$stat
  t2[i] <- t.test(trt1.samp[, 2], trt2.samp[, 2])$stat
    
}
rho[j] <- cor(t1, t2)
}

mean(rho)













## plot t.cor, sample.cor and true cor
plot.cor.t.test <- function(mu1, mu2, sig1, sig2, sig3, sig4, sigqua){
  
  true.cor  <- seq(-0.9, 0.9, 0.01) 
  t.cor <- samp.cor <- c()
  
  for ( j in 1: length(true.cor))
  {
    permu.B <- 100
    rho <- true.cor[j]
    
    # covariance 
    sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
    sigma2 <- matrix(c(sig3^2, rho*sig3*sig4, rho*sig3*sig4, sig4^2), 2, 2)
    nsim1 <- 20
    nsim2 <- 30
    
    t.stat <- matrix(NA, permu.B, 2)
    s.cor <- c()
    for ( i in 1:permu.B)
    {
      
      x1 <- mvrnorm(nsim1, mu = mu1, sigma1)
      x2 <- mvrnorm(nsim2, mu = mu2, sigma2)
      
      t1 <- t.test(x1[, 1], x2[, 1])$stat
      t2 <- t.test(x1[, 2], x2[, 2])$stat
      
      ## I use mean of sample correlation as the sample correlation.
      rho1 <- cor(x1)[1,2]
      rho2 <- cor(x2)[1,2]
      s.cor[i] <- mean(c(rho1, rho2))
      # s.cor[i] <- cor(c(rm1, rm3), c(rm2, rm4)) # sample corr
      
      t.stat[i, ] <- c(t1, t2)
      
    }
    
    t.cor[j] <- cor(t.stat)[1, 2]
    samp.cor[j] <- mean(s.cor)
    
    
  }
  
  plot(true.cor, samp.cor, pch=3, cex=0.4, main=sigqua)
  abline(0, 1)
  points(true.cor, t.cor, pch=20, cex=0.5, col='red')
  legend("topleft", legend=c("sample", "test stat"), pch=c(3, 20), col=c("black", "red"))
}

mu1 <- c(20, 60)
mu2 <- c(5, 30)

plot.cor.t.test(mu1, mu2, sig1=11, sig2= 3,sig3=0.2, sig4=0.2, sigqua="DE")

mu1 <- c(20, 60)
mu2 <- c(200, 120)
plot.cor.z.test(mu1, mu2, sig1=100, sig2= 5,sig3=10, sig4=0.5, sigqua="rho constant")
plot.cor.z.test(mu1, mu2, sig1=10, sig2= 10, sig3=10, sig4=10, sigqua="Null") 
# for constant rho,  sig1*sig4 must equal to sig2*sig3
plot.cor.z.test(mu1, mu2, sig1=0.2, sig2= 0.3,sig3=0.2, sig4=0.3, sigqua="DE")
plot.cor.z.test(mu1, mu2, sig1=50, sig2= 0.3,sig3=0.03, sig4=10, sigqua="DE")




