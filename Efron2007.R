library(MASS)



##   simple linear regression  correlation of estimate beta
nsim <- 20

x <- cbind( rep(1, 2*nsim), c(rep(0, nsim ), rep(1, nsim)))

beta <- c(10, 2) # gene1  control = 10, trt= 12
mu1 <- t(unique(x%*%beta))

delta <- c(10, 5) # gene2  control = 10 , trt = 15
mu2 <-  t(unique(x%*%delta))

x.inv <- solve(t(x)%*%x)






# generate correlated poisson random variable

cor.pois <- function(n, mu1, mu2, mu3)
{
  y1 <- rpois(n, mu1)
  y2 <- rpois(n, mu2)
  y3 <- rpois(n, mu3)
  
  x1 <- y1  + y2
  x2 <- y1  + y3
  
 # the correlation between x1 and x2 
 # rho <- mu1 /sqrt( (mu1 + mu3)*(mu1 + mu2))
#  print(round(rho, 4))
  return(cbind(x1, x2))
}



t.stat <- matrix(NA, 10000, 3)
z= rho = matrix(NA, 10000, 2)
num <- 20  # how many samples in each group 


for (i in 1:10000)
{
  # in this case, gene 1 is not DE, but gene 2 is DE
  # for treatment 1
  exp1 <- cor.pois(num, 10, 30, 40) # pi1 = 4, p2= 5, rho = 0.2236
  # for treatment 2
  exp2 <- cor.pois(num, 20, 20, 180) #pi1 = 4, p2= 20 rho = 0.2236
  
  t1 <- as.numeric(t.test(exp1[, 1], exp2[, 1])$stat)
  t2 <- as.numeric(t.test(exp1[, 2], exp2[, 2])$stat, alternative="less")
  
  t3 <- as.numeric(t.test(exp1[, 2], exp2[, 2]-150)$stat)
  
  
  z[i,] = colMeans(exp2) - colMeans(exp1)
  
  rho[i, ] <- c(cor(exp1)[1, 2], cor(exp2)[1, 2])
  t.stat[i, ] <- c(t1, t2, t3)
}

cor(t.stat) # correlation of t.stat is the same as cor(zval)
apply(rho, 2, mean)
pval <- pt(t.stat, 2*num-2) 
cor(pval)    #  but not for cor(pval)
zval <- qnorm(pval)
cor(zval)


# 
mu3 <- 3
mu4 <- 4
mu5 <- 6
mu1 <- 0.5*mu4
mu2 <- 0.5*mu4 + mu5
mu6 <- mu4 + 4*mu3

c(mu1, mu2, mu3, mu4, mu5, mu6)
rho1 <- mu1 /sqrt( (mu1 + mu3)*(mu1 + mu2))
rho2 <- mu4 /sqrt( (mu4 + mu5)*(mu4 + mu6))
c(rho1, rho2)



# How close is a t distribution to normal 
df <- 10
va <- df/(df-2)
pi <- seq(-2, 2, length=1000)

f1 <- dnorm(pi, 0, sqrt(va))
f2 <- dt(pi, df)

plot(pi, f1, pch=20, cex=0.3)
lines(pi, f2, pch=3, cex=0.3, col="red")


## how is 77 obtained in page 100 of Efron's paper
N <- 3226
A <- 0.57
x0 <- 2.5
bar.phi <- 1-pnorm(x0)
phi <- dnorm(x0)
cond <- N*bar.phi*(1 + A*x0*phi/(sqrt(2)*bar.phi))
cond

## reproduce figure 2 #####

FDP <- c()
A <- c()
A.hat <-c()

for ( i in 1:1000){
  
N <- 3000
breaks <- 50

z1 <- rnorm(0.95*N, 0, 1)
z2 <- rnorm(0.05*N, 2.5, sqrt(1.25))

z <- c(z1, z2)

out <- cut(z, quantile(z, seq(0, 1, len = 11)), include.lowest = TRUE) 

bins <- seq(min(z),max(z), length=breaks)
bin.center <- 0.5*(bins[-1] + bins[-breaks])
delta <- diff(bins)[1]

# define the phi funciton and bar.phi

v <- N*delta*dnorm(bin.center)

w <- N*delta*dnorm(bin.center)*(bin.center^2-1)/sqrt(2)
A[i] <- rnorm(1, 0, 0.15)

id1 <- which(abs(bins-(-1))==min(abs(bins-(-1)))) 
id2 <- which(abs(bins-(1))==min(abs(bins-(1)))) 
u <- v + A[i]*w
u[u<0] <- 1e-3
y.null <- rpois(length(u), u)

P0 <- 2*pnorm(1)-1
Q0 <- sqrt(2)*x0*dnorm(1)

Y0 <- sum(y.null[id1:id2])
#Y0 = sum(z1>-1 & z1<1)
P0.hat <- Y0/N
A.hat[i] <- (P0-P0.hat)/Q0

E.A <- N*(1-pnorm(x0)) *(1 + A.hat[i] *x0 *dnorm(x0)/(sqrt(2)*(1-pnorm(x0))))

T.x <- sum(z>x0)
FDP[i] <- E.A/T.x

}

plot(A, FDP, cex=0.5)


hist(z, breaks=30, freq=F)
lines(seq(-4, 4, by=0.001), dnorm(seq(-4, 4, 0.001)))

# the z transformation will preserve the correlation structure
z2 <- mvrnorm(100, mu=c(2, 0), sigma)
p <- pnorm(z2)
invz <- qnorm(p)
cor(z2)
cor(invz)



## exploring t.test: 
## the default t.test is two sided
nsim <- 20
x1 <- rnorm(nsim, 0, 1)
x2 <- rnorm(nsim, 1.2, 0.8)
t.test(x1, x2)$p.va

# by hand
s1 <- var(x1)
s2 <- var(x2)
t <- (mean(x1)-mean(x2))/ (sqrt(s1/nsim + s2/nsim))
df <- (s1/nsim + s2/nsim)^2/( 
  (s1/nsim)^2/(nsim-1) + (s2/nsim)^2/(nsim-1))
pt(t, df)*2 



## t to z
# the conclusion here is that, if one of the test is non-null, then
# the correlation is not consistent

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
cor(stat.t)





## plot z.cor, sample.cor and true cor
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
      #s.cor[i] <- cor(c(x1[, 1], x2[, 1]), c(x1[, 2], x2[, 2]))
#       p1 <- sum(x1[ ,1]*x1[ ,2])- nsim1*mean(x1[, 1])*mean(x1[, 2])
#       p2 <- sum(x2[ ,1]*x2[ ,2])- nsim2*mean(x2[, 1])*mean(x2[, 2])
#       d1 <- sqrt(sum(x1[, 1]^2) + sum(x2[, 1]^2)- nsim1*mean(x1[, 1])^2 - nsim2*mean(x2[, 1])^2)
#       d2 <- sqrt(sum(x1[, 2]^2) + sum(x2[, 2]^2)- nsim1*mean(x1[, 2])^2 - nsim2*mean(x2[, 2])^2)
#       r <- (p1 + p2)/(d1*d2)
        
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
plot.cor.z.test(mu1, mu2, sig1=100, sig2= 5,sig3=10, sig4=0.5, sigqua="rho constant")
plot.cor.z.test(mu1, mu2, sig1=10, sig2= 10, sig3=10, sig4=10, sigqua="Null") 
# for constant rho,  sig1*sig4 must equal to sig2*sig3
plot.cor.z.test(mu1, mu2, sig1=0.2, sig2= 0.3,sig3=0.2, sig4=0.3, sigqua="DE")
plot.cor.z.test(mu1, mu2, sig1=50, sig2= 0.3,sig3=0.03, sig4=10, sigqua="DE")



### use pooled correlation? ??
## be cautious 

rho <- 0.3
sig1=50
sig2= 0.1*sig1
mu1 <- mu2 <-  c(0, 0)
sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
rand1 <- mvrnorm(1000, mu1, sigma1)
cor(rand1)

sig4=50
sig3=0.1*sig4

sigma2 <- matrix(c(sig3^2, rho*sig3*sig4, rho*sig3*sig4, sig4^2), 2, 2)
rand2 <- mvrnorm(1000, mu2, sigma2)
cor(rand2)
rand <- rbind(rand1, rand2)
cor(rand)
tt1 <- c(rand1[, 1], rand2[, 1])
tt2 <- c(rand1[, 2], rand2[, 2])
cor(tt1, tt2)





## explore the correlations of two sample t test 
preseve.cor <- function(mu1, mu2, coeff, sig, nsim)  # coeff is the correlation between two genes
{
  
  sigma <- sig*matrix(c(1, coeff, coeff, 1), 2, 2)
  
  #
  t.stat1 <- matrix(NA, nsim, 2) # store the test stats
  # rho1 <- t.stat1 # store the sample corr
  rho1 <- c()
  num1 <- 20
  
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
  
  return(c(t.cor, mean.coeff, coeff))
}



## the correlation cannot be preserved if the null are not true
## even if we consider the test statistics
coeff <- seq(-0.9, 0.9, by = 0.01)
mu1 <- c(20, 30)
mu2 <- c(10, 40)
sig <- 10  # if the variance is too large, the correlation preserves 
corr.matrix <- data.frame(matrix(NA, length(coeff), 3))
colnames(corr.matrix) <- c("t","rho", "rho.true")

for (i in 1:length(coeff))
{
  corr.matrix[i, ] <- preseve.cor(mu1, mu2, coeff[i], sig, nsim = 100)
  
}

plot(corr.matrix$rho.true, corr.matrix$rho, pch=20, cex=0.3, 
     xlab="rho.true", ylab="rho",
     main = paste("mu1 = ", mu1, "mu2 = ", mu2, "sigma =", sig))
abline(0, 1)
points(corr.matrix$rho.true, corr.matrix$t, cex= 0.5, pch =3, col="magenta")


s1 <- matrix(NA, 10000, 2)

for ( i in 1: 10000)
{
  
  rho <- 0.6
  sig1=50
  sig2= 0.8*sig1
  mu1 <- mu2 <-  c(0, 0)
  sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
  rand1 <- mvrnorm(1000, mu1, sigma1)
  v <- apply(rand1, 2, var)
  s0 <- v[1]/sig1^2
  s2 <- v[2]/sig2^2
  s1[i, ] <- c(s0, s2)
}
cor(s1)




#########
### my point: null/nonnull affects the location of mean(z values)
### while correlation affects the variance of z values

##  if the correlation increases, the variance of z decreases.
library(MASS)
Cov.1 <- c()

rho <- seq(0.1, 0.99, by = 0.01)
for ( k in 1:length(rho))
{
  ndim <- 100
  mu <- rep(0, ndim)
  V <- matrix(rho[k], ndim, ndim)
  diag(V) <- 1
  sigma <- V * 1
  rand.effect <- mvrnorm(1, mu, V)
  Cov.1[k] <- var(rand.effect)
  
}

plot(rho, Cov.1, pch=3, cex=0.5)

abline(lm(Cov.1 ~rho))

### under the null, correlation affects enrichment analysis
rho <- 0.8
ndim <- 200

mu <- rep(0, ndim)
V <- matrix(rho, ndim, ndim)
diag(V) <- 1
sigma <- V*1

correlated.z <- mvrnorm(1, mu, V)
uncor.z <- rnorm(ndim)

beta0 <- 0
pi.uncor <- exp(beta0 + uncor.z)/( 1 + exp(beta0 + uncor.z) )
pi.cor <- exp(beta0 + correlated.z) / ( 1 + exp(beta0 + correlated.z))

GO.uncor <- rbinom(ndim, size=1, pi.uncor)
sum(GO.uncor)
GO.cor <- rbinom(ndim, size=1, pi.cor)
sum(GO.cor)


x <- rnorm(100,2,1)
y <- 1 + 2*x + rand.effect
y2 <- 1 + 2*x + rnorm(100)

summary(lm(y2~ x))




