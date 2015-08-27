
## generate correlated Normal data

## input:  
#    mu:  a vector of 2, specifying the mean level for each gene
#   rho:  the true correlation we want to simulate 
# sigma:  a vector of 2, specifying the variance 
#     n:  the number of samples to be generated

## output: a matrix of 2 by n
#

correlated.norm <- function(mu, rho, sigma, n)
{
  library(MASS)
  sig.element <- c(sigma[1], rho*sqrt(sigma[1]*sigma[2]),     ## the elements in Sigma
                   rho*sqrt(sigma[1]*sigma[2]), sigma[2])
  sig <- matrix(sig.element, 2, 2)                            # the covariance matrix
  dat <- mvrnorm(n=n, mu = mu, Sigma = sig)                   # simulate the correlated normal data
  return(t(dat))
}


ttest.stat <- function(y)
{
  n <- dim(y)[2]                              # get the number of samples from each gene
  mean1 <- apply(y[, 1:(n/2)], 1, mean)       # the first half is from treatment
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)    # the second half is from control
  
  s1 <- apply(y[, 1:(n/2)], 1, var)           # the variance for the first half
  s2 <- apply(y[, -(1:(n/2))], 1, var)        # the variance for the second half
  
#   t1 <- t.test(y[1, 1:(n/2)], y[1, -(1:(n/2))])$stat
#   t2 <- t.test(y[2, 1:(n/2)], y[2, -(1:(n/2))])$stat
#   t.stat <- c(t1, t2)
#   
  t.stat <- (mean1 - mean2)/sqrt(s1/(n/2) + s2/(n/2))  # the test statistics
  return(t.stat)
}


ztest.stat <- function(y, sigma1, sigma2)
{
  n <- dim(y)[2]                              # get the number of samples from each gene
  mean1 <- apply(y[, 1:(n/2)], 1, mean)       # the first half is from treatment
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)    # the second half is from control
  
  z.stat <- (mean1-mean2)/sqrt(sigma1/(n/2) + sigma2/(n/2))
  return(z.stat)
}


sample.stat.correlation <- function(mu1, sigma1, mu2, sigma2, rho, n, nreps)
{
  
  stat <- matrix(NA, nrow=nreps, ncol=2)                     # matrix to store the t.stat
  rho.sample <- c()
 
  for ( k in 1: nreps)                                       # create a vector of score stat for each gene
  {
    y.t1 <- correlated.norm(mu1, rho, sigma1, n/2)           # treatment 
    y.t2 <- correlated.norm(mu2, rho, sigma2, n/2)           # control
    
    y <- cbind(y.t1, y.t2)                                   # expression data of two genes
    rho1 <- cor(t(y.t1))[1, 2]; rho2 <- cor(t(y.t2))[1, 2]   # sample correlation for each group
    rho.sample[k] <- mean(c(rho1, rho2))                     # mean of sample correlation
 
#    stat[k, ] <- ztest.stat(y, sigma1, sigma2)              # calculate the z stat
    stat[k, ] <- ttest.stat(y)                               #  calculate the t stat
  }
  
  rho.ave <- mean(rho.sample)                                # calculate the mean of the sampel rhos
  rho.t <- cor(stat)[1, 2]                                   # calculate the correlation of statistics
  rho.list <- list(sample.cor= rho.ave, stat.cor = rho.t)
    
  return(rho.list)  
}


mu1 <- c(10, -10); mu2 <- c(10, -20)
sigma1 <- c(.1, .3); sigma2 <- c(.2, .6)
rho <- 0.5; n <- 1000; nreps <- 100;

rho <- seq(0, 0.99, by = 0.01)
store.corr <- data.frame(matrix(NA, length(rho), 3))
colnames(store.corr) <- c("sample", "test.stat", "true")


for ( i in 1:length(rho))
{
  a <- sample.stat.correlation(mu1, sigma1, mu2, sigma2, rho[i], n, nreps)
  
  store.corr[i, 1] <- a$sample.cor
  store.corr[i, 2] <- a$stat.cor
  store.corr[i, 3] <- rho[i]
  
}


library(reshape2)
store.corr1 <- melt(store.corr, id="true")
library(ggplot2)
ggplot(data = store.corr1, aes(x=true, y=value) ) +
  geom_point(aes(colour=variable)) + 
  labs(x="true", y="correlation")




##########################
# If we look at the $z$-test statistic and $t$-test statistic, 
# \item if $m$ and $n$ are small, then $S_{X_1}^2$ and $S_{Y_1}^2$ cannot be accurately estimated
# \item if $m$ and $n$ are large, then the denominator will be very small, in which case a slight difference between the sample variance  $S_{X_1}^2$ and true variance  $\sigma^2_1$ will augment the difference of test statistics. 


nsim <- 100
t.stat  <- s.mat<- m.mat <- matrix(NA, nsim, 4)
t.stat <- s.mat <- data.frame(t.stat)
colnames(t.stat) <- c("t.test1", "t.test2", "z.test1", "z.test2")
colnames(s.mat) <- c("s11", "s21", "s12", "s22")

for ( i in 1:nsim){
  
mu1 <- c(15, 10); mu2 <- c(10, 10)
sigma1 <- c(1, 1); sigma2 <- c(1, 1)
rho <- 0.5; n <- 40;

y.t1 <- correlated.norm(mu1, rho, sigma1, n/2)           # treatment 
y.t2 <- correlated.norm(mu2, rho, sigma2, n/2)           # control

y <- cbind(y.t1, y.t2)   

mean1 <- apply(y[, 1:(n/2)], 1, mean)       # the first half is from treatment
mean2 <- apply(y[, -(1:(n/2))], 1, mean)    # the second half is from control

s1 <- apply(y[, 1:(n/2)], 1, var)           # the variance for the first half
s2 <- apply(y[, -(1:(n/2))], 1, var)        # the variance for the second half

m.mat[i, ] <- c(mean1, mean2)
s.mat[i, ] <- c(s1, s2)

t.stat[i,1:2 ] <- (mean1 - mean2)/sqrt(s1/(n/2) + s2/(n/2))  # the test statistics
t.stat[i,3:4 ] <- (mean1-mean2)/sqrt(sigma1/(n/2) + sigma2/(n/2))
}

par(mfrow=c(1, 2))
plot(t.stat[, 1], t.stat[, 3], pch=20, xlab="t-test stat", ylab="z-test stat", main="DE")
plot(t.stat[, 2], t.stat[, 4], pch=20, xlab="t-test stat", ylab="z-test stat", main="Non-DE")
cor(t.stat)


x1 <- rnorm(100, 15, 1)
x2 <- 15 + rnorm(100)

cor(x1,x2)

plot(x1, x2, pch=20)

## Why is this happening

barx <- varx <- c()
for ( i in 1:100)
{
  x <- rnorm(100)
  
  barx[i] <- mean(x)
  varx[i] <- var(x)
  
}

cor(barx, varx)


rho.vector <- seq(0.05, 0.95, by = 0.01)                      # the true correlation
sd.mean.corr <- c()                                           # the sample correlation of std of g1 and mean of g2
mean.cor <- sd.corr <- c()

for ( k0 in 1:length(rho.vector))
{
  rho <- rho.vector[k0]

  sd.ave <- mu.ave <- matrix(NA, 100, 2)
  for ( k in 1:100)
    {
    mu <- c(0, 0);  sigma <- matrix(c(1, rho, rho, 1), ncol=2, nrow=2)
    xx <- mvrnorm(1000, mu, sigma )
    var1 <- apply(xx, 2, var)
  
    mu2 <- c(0, 0);  sigma2 <- matrix(c(1, rho, rho, 1), ncol=2, nrow=2)
    xx2 <- mvrnorm(1000, mu2, sigma2 )
    var2 <- apply(xx2, 2, var)

    xy <- rbind(xx, xx2)
  
    sd.ave[k, ] <- 1/sqrt((var1 + var2)/2)                         # mean sd as the pooled sd

    mu.ave[k, ] <- apply(xx, 2, mean) - apply(xx2, 2, mean)
    }

   sd.mean.corr[k0] <- cor(sd.ave[, 1], mu.ave[, 2])
   sd.corr[k0] <- cor(sd.ave)[1, 2]
   mean.cor[k0] <- cor(mu.ave)[1, 2]
}

plot(rho.vector, sd.mean.corr)
abline(0, 0)
plot(rho.vector, sd.corr)
abline(0, 1)
plot(rho.vector, mean.cor)
abline(0, 1)


n <- 10:170

a <- sqrt(n-1)*(gamma(n-1.5))/(gamma(n-1))
b <- sqrt((n-1)/(n-2))
#c(a, b)
plot(n, a, pch=20, cex=0.5)
points(n, b, cex=0.5, pch=3)
plot(n, n*(-a^2+b^2))

n <- 20:500

a <- sqrt((n-1)/(n-2))
b <-  choose(2*n-4, n-2)/4^(n-2)*sqrt((n-1)*pi)
which(is.nan(b))
p <-  a^2 -b^2


plot(n, n*p, pch=20, cex=0.3)
range(n*p)

############ now the correlation between test statistic is related to both 
###########  rho0 = cor(z1, z2) and rho1 = cor(1/s1, 1/s2)

# simulate the correlation of 1/s1 and 1/s2

##  calculate the inverse of standard errors
inv_sd <- function(y)
{
  n <- dim(y)[2]                              # get the number of samples from each gene
  var1 <- apply(y[, 1:(n/2)], 1, var)         # the variance for the first half
  var2 <- apply(y[, -(1:(n/2))], 1, var)      # the variance for the second half
  
  inv_s1 <- 1/sqrt((var1[1] + var2[1])/2)                 # pooled variance for gene 1
  inv_s2 <- 1/sqrt((var1[2] + var2[2])/2)                 # pooled variance for gene 2
  
#  inv_s1 <- (var1[1] + var2[1])/2                 # pooled variance for gene 1
#  inv_s2 <- (var1[2] + var2[2])/2                 # pooled variance for gene 2
  
  return(c(inv_s1, inv_s2))
}



# the correlation of sample and the correlation of inverse_standard_errors
sample.inv_sd.correlation <- function(mu1, sigma1, mu2, sigma2, rho, n, nreps)
{
  
  stat <- matrix(NA, nrow=nreps, ncol=2)                     # matrix to store the t.stat
  rho.sample <- c()
  
  for ( k in 1: nreps)                                       # create a vector of score stat for each gene
  {
    y.t1 <- correlated.norm(mu1, rho, sigma1, n/2)           # treatment 
    y.t2 <- correlated.norm(mu2, rho, sigma2, n/2)           # control
    
    y <- cbind(y.t1, y.t2)                                   # expression data of two genes
    rho1 <- cor(t(y.t1))[1, 2]; rho2 <- cor(t(y.t2))[1, 2]   # sample correlation for each group
    rho.sample[k] <- mean(c(rho1, rho2))                     # mean of sample correlation
    
    #    stat[k, ] <- ztest.stat(y, sigma1, sigma2)              # calculate the z stat
    stat[k, ] <- inv_sd(y)
  }
  
  rho.ave <- mean(rho.sample)                                # calculate the mean of the sampel rhos
  rho.s <- cor(stat)[1, 2]                                   # calculate the correlation of statistics
  rho.list <- list(sample.cor= rho.ave, inv_sd.cor = rho.s)
  
  return(rho.list)  
}



mu1 <- c(10, 10); mu2 <- c(10, 10)
sigma1 <- c(1, 3); sigma2 <- c(1, 3)
rho <- 0.5; n <- 1000; nreps <- 100;

rho <- seq(0, 0.99, by = 0.01)
store.corr <- data.frame(matrix(NA, length(rho), 3))
colnames(store.corr) <- c("sample", "test.stat", "true")


for ( i in 1:length(rho))
{
  a <- sample.inv_sd.correlation(mu1, sigma1, mu2, sigma2, rho[i], n, nreps)
  
  store.corr[i, 1] <- a$sample.cor
  store.corr[i, 2] <- a$inv_sd.cor
  store.corr[i, 3] <- rho[i]
  
}


library(reshape2)
store.corr1 <- melt(store.corr, id="true")
library(ggplot2)
ggplot(data = store.corr1, aes(x=true, y=value) ) +
  geom_point(aes(colour=variable)) + 
  labs(x="true", y="correlation")




