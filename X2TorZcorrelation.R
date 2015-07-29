
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
  
  t.stat <- (mean1 - mean2)/sqrt(s1/(n/2) + s2/(n/2)) # the test statistics
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


mu1 <- c(0, 0); mu2 <- c(0, 0)
sigma1 <- c(23, 3); sigma2 <- c(23, 3)
rho <- 0.5; n <- 100; nreps <- 100;

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


