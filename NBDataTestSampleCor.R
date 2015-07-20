
###  generate correlated negative binomial data

## input
##  r0 and p are parameters for NB distribution
##  rho is the correlation between two NB distributions
##  n is the sample size


### under this parameterization, 
#   mu = r(1-p)/p, var = r(1-p)/p^2 = mu + k*mu^2
#   p = 1/(1 + k*mu),  r = 1/k

## this will give me identically correlated NB data
## i.e., y1 and y2 have the same mu and k.



correlated.nb <- function(r, p, rho, n)
{
  r1 <- r2 <-  (1-rho)*r
  r0 <- rho*r
  
  ## k*mu is held constant under this simulation 
  # the size parameter need not to be integer.
  x0 <- rnbinom(n, size = r0, prob=p)
  x1 <- rnbinom(n, size = r1, prob=p)
  x2 <- rnbinom(n, size = r2, prob=p)
  
  # y1 ~ NB(r, p), y2 ~ NB(r, p)  they have the same distribution.
  y1 <- x0 + x1
  y2 <- x0 + x2
  
  return(rbind(y1, y2))
}

########### TRY different way of simulating correlated NB data if 


## input: y,  matrix (both control and case are included)
#         k,  common dispersion parameter
#  output: vector of 2, test statistics for gene1 and gene2

score.stat <- function(y, k)
{
  n <- dim(y)[2]
  mean1 <- apply(y[, 1:(n/2)], 1, mean)
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)
  
  a1 <- (mean1-mean2)^2 ## denominator
  
  a2 <- (mean1 + mean2)
  a3 <- k*mean1 + 1
  a4 <- k*mean2 + 1
  a5 <- a3 + a4
  
  ## numerator part
  numerat <- 4/n * a2 *a3 *a4/a5
  
  # the test statistic 
  # u <- a1/numerat
  ## square root of test statistic
  u <- (mean1-mean2)/sqrt(numerat)
  
  return(u)
  # return(sqrt(u))  # vector of 2
  
}



## sample correlation, statistics correlation and true correlation
## if you want DE effects, just set mu1 != mu2
# input:  k, dispersion parameter 
#         rho, true correlation of genes
#         mu1 mean for case, mu2 for control
#         n, total sample size
#         nreps, number of score statistics to be generated

# output:  two numbers

sample.stat.correlation <- function(k, rho,  mu1, mu2, n, nreps)

  {
  
  r <- 1/k    # size parameter
  p1 <- 1/(1 + k*mu1) # probability parameter
  p2 <- 1/(1 + k*mu2)
  
  #  store the score test statistics
  u.stat <- matrix(NA, nrow=nreps, ncol=2)
  
  rho.sample <- c()
  
  for ( i in 1:nreps)
  {
    ## treatment 
    y.t <- correlated.nb(r, p = p1, rho = rho, n/2)
    ## control
    y.c <- correlated.nb(r, p = p2, rho = rho, n/2)
    
    # expression data of two genes
    y <- cbind(y.t, y.c)
    
    ## correlations of gene 1 and 2, for case/control respectively
    rho1 <- cor(t(y.t))[1, 2]; rho2 <- cor(t(y.c))[1, 2]
    
    # use average correlation as sample correlation
    rho.sample[i] <- mean(rho1, rho2)
    
    ## score test statistic
    u.stat[i, ] <- score.stat(y, k)
  }
  
  ## average sample correlation
  rho.ave <- mean(rho.sample)
  
  ## correlation between test statistics
  rho.score <- cor(u.stat)[1, 2]
  
  rho.list <- list(sample.cor= rho.ave, score.stat.cor = rho.score)
  return(rho.list)  # two numbers 
  
}



### THIS SIMULATION 
# 1.  For the same gene, treatment and control have the SAME k
#     (or r under this parameterization) because I derived test
#     statisitc under this assumption.
# 2.  If there is DE, allow p (or mu) to be different between Case/Control.
# 3.  we set within case (or control), mu and k are the same

k <- 0.001  ## dispersion 
mu1 <- 120  # group mean for case
mu2 <- 100  # group mean for control
n <- 100  ## total sample size (two groups)

nreps <- 100




rho <- seq(0.01, 0.99, by = 0.01)

# store the correlations in a data frame
store.corr <- data.frame(matrix(NA, length(rho), 3))
colnames(store.corr) <- c("sample", "score.stat", "true")


for ( i in 1:length(rho))
{
  a <- sample.stat.correlation(k, rho[i], mu1, mu2, n, nreps)
  
  store.corr[i, 1] <- a$sample.cor
  store.corr[i, 2] <- a$score.stat.cor
  store.corr[i, 3] <- rho[i]
  
}

library(reshape2)
store.corr1 <- melt(store.corr, id="true")
library(ggplot2)
ggplot(data = store.corr1, aes(x=true, y=value) ) +
  geom_point(aes(colour=variable)) + 
  labs(x="true", y="correlation")




