
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
#         k,  vector of 2, dispersion parameter vector, k[1] is for gene1 and k[2] is for gene 2
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

k <- rep(0.001, 2)  ## dispersion 
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




## set up parameters for desired correlation
##-------------------------------------------------------------------------
# the correlation is based on the fact that, x0, x1 and x2 are independent 
# poisson, then Y1 = x0 + x1 and Y2 = x0 + x2 are correlated.
##-------------------------------------------------------------------------

                                  
param.setup <- function(lambda0, lambda2=seq(1, 10000, by =1), rate=1/10)
{
  
  # lambda0 <- 100      
  # lambda2=seq(1, 10000, by =1)
  # rate=20
  lambda1 <- rexp(length(lambda2), rate)                             # lambda1 is simulated from exponential distribution
  true.rho <- lambda0/sqrt((lambda0 + lambda1)*(lambda0 + lambda2))  # the true correlation is calculated
  
  desired.rho <- seq(0.01, 0.99, by = 0.01)                          # we want these
  
  ids <- c()                                                       
  for( i in 1:length(desired.rho))
  {
    ids[i] <- which.min(abs(true.rho - desired.rho[i]))             # find the which true correlation is closest to the desired one
  }
  true.rho.new <- true.rho[ids]                                     # the true correlations
  lambda0.new <- lambda0*rep(1, length(desired.rho))                # vector of lambda0
  lambda1.new <- lambda1[ids]                                       # vector of lambda1
  lambda2.new <- lambda2[ids]                                       # vector of lambda2
  param <- data.frame(lambda0 = lambda0.new, lambda1 = lambda1.new,
                      lambda2 =lambda2.new, true.rho= true.rho.new) # the desired parameter 
  
  rm.dup <- param[!duplicated(param),]                              # remove the duplicated rows
}


## simulate correlated NB data with desired correlation
# if we know x0 ~ NB(lambda0, p), then kappa0 = 1/lambda0, mu0 = (1/p -1)/kappa0
# We want all the p's be the same for x0, x1 and x2
correlated.nb2 <- function(lambda0, lambda1, lambda2, mu0, n)
{
  kappa <- 1/c(lambda0, lambda1, lambda2)             # change the parameterization
  mu1 <- kappa[1]/kappa[2]*mu0                        # mean parameter for x1
  mu2 <- kappa[1]/kappa[3]*mu0                        # mean parameter for x2
  
  x0 <- rnbinom(n, mu = mu0, size=1/kappa[1])
  x1 <- rnbinom(n, mu = mu1, size=1/kappa[2])
  x2 <- rnbinom(n, mu = mu2, size=1/kappa[3])
  
  y1 <- x0 + x1
  y2 <- x0 + x2
  return(rbind(y1, y2))
}

lambda <- dat[66, 1:3]

##  calculate the sample and test stattistic correlation
sample.stat.correlation2 <- function(lambda0, lambda1, lambda2, mu0, n, nreps, kappa)
{
  
  
  u.stat <- matrix(NA, nrow=nreps, ncol=2)                           #  store the score test statistics
  rho.sample <- c()                                                  # store the sample correlation
  
  for ( k in 1: nreps)
  {
    ##  -----------------------------------------------------------------------
    # MY score test statistic is good only under the NULL 
    
    y.t1 <- correlated.nb2(lambda0, lambda1, lambda2, mu0, n/2)      # treatment 
    y.t2 <- correlated.nb2(lambda0, lambda1, lambda2, mu0, n/2)      # control, DE in 
    
    
    y <- cbind(y.t1, y.t2)                                           # expression data of two genes
    rho1 <- cor(t(y.t1))[1, 2]; rho2 <- cor(t(y.t2))[1, 2]           # sample correlation within each group
    rho.sample[k] <- mean(c(rho1, rho2))                             # average them
    
    u.stat[k, ] <- score.stat(y, kappa)                              # the score test statistics
    
  }
  
  rho.ave <- mean(rho.sample)                                        # the overall mean as the sample mean
  rho.score <- cor(u.stat)[1, 2]                                     # the correlation between the test stats.
  rho.list <- list(sample.cor= rho.ave, score.stat.cor = rho.score)
  
  return(rho.list)
}


## ----------------------------------------  SIMULATION Study ---------------------



param <- param.setup(lambda0 = 100)                                  # the true correlation is in the 4th column
param <- as.matrix(param)
dim.param <- dim(param)[1]                                           # the number of true correlations to be simulated

store.corr <- data.frame(matrix(NA,dim.param, 3))                    # store the three correlations
colnames(store.corr) <- c("sample", "test.stat", "true")             # give them the names


for ( i in 1:dim.param)
{
  kappa <- c( 1/(param[i, 1] + param[i, 2]), 1/(param[i, 1] + param[i, 3]))
  a <- sample.stat.correlation2(param[i, 1], param[i, 2], param[i, 3], mu0 = 100, n=1000, nreps=100, kappa)
  store.corr[i, ] <- c(a$sample.cor, a$score.stat.cor, param[i, 4])
}


library(reshape2)
store.corr1 <- melt(store.corr, id="true")
library(ggplot2)
ggplot(data = store.corr1, aes(x=true, y=value) ) +
  geom_point(aes(colour=variable)) + 
  labs(x="true", y="correlation")


## this small experiment shows that the socre.stat is basically the same as a test stat from glm.nb() function
n <- 100
mu0 <- 20
mu1 <- 30
kappa <- c(0.01, 0.03)
x0 <- rnbinom(n, mu = mu0, size=1/kappa[1])
x1 <- rnbinom(n, mu = mu1, size=1/kappa[2])
c(mean(x1), var(x1))
c(mean(x0), var(x0))

library(MASS)
ids <- rep(c(1, 0), each=n/2)
t1 <- summary(glm.nb(x0 ~ ids, link = log))$coefficients[2, 3]
t2 <- summary(glm.nb(x1 ~ ids, link = log))$coefficients[2, 3]
c(t1, t2)

y <- rbind(x0, x1)
score.stat(y, kappa)










