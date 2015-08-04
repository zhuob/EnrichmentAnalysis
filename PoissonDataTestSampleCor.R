

# ## pooled correlation coefficient or separate correlation?
# 
# ## they look the same
# n <- 100
# lambda <- 100
# r <- 0.5
# delta <- 0
# r1 <- r2 <- c()
# for ( i in 1:100)
# {
#   # case
#   y.t1 <- correlated.poisson(lambda, r, n/2)
#   # control
#   y.t2 <- correlated.poisson(lambda + delta, r, n/2)
#   
#   rm1 <- y.t1[1, ] - mean(y.t1[1, ])
#   rm2 <- y.t1[2, ] - mean(y.t1[2, ])
#   
#   rm3 <- y.t2[1, ] - mean(y.t2[1, ])
#   rm4 <- y.t2[2, ] - mean(y.t2[2, ])
#   
#   ## correlation by pooled data  
#   r1[i] <- cor(c(rm1, rm3), c(rm2, rm4))
#   
#   # correlation calculated by separate groups
#   rho1 <- cor(t(y.t1))[1,2]
#   rho2 <- cor(t(y.t2))[1,2]
#   r2[i] <- mean(c(rho1, rho2))
# }
# 
# plot(r1, r2)
# 



#######################################################
##### Poisson Regression: Score test statistic #####
#######################################################

correlated.poisson <- function(lambda, r, n)
 {
 lambda1 <- lambda2 <-  (1-r)*lambda 
 lambda0 <- r*lambda
 
 x0 <- rpois(n, lambda0)
 x1 <- rpois(n, lambda1)
 x2 <- rpois(n, lambda2)
 
 y1 <- x1 + x0
 y2 <- x2 + x0
 
 return(rbind(y1, y2))
}



## calculate score test stat

## square root of score test statistic does not preserve the correlation !!
score.stat1 <- function(y)
{
  n <- dim(y)[2]
  mean1 <- apply(y[, 1:(n/2)], 1, mean)
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)
  
  ## if I use this expression, the positiveness of u is messed up????
   u <- sqrt( n/2*(mean1-mean2)^2/(mean1+ mean2) )
  
  return(u)
  
}


## or we can use this expression
score.stat <- function(y)
{
  n <- dim(y)[2]
  mean1 <- apply(y[, 1:(n/2)], 1, mean)
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)
  
  ## if I use this expression, the positiveness of u is messed up????
 # u <- sqrt( n/2*(mean1-mean2)^2/(mean1+ mean2) )
  
  ## only this work if want perfect correlation
  u <- sqrt(n/2)*(mean1-mean2)/(sqrt(mean1+ mean2))
  ###########################################################################
#  It looks it doesn't matter if I change the denominator 
  ###########################################################################
#  u <- sqrt(n/2)*(mean1-mean2)/((mean1 ^2+ mean2)^2 - mean1*mean2 + 12)^(1/3)
  
  return(u)
  
}


## wald test statistic 
wald.stat <- function(y)
{
  n <- dim(y)[2]
  mean1 <- apply(y[, 1:(n/2)], 1, mean)
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)
  
  sqrt(n/2 * mean1)*(log(mean1)- log(mean2))
  
}

## likelihood ratio test stat
lr.stat <- function(y)
{
  n <- dim(y)[2]
  mean1 <- apply(y[, 1:(n/2)], 1, mean)
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)
  
  a1 <- mean1*log( (mean1 + mean2)/(2*mean1) )
  a2 <- mean2*log( (mean1 + mean2)/(2*mean2) )
  
  n*(a1 + a2)
}

## sample correlation, statistics correlation and true correlation

# if DE, the expression level difference between two groups is delta
sample.stat.correlation <- function(lambda, r, n, nreps, delta, method = "Score")
{
  
  #  store the score test statistics
  u.stat <- matrix(NA, nrow=nreps, ncol=2)
  
  rho.sample <- c()
  # create a vector of score stat for each gene
  for ( k in 1: nreps)
    {
  ## treatment 
  y.t1 <- correlated.poisson(lambda, r, n/2)
  # control
  y.t2 <- correlated.poisson(lambda + delta, r, n/2)
  
  ## expression data of two genes
  y <- cbind(y.t1, y.t2)
  rho1 <- cor(t(y.t1))[1, 2]; rho2 <- cor(t(y.t2))[1, 2]
 # print(rho1)
  rho.sample[k] <- mean(c(rho1, rho2))
 # print(rho.sample[k])
  #
  if (method=="Wald")
    {  u.stat[k, ] <- wald.stat(y)
     }
    else if (method== "LR")
    {
      u.stat[k, ] <- lr.stat(y)
      
    }
      else {u.stat[k, ] <- score.stat(y)}
  
  }
  
  rho.ave <- mean(rho.sample)
  rho.score <- cor(u.stat)[1, 2]
  rho.list <- list(sample.cor= rho.ave, score.stat.cor = rho.score)
  
  return(rho.list)
}
  

# correlations


lambda <- 100 # mean 
n <- 100  # total sample size (two groups)

# score stat for each rho
nreps <- 100
## differential expression

 # correlation
rho <- seq(0, 0.99, by = 0.01)
#set.seed(123)
delta <- rnorm(length(rho), 0, 0.2*lambda)
#delta <- rep(0, length(rho))

test.method <- "Score"

# store the correlations in a data frame
store.corr <- data.frame(matrix(NA, length(rho), 3))
colnames(store.corr) <- c("sample", "test.stat", "true")

for ( i in 1:length(rho))
{
  a <- sample.stat.correlation(lambda, rho[i], n, nreps, delta[i], method=test.method)
  
  store.corr[i, 1] <- a$sample.cor
  store.corr[i, 2] <- a$score.stat.cor
  store.corr[i, 3] <- rho[i]
  
}


library(reshape2)
store.corr1 <- melt(store.corr, id="true")
library(ggplot2)
ggplot(data = store.corr1, aes(x=true, y=value) ) +
  geom_point(aes(colour=variable)) + 
  labs(x="true", y="correlation", title=test.method)
 





###########################################################################
#  SIMULATE TWO GENES, BUT HAVING DIFFERENT EXPRESSION LEVELS FOR CONTROL GROUP
###########################################################################



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
  
  desired.rho <- seq(0.01, 0.99, by = 0.01)                           # we want these
 
  ids <- c()                                                       
   for( i in 1:length(desired.rho))
   {
     ids[i] <- which.min(abs(true.rho - desired.rho[i]))                 # find the which true correlation is closest to the desired one
   }
   true.rho.new <- true.rho[ids]                                         # the true correlations
   lambda0.new <- lambda0*rep(1, length(desired.rho))                # vector of lambda0
   lambda1.new <- lambda1[ids]                                       # vector of lambda1
   lambda2.new <- lambda2[ids]                                       # vector of lambda2
   param <- data.frame(lambda0 = lambda0.new, lambda1 = lambda1.new,
                       lambda2 =lambda2.new, true.rho= true.rho.new) # the desired parameter 
  
   rm.dup <- param[!duplicated(param),]                               # remove the duplicated rows
}

## --------------------------------------  simulate the correlated data

correlated.poisson2 <- function(lambda0, lambda1, lambda2, n)        # correlated poisson
{

  x0 <- rpois(n, lambda0)
  x1 <- rpois(n, lambda1)
  x2 <- rpois(n, lambda2)
  
  y1 <- x1 + x0
  y2 <- x2 + x0
  
  return(rbind(y1, y2))
}




# -------------------------------------  calculate the correlation under the null of both genes.

sample.stat.correlation2 <- function(lambda0, lambda1, lambda2, n, nreps)
{
  
  
  u.stat <- matrix(NA, nrow=nreps, ncol=2)                           #  store the score test statistics
  rho.sample <- c()                                                  # store the sample correlation
 
   for ( k in 1: nreps)
  {
    y.t1 <- correlated.poisson2(lambda0, lambda1, lambda2, n/2)      # treatment 
    y.t2 <- correlated.poisson2(lambda0, lambda1, lambda2, n/2)      # control, DE in 
    
   
    y <- cbind(y.t1, y.t2)                                           # expression data of two genes
    rho1 <- cor(t(y.t1))[1, 2]; rho2 <- cor(t(y.t2))[1, 2]           # sample correlation within each group
    rho.sample[k] <- mean(c(rho1, rho2))                             # average them
    
    u.stat[k, ] <- score.stat(y)                                     # the score test statistics
    
  }
  
  rho.ave <- mean(rho.sample)                                        # the overall mean as the sample mean
  rho.score <- cor(u.stat)[1, 2]                                     # the correlation between the test stats.
  rho.list <- list(sample.cor= rho.ave, score.stat.cor = rho.score)
  
  return(rho.list)
}


## ----------------------------------------  SIMULATION Study -------------------------------


param <- param.setup(lambda0 = 10)                                  # the true correlation is in the 4th column
param <- as.matrix(param)
dim.param <- dim(param)[1]                                           # the number of true correlations to be simulated

store.corr <- data.frame(matrix(NA,dim.param, 3))                    # store the three correlations
colnames(store.corr) <- c("sample", "test.stat", "true")             # give them the names
 

for ( i in 1:dim.param)
{
  a <- sample.stat.correlation2(param[i, 1], param[i, 2], param[i, 3], n=100, nreps=1000)
  store.corr[i, ] <- c(a$sample.cor, a$score.stat.cor, param[i, 4])
}
 

library(reshape2)
store.corr1 <- melt(store.corr, id="true")
library(ggplot2)
ggplot(data = store.corr1, aes(x=true, y=value) ) +
  geom_point(aes(colour=variable)) + 
  labs(x="true", y="correlation")

#-------------  CONCLUSION --------------------------------
# The correlation remains 1-1 for all null cases




