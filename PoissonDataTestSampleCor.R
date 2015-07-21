#  poisson regression statistic correlation VS sample correlation
#

## HOW IS WALD STATISTIC CALCULATED? ####
#########################################
N <- 100
mu <- rep(c(5,10), each=N/2)
x <- rep(c(1,0), each=N/2)

y <- rpois(N, mu)

pois.mod <- glm(y~x, family = poisson(link = "log"))
beta0 <- summary(pois.mod)$coe[1, 1]
beta1 <- summary(pois.mod)$coe[2, 1]

### calculate the hessian matrix 
var.beta1 <- sum(exp(beta0 + beta1*x)*x^2)
var.beta0 <- sum(exp(beta0 + beta1*x))
cov.beta0beta1 <- sum(exp(beta0+ beta1*x)*x)

hessian <- matrix(c(var.beta0, cov.beta0beta1, cov.beta0beta1, var.beta1), 2, 2)
## calculate the standard deviation for each parameter. 
## It's the same as the output of glm function.
sqrt(diag(solve(hessian)))
summary(pois.mod)
  
cor.pois <- function(n, lambda1, r)  ## generate correlated poisson random variables
{
  
lambda2 <- (1-r^2)/(0.1+r^2)*lambda1 # this makes lambda3 always positive
lambda3 <- lambda1^2/(r^2*(lambda1+lambda2))-lambda1

cor <- lambda1/(sqrt((lambda1+lambda3)*(lambda1+lambda2)))

  y1 <- rpois(n, lambda1)
  y2 <- rpois(n, lambda2)
  y3 <- rpois(n, lambda3)
  
  x1 <- y1  + y2
  x2 <- y1  + y3
  
  # the correlation between x1 and x2 
  # rho <- mu1 /sqrt( (mu1 + mu3)*(mu1 + mu2))
  #  print(round(rho, 4))
  return(cbind(x1, x2))
#  return(c(lambda1, lambda2, lambda3, r))
}



num <- 20  # how many samples in each group
# in this case, gene 1 is not DE, but gene 2 is DE
# for treatment 1


source("SimulateLabData.R")

nsim <- 1000
test.correlation <- sample.correlation <- c()
r <- seq(0.1, 0.95, length=50)
trt <- rep(c(0, 1), each=num)


calculate.lambda <- function(lambda1, r)
{
  lambda2 <- (1-r^2)/(0.1+r^2)*lambda1 # this makes lambda3 always positive
  lambda3 <- lambda1^2/(r^2*(lambda1+lambda2))-lambda1
  return(c(lambda1, lambda2, lambda3, r))
}

lambda.20 <- matrix(NA, ncol=4, nrow=length(r))
for ( j in 1:length(r))
{
  lambda.20[j, ] <- calculate.lambda(37, r[j])
  
}

for( k in 1: length(r))
{
  
  t1 <- t2 <- pois.cor <-  c()
  
for ( i in 1:nsim)
{
  
  exp1 <- cor.pois(num, 100, r[k]) # pi1 = 4, p2= 5, rho = 0.2236
  # for treatment 2
  exp2 <- cor.pois(num, 37, r[k]) #pi1 = 4, p2= 20 rho = 0.2236
  
  poisdata <- rbind(exp1, exp2)
  pois.cor[i] <- calculate.cor(t(poisdata))[1, 2]
  
  trt <- rep(c(0, 1), each=num)
  pois1 <- glm(poisdata[, 1]~trt, family = poisson(link = "log"))
  #pois1 <- lm(poisdata[, 1]~trt)
  t1[i] <- summary(pois1)$coe[2, 3] # wald statistic
  
  # pois2 <- lm(poisdata[, 2]~trt)
  pois2 <- glm(poisdata[, 2]~trt, family = poisson(link = "log"))
  t2[i] <- summary(pois2)$coe[2, 3] # wald statistic
  
}
  
sample.correlation[k] <- mean(pois.cor)
test.correlation[k] <- cor(t1, t2)
}

plot(sample.correlation, r, pch=20)
points(sample.correlation, test.correlation, col="red", pch=3)

pois.wald.stat.sample.cor <- data.frame(
  wald.stat= test.correlation, samp.cor= sample.correlation, true=r
)

saveRDS(pois.wald.stat.sample.cor, "pois.regression.wald.stat.sample.cor.rds")




## pooled correlation coefficient or separate correlation?

## they look the same
n <- 100
lambda <- 100
r <- 0.5
delta <- 0
r1 <- r2 <- c()
for ( i in 1:100)
{
  # case
  y.t1 <- correlated.poisson(lambda, r, n/2)
  # control
  y.t2 <- correlated.poisson(lambda + delta, r, n/2)
  
  rm1 <- y.t1[1, ] - mean(y.t1[1, ])
  rm2 <- y.t1[2, ] - mean(y.t1[2, ])
  
  rm3 <- y.t2[1, ] - mean(y.t2[1, ])
  rm4 <- y.t2[2, ] - mean(y.t2[2, ])
  
  ## correlation by pooled data  
  r1[i] <- cor(c(rm1, rm3), c(rm2, rm4))
  
  # correlation calculated by separate groups
  rho1 <- cor(t(y.t1))[1,2]
  rho2 <- cor(t(y.t2))[1,2]
  r2[i] <- mean(c(rho1, rho2))
}

plot(r1, r2)








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


score.stat1 <- function(y) 
## y is the read count matrix, containing both case and control
{
  # score stat
  sum1 <- apply(y[, 1:(n/2)], 1, sum)
  sum2 <- apply(y[, -(1:(n/2))], 1, sum)
  
  sum3 <- apply(y, 1, sum)

  u <- (sum1-sum2)/sqrt(sum3)
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
  return(u)
  
}


## wald statistic
wald.stat <- function(y)
{
  n <- dim(y)[2]
  mean1 <- apply(y[, 1:(n/2)], 1, mean)
  mean2 <- apply(y[, -(1:(n/2))], 1, mean)
  
  sqrt(n/2 * mean1)*(log(mean1)- log(mean2))
  
}

## likelihood ratio stat
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
  u.stat[k, ] <- score.stat(y)
  
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
set.seed(123)
delta <- rnorm(length(rho), 0, 0.2*lambda)
delta <- rep(0, length(rho))

test.method <- "LR"

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
 

