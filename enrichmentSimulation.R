## 
source("SimulateLabData.R")

size <- 100
N <- 500
rho <- 0.2
prob <- 0.1
DE <- N*prob
dat <- simu.microarray(size = size, N = N,prob = prob, rho=rho)
lab1 <- dat$data


#  for N=500 and 50 DE genes, 60 genes are in the go term with 40 True DE
#  so it is enriched
goterm <- c(sample(1:DE, 0.8*DE), sample((DE + 1):N, 0.4*DE))
GO <- rep(0, 500); GO[goterm] <- 1


t.result <- function(data, quantity= "stat")
{
  n <- length(data)
  t1 <- data[1:n/2]
  t2 <- data[-(1:n/2)]
    t.test(t1, t2)$stat
}

# calculate the test stat
T1 <- apply(lab1, 1, t.result)
 id <- which(T1< qt(0.05, 98))
T1[id] <- -1*T[id]


trt <- rep(c(1, 2), each=size/2)
demean <- lab1
for ( i in 1:500)
{
  demean[i, ] <- lab1[i, ]- rep(tapply(lab1[i, ], trt, mean), each= size/2)
  
}

cor.est <- cor(t(demean))

library(regress)

summary(lm(T1~GO))
model.sim2 <- regress(T1~GO, ~cor.est)
cor.struct <- matrix(rho, N, N);diag(cor.struct) <- 1
model.sim3 <- regress(T1~GO, ~cor.struct)
model.sim2
model.sim3

names(model.sim2)
model.sim3$beta[1]

df2 <- df3 <- data.frame(matrix(ncol=6, nrow=100))
colnames(df2) <- c("int", "int.sd", "go", "go.sd", "sigma1", "In")
colnames(df3) <- c("int", "int.sd", "go", "go.sd", "sigma1", "In")

for ( k in 1:size)
{
  T1 <- lab1[, k]
  model.sim2 <- regress(T1~GO, ~cor.est)
  value2 <- c(model.sim2$beta[1], model.sim2$beta.se[1], 
              model.sim2$beta[2], model.sim2$beta.se[2],
              model.sim2$sigma[1], model.sim2$sigma[2])
  
  df2[k, ] <- value2
  cor.struct <- matrix(rho, N, N);diag(cor.struct) <- 1
  model.sim3 <- regress(T1~GO, ~cor.struct)
  
  value3 <- c(model.sim3$beta[1], model.sim3$beta.se[1], 
              model.sim3$beta[2], model.sim3$beta.se[2],
              model.sim3$sigma[1], model.sim3$sigma[2])
  df3[k, ] <- value3
}





######  try 100 simulations of estimating coefficients and variance components 
library(MASS)

mu <- GO*0.8 + 1.5 # the true mean

# covariance
correlation.true <- matrix(0.4, N, N); diag(correlation.true) <- 1
sigma <- 1.5*correlation.true

##  expression matrix 

expression.mat <- t(mvrnorm(size, mu, sigma)) 

## estimated correlation matrix
correlation.est <- cor(t(expression.mat))


df4 <- df5 <- data.frame(matrix(ncol=6, nrow=100))
colnames(df4) <- c("int", "int.sd", "go", "go.sd", "sigma1", "In")
colnames(df5) <- c("int", "int.sd", "go", "go.sd", "sigma1", "In")

for ( k in 1:size)
{
  y <- expression.mat[, k]
  model.sim2 <- regress(y~GO, ~correlation.est)
  value2 <- c(model.sim2$beta[1], model.sim2$beta.se[1], 
              model.sim2$beta[2], model.sim2$beta.se[2],
              model.sim2$sigma[1], model.sim2$sigma[2])
  
  df4[k, ] <- value2
  model.sim3 <- regress(y~GO, ~correlation.true)
  
  value3 <- c(model.sim3$beta[1], model.sim3$beta.se[1], 
              model.sim3$beta[2], model.sim3$beta.se[2],
              model.sim3$sigma[1], model.sim3$sigma[2])
  df5[k, ] <- value3
}


### conclusion::::::::::::::::::::::::::::::::::::::::::::::::::::
##  there might be computational issues for the case of true correlation matrix


###  MCMC is not an option, for 500 genes with just 1000 iterations,
##  it took alomst 8GB of memeory (more than 30 minutes)



# 
# y <- lab1[, 1]
# 
# dtlist <- list(N=length(y), y=y, GO=GO, corr= cor.struct)
# library(rstan)
# fit <- stan(file = "LMMcovarianceMatrix.stan", data=dtlist, iter=1000, chains=4)
# 
# print(fit, pars= c("beta", "tau"), probs=c(0.025, 0.975))
# fit2 <- stan(fit=fit, data=dtlist, iter=10000, chains=4)




# this estimate is quite good
################################################
# Inference for Stan model: LMMcovarianceMatrix.
# 4 chains, each with iter=1000; warmup=500; thin=1; 
# post-warmup draws per chain=500, total post-warmup draws=2000.
# 
#       mean se_mean   sd  2.5% 97.5% n_eff Rhat
# beta -0.02       0 0.12 -0.25  0.22  1471    1
# tau   1.00       0 0.03  0.94  1.07  1681    1
# 
# Samples were drawn using NUTS(diag_e) at Thu Jun  4 15:02:04 2015.
# For each parameter, n_eff is a crude measure of effective sample size,
# and Rhat is the potential scale reduction factor on split chains (at 
#                                                                   convergence, Rhat=1).



# ## regress gives me this 
################################################
# model.sim2 <- regress(y~GO, ~cor.struct)
# 
# Likelihood kernel: K = (Intercept)+GO
# 
# Maximized log likelihood with kernel K is  -193.215 
# 
# Linear Coefficients:
#   Estimate Std. Error
# (Intercept)    0.228      0.301
# GO            -0.022      0.123
# 
# Variance Coefficients:
#   Estimate Std. Error
# cor.struct    0.444      0.025
# In            0.444      0.031
# 
# 
