## solve the equation systems ##




#' @title Get the inverse of Sigma
#' 
#' @param sample.rho   sample correlation of genes
#' @param xi     what is the xi
#' @return the updated inverse sigma

# 
# solve.equations <- function(t.val, samp.rho)
# {
#   n <- length(t.val)                            # number of observations
#   ones <- rep(1, n)                             #  create a column of 1s
# 
#   
#   ## eigen decomposition of sample correlation matrix
#   ###------- Calculated ONLY ONCE, used for all the steps that follows ---------
#   inv.samp <- eigen(samp.rho)
# 
#   
#   eig.value <- inv.samp$values                  # the original eigen values
#   u <- inv.samp$vectors                         # the original eigen vectors
#   
#   ## calculating the 4 matrices defined in manuscript
#   AA <- t(ones) %*% u 
#   BB <- t(t.val) %*% u
#   DD <- diag(eig.value - ones)
#   HH <- 1/n* matrix(1, nrow=n, ncol=n) 
#   EE <- diag(n)- t(u) %*% HH %*% u
#   
#   ## objective function, to be optimized for xi ------------------------------
#   obj.xi <- function(xi)
#   {
#     ##  count how many times this function has been called
#     # count <<- count+1
#     
#     ## update the kernel inverse gamma
#     new.eig.value <- 1 / ( xi*eig.value +  1- xi) # the updated eigen values
#     print(xi)
#     inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
#     
#     ## update beta0
#     denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
#     numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
#     beta.temp <- numerat / denom                  # update new beta0 
#     
#     ## caculate matrix C, defined in manuscript
#     CC.temp <- BB - beta.temp*AA
#     sigma2.t <- drop(CC.temp %*% inv.gamma %*% t(CC.temp)/(n-1))      # update sigma2^2_T
#  
#     ## optimize objective function of xi
#     # term on the left hand side of equation
#     a1 <- CC.temp %*% inv.gamma %*% DD %*% inv.gamma %*% t(CC.temp) 
#     
#     # right hand side of equation 
#     gamma <- diag(xi*eig.value + 1-xi)
#     sigma <- u %*% gamma %*% t(u)
#     I_H <- diag(n) - matrix(1/n, n, n)
#     sigma.inv <- u %*% inv.gamma %*% t(u)
#     a2 <- sum(diag(sigma %*% I_H %*% sigma.inv %*% (samp.rho- diag(n)) %*% sigma.inv %*% I_H  ))
#     a2.1 <- diag( xi*eig.value +  1- xi)          # original GAMMA
#     a2.2 <- a2.1 %*% EE %*% inv.gamma %*% DD %*% inv.gamma %*% EE 
#     a2 <- sum( diag( a2.2) )                      # the trace, right hand side of the equation
#     abs(a1/sigma2.t - a2)                         # return the objective function 
#   }
#   
#   ## the initial value of 0.381966 is the golden ratio  (sqrt(5)-1)/2
#   xi.update <- optimize(obj.xi, c(0, 1))         # use optimize()
#   xi.new <- xi.update$minimum
#   # previous optim
#   # xi.update <- optim(ini.xi, obj.xi, method="Brent", lower=0, upper=1)
#   
# #####   calculating  beta0 and sigma^2  ----------------------------------
#   
#   ## update the kernel inverse gamma with xi estimate 
#   new.eig.value <- 1 / ( xi.new*eig.value +  1- xi.new ) 
#                                                 # the updated eigen values
#   inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
#   
#   ### update beta0
#   denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
#   numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
#   beta0.new <-  numerat / denom                 # update new beta0 
#   
#   
#   ## update new sigma^2
#   CC <- BB - beta0.new*AA
#   sigma2.t.new <- CC %*% inv.gamma %*% t(CC)/n
#   
# ##  return all estimates 
#   return(list(beta0 = beta0.new, sigma2.t = sigma2.t.new, xi = xi.new))
# }



solve.equations <- function(t.val, samp.rho, start, end, n.or.nminus1)
{
  n <- length(t.val)                            # number of observations
  ones <- rep(1, n)                             #  create a column of 1s
  
  
  ## eigen decomposition of sample correlation matrix
  ###------- Calculated ONLY ONCE, used for all the steps that follows ---------
  inv.samp <- eigen(samp.rho)
  
  eig.value <- inv.samp$values                  # the original eigen values
  u <- inv.samp$vectors                         # the original eigen vectors
  
  ## calculating the 4 matrices defined in manuscript
  AA <- t(ones) %*% u 
  BB <- t(t.val) %*% u
  DD <- diag(eig.value - ones)
  
  ## objective function, to be optimized for xi ------------------------------
  obj.xi <- function(xi)
  {

    ## update the kernel inverse gamma
    new.eig.value <- 1 / ( xi*eig.value +  1- xi) # the updated eigen values
    inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
    
    ## update beta0
    denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
    numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
    beta.temp <- numerat / denom                  # update new beta0 

    
    ## caculate matrix C, defined in manuscript
    CC.temp <- BB - beta.temp*AA
    sigma2.t <- drop(CC.temp %*% inv.gamma %*% t(CC.temp)/(n.or.nminus1))          # update sigma2^2_T
    # what if I replace with true value
    # sigma2.t <- sigma.t
     #print(sigma2.t)
    
    ## optimize objective function of xi
    a1 <- CC.temp %*% inv.gamma %*% DD %*% inv.gamma %*% t(CC.temp)               # term on the left part of equation
    
    # terms on the right hand side
    gamma <- diag(xi*eig.value +  1- xi)
    sigma=u%*%gamma%*%t(u)
    sigma_inv=u%*%(diag(1/(xi*eig.value +  1- xi)))%*%t(u)
    I_H=diag(n)-matrix(1,ncol=n,nrow=n)/n
    a2=sum(diag(sigma%*%I_H%*%sigma_inv%*%(samp.rho-diag(n))%*%sigma_inv%*%I_H))   
 #   print( c(xi, beta.temp, sigma2.t, abs(a1/sigma2.t - a2))  )
    
       abs(a1/sigma2.t - a2)                                 # return the objective function 
    }
  
  
  ## the initial value of 0.381966 is the golden ratio  (sqrt(5)-1)/2
  xi.update <- optimize(obj.xi, c(start, end))               # use optimize()
  xi.new <- xi.update$minimum

  #####   calculating  beta0 and sigma^2  ----------------------------------
  
  ## update the kernel inverse gamma with xi estimate 
  new.eig.value <- 1 / ( xi.new*eig.value +  1- xi.new )     # the updated eigen values
  inv.gamma <- diag(new.eig.value)                           # the updated inverse gamma
  
  ### update beta0
  denom <-  drop(AA %*% inv.gamma %*% t(AA))                 # calculate denominator
  numerat <- drop(BB %*% inv.gamma %*% t(AA))                # calculate numerator
  beta0.new <-  numerat / denom                              # update new beta0 
  
  ## update new sigma^2
  CC <- BB - beta0.new*AA
  sigma2.t.new <- CC %*% inv.gamma %*% t(CC)/n.or.nminus1
  
  ##  return all estimates 
  return(list(beta0 = beta0.new, sigma2.t = sigma2.t.new, xi = xi.new))
}




#####################  a simulation study  -----------------------------


library(doParallel)                             # parallel computing packages
library(foreach)

cl <- makeCluster(3)                            # request 3 cores to do the simulation
registerDoParallel(cl)                          # Register cluster
getDoParWorkers()                               # Find out how many cores are being used

system.time(
  fti <- foreach(i = 1:2, .combine = rbind) %dopar% {
    
    N <- m <- 500 # number of genes
   # n <- 20 # biological sample in each group
    library(MASS)
    
    ## generate correlation matrix
    rho <- 0.4
    cor.struct <- matrix(rho, N, N);diag(cor.struct) <- 1
    sigma <- 1.5*cor.struct
    dat <- mvrnorm(1000, mu=rep(0, m), sigma)
    
    sigma1 <- cor(dat)
    sigma2 <- 1.5*sigma1
    dat <- mvrnorm(1000, mu=rep(0, m), sigma2)
    
    sim.sigma <- cor(dat)
    
    rho <- 0.5
    n <- 500 
    
    # r1 <- matrix(0.5, n, n)                         # true correlation structure
    r1 <- sim.sigma
    diag(r1) <- 1
    
    xi <- 0.8                                       # true xi
    sigma.t <- 4                                    # true sigma2.t
    beta0 <- rep(-2, n)                              # true beta0
    
    SIGMA <- sigma.t * ( (1- xi)*diag(1, n)  + xi* r1) # the covariance
    # set.seed(100)
    library(MASS)
    t.val <- mvrnorm(n=1, mu=beta0, SIGMA)           # generate the t values
    
    answer <- solve.equations(t.val, r1, 0, 1, n)
    return(unlist(answer))
  }
)





# ## convergence == 0 does not mean it's converged. 
# fy <- function(y){ (y^2- 3*y + 2)^2}
# optim(0.6, fy,  method="Brent", lower=0, upper=0.8)
# optimize(fy, c(0, 1))

BetaSigmaXi <- read.table("BetaSigmaXi3.txt", header=T)
#View(BetaSigmaXi)

library(ggplot2)
library(reshape2)
library(dplyr)

BetaSigmaXi2 <- melt(BetaSigmaXi)
ggplot(BetaSigmaXi2, aes(x = value)) +
  geom_histogram() +
  facet_wrap(~variable, scale= "free")

BetaSigmaXi2 %>%
  group_by(variable)  %>%
  summarize(mean(value), sd(value))

sigma <- seq(0.1, 10, by = 0.1)
plot(sigma, BetaSigmaXi$sigma2.t, pch=20, cex=0.5)
abline(0 ,1)















