
solve.equations <- function(t.val, samp.rho)
{
  n <- length(t.val)                            # number of observations
  ones <- rep(1, n)                             #  create a column of 1s
  inv.samp <- eigen(samp.rho)
  eig.value <- inv.samp$values                  # the original eigen values
  u <- inv.samp$vectors                         # the original eigen vectors
  AA <- t(ones) %*% u 
  BB <- t(t.val) %*% u
  DD <- diag(eig.value - ones)

  obj.xi <- function(xi)
  {
    new.eig.value <- 1 / ( xi*eig.value +  1- xi) # the updated eigen values
    print(xi)
    inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
    denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
    numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
    beta.temp <- numerat / denom                  # update new beta0
    CC.temp <- BB - beta.temp*AA
    sigma2.t <- drop(CC.temp %*% inv.gamma %*% t(CC.temp)/n)      # update sigma2^2_T
    a1 <- CC.temp %*% inv.gamma %*% DD %*% inv.gamma %*% t(CC.temp) 
    gamma <- diag(xi*eig.value +  1- xi)
    sigma=u%*%gamma%*%t(u)
    sigma_inv=u%*%(diag(1/(xi*eig.value +  1- xi)))%*%t(u)
    I_H=diag(n)-matrix(1,ncol=n,nrow=n)/n
    a2=sum(diag(sigma%*%I_H%*%sigma_inv%*%(samp.rho-diag(n))%*%sigma_inv%*%I_H))   
    abs(a1/sigma2.t-a2)
  }
  
  xi.update <- optimize(obj.xi, c(1e-4, 1))         # use optimize()
  xi.new <- xi.update$minimum
  new.eig.value <- 1 / ( xi.new*eig.value +  1- xi.new ) 
  inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
  denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
  numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
  beta0.new <-  numerat / denom                 # update new beta0 
  CC <- BB - beta0.new*AA
  sigma2.t.new <- CC %*% inv.gamma %*% t(CC)/n
  return(list(beta0 = beta0.new, sigma2.t = sigma2.t.new, xi = xi.new))
}


rho <- 0.5
n <- 500 
r1 <- matrix(0.5, n, n)                         # true correlation structure
diag(r1) <- 1

xi <- 0.8                                       # true xi
sigma.t <- 1                                    # true sigma2.t
beta0 <- rep(-2, n)                              # true beta0

SIGMA <- sigma.t * ( (1- xi)*diag(1, n)  + xi* r1) # the covariance
set.seed(100)
library(MASS)
tval <- mvrnorm(n=1, mu=beta0, SIGMA)           # generate the t values


system.time(answer <- solve.equations(tval, r1))
unlist(answer) 

