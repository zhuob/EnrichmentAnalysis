## solve the equation systems ##



#' 
#' @title Get the inverse of Sigma
#'
#' @param sample.rho   sample correlation of genes
#' @param xi     what is the xi
#' @return the updated inverse sigma

#' \item{mu1} group mean from which the treatment 1 is simulated.
#' \item{mu2} group mean from which the treatment 2 is simulated
#'




solve.equations <- function(t.val, samp.rho)
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
  HH <- 1/n* matrix(1, nrow=n, ncol=n) 
  EE <- diag(1, n)- t(u) %*% HH %*% u
  
  ## objective function, to be optimized for xi ------------------------------
  obj.xi <- function(xi)
  {
    ##  count how many times this function has been called
    # count <<- count+1
    
    ## update the kernel inverse gamma
    new.eig.value <- 1 / ( xi*eig.value +  1- xi) # the updated eigen values
   # print(xi)
    inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
    
    ## update beta0
    denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
    numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
    beta.temp <- numerat / denom                  # update new beta0 
    beta.temp <- 1
    #print(paste("beta=", beta.temp))
    
    ## caculate matrix C, defined in manuscript
    CC.temp <- BB - beta.temp*AA
    sigma2.t <- drop(CC.temp %*% inv.gamma %*% t(CC.temp)/n)      # update sigma2^2_T
    #print(sigma2.t)
    
    ## optimize objective function of xi
    # term on the left hand side of equation
    a1 <- CC.temp %*% inv.gamma %*% DD %*% inv.gamma %*% t(CC.temp) 
    # right hand side of equation 
    a2.1 <- diag( xi*eig.value +  1- xi)          # original GAMMA
    a2.2 <- EE %*% inv.gamma %*% DD %*% inv.gamma %*% EE %*% a2.1
    a2 <- sum( diag( a2.2) )                      # the trace, right hand side of the equation
    (a1 - a2)^2                          # return the objective function 
  }
  
  
  ## the initial value of 0.381966 is the golden ratio  (sqrt(5)-1)/2
  xi.update <- optimize(obj.xi, c(0, 1))         # use optimize()
  xi.new <- xi.update$minimum
  # previous optim
  # xi.update <- optim(ini.xi, obj.xi, method="Brent", lower=0, upper=1)
  
#####   calculating  beta0 and sigma^2  ----------------------------------
  
  ## update the kernel inverse gamma with xi estimate 
  new.eig.value <- 1 / ( xi.new*eig.value +  1- xi.new ) 
                                                # the updated eigen values
  inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
  
  ### update beta0
  denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
  numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
  beta0.new <-  numerat / denom                 # update new beta0 
  
  
  ## update new sigma^2
  CC <- BB - beta0.new*AA
  sigma2.t.new <- CC %*% inv.gamma %*% t(CC)/n
  
##  return all estimates 
  return(list(beta0 = beta0.new, sigma2.t = sigma2.t.new, xi = xi.new))
}

#system.time(solve.equations(t.val, samp.rho))
# solve.equations(t.val, samp.rho)

#####################  a simulation study  -----------------------------
rho <- 0.5
n <- 500 
r1 <- matrix(0.5, n, n)                         # true correlation structure
diag(r1) <- 1

xi <- 0.8                                       # true xi
sigma.t <- 1                                    # true sigma2.t
beta0 <- rep(1, n)                              # true beta0

SIGMA <- sigma.t * ( (1- xi)*diag(1, n)  + xi* r1) # the covariance
#set.seed(100)
library(MASS)
tval <- mvrnorm(n=1, mu=beta0, SIGMA)           # generate the t values


solve.equations(tval, r1)

## convergence == 0 does not mean it's converged. 
fy <- function(y){ (y^2- 3*y + 2)^2}
optim(0.6, fy,  method="Brent", lower=0, upper=0.8)
optimize(fy, c(0, 1))


