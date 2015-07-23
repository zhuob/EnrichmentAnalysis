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

samp.rho <- sample.correlation
t.val <- t.stat



solve.equations <- function(t.val, samp.rho)
{
  n <- length(t.val)                            # number of observations
  ones <- rep(1, n)                             #  create a column of 1s
#  beta0 <- mean(t.val)                          # initial value for beta0 

  
  ## eigen decomposition of sample correlation matrix
  ###------- Calculated ONLY ONCE, used for all the steps that follows ---------
  inv.samp <- eigen(samp.rho)

  
  eig.value <- inv.samp$values                  # the original eigen values
  u <- inv.samp$vectors                         # the original eigen vectors
  
  ## calculating the 4 matrix defined in manuscript
  AA <- t(ones) %*% u 
  BB <- t(t.val) %*% u
  DD <- diag(eig.value - ones)
  
  
  ## objective function for xi, to be optimized for xi ------------------------
  obj.xi <- function(xi)
  {
    ##  count how many times this function has been called
    # count <<- count+1
    
    ## update the kernel inverse gamma
    new.eig.value <- 1 / ( xi*eig.value +  1- xi) # the updated eigen values
#print(xi)
    inv.gamma <- diag(new.eig.value)              # the updated inverse gamma
    
    ## update beta0
    denom <-  drop(AA %*% inv.gamma %*% t(AA))    # calculate denominator
    numerat <- drop(BB %*% inv.gamma %*% t(AA))   # calculate numerator
    beta.temp <- numerat / denom                  # update new beta0 
#print(paste("beta=", beta.temp))
    
    ## caculate matrix C, defined in manuscript
    CC <- BB - beta.temp*AA
    a0 <- CC %*% inv.gamma                        # update sigma2^2_T
    sigma2.t <- drop(a0 %*% t(CC)/n)
#print(sigma2.t)
    ## optimize objective function of xi
    a1 <- a0 %*% DD %*% inv.gamma %*% t(a0)       # term on the left part of equation
    a2 <- sum( diag(inv.gamma %*% DD ) )          # the trace, right part of the equation
    (1/sigma2.t*a1 - a2)^2                        # the objective function 
  }
  
  
  xi.update <- optimize(obj.xi, c(0, 1))        # use optimize()
  xi.new <- xi.update$minimum
  
  
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
solve.equations(t.val, samp.rho)


rho <- 0.5
n <- 500 
r1 <- matrix(0.5, n, n)
diag(r1) <- 1

xi <- 0.5                                       # true xi
sigma.t <- 3                                    # true sigma2.t
beta0 <- rep(2, n)                              # true beta0

SIGMA <- sigma.t * ( (1- xi)*diag(1, n)  + xi* r1)
set.seed(100)
tval <- mvrnorm(n=1, mu=beta0, SIGMA)


solve.equations(tval, r1)








