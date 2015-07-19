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

sample.rho <- sample.correlation

xi <- 0.5

### maybe I don't need this one
sigma.matrix <- function(samp.rho, xi)
{
  ## obtain the dimension of sample
  n <- dim(samp.rho)[1]
  return( (1-xi)* diag(1, n) + xi*samp.rho )
}


## update inv.sigma for each new xi
update.inv.sigma <- function(inv.samp, xi)
{
  # the original eigen values
  eig.value <- inv.samp$values  
  # the original eigen vectors, this is fixed at each iteration
  u <- inv.samp$vectors
  
  ## the updated eigen values
  new.eig.value <- 1 / ( xi*eig.value +  1- xi) 
  ## the updated inverse sigma
  new.inv.sigma <- u %*% diag(new.eig.value) %*% t(u)
  
  return(new.inv.sigma)
}


## solve for beta0

update.beta0 <- function(xi)
{
  ### the inv.samp is the inverse of sample correlation
  inv.sigma <- update.inv.sigma(inv.samp, xi)
  
  # calculate (1^T sig^{-1} 1)
  denom <-  drop(ones %*% inv.sigma %*% ones)
  # calculate t^T sig^{-1} 1
  numerat <- drop(t.val %*% inv.sigma %*% ones)
  ## update new beta0 
  beta0 <- 1/denom * numerat 
  
  return(beta0)
}



## objective function for xi, to be optimized for xi
obj.xi <- function(xi)
{
  ## get the inverse of SIGMA ------------------
  inv.sigma <- update.inv.sigma(inv.samp, xi)
  
  ## vectorize beta0
  beta0 <- beta0*ones
  
  a1 <- (t.val- beta0)
  a2 <- samp.rho - diag(1, n)
  a3 <- inv.sigma %*% a2
  
  ## calculate the trace of a3
  tr.a3 <- sum(diag(a3))
  
  ## the objective function 
  (1/sigma2.t * a1 %*% inv.sigma %*% a2 %*% inv.sigma %*% (a1) - tr.a3)^2
}


###  if I know beta0 and xi, I can calculate sigma^2.t
update.sigma2.t <- function(beta0, xi)
{
  ## get the inverse of SIGMA ------------------
  inv.sigma <- update.inv.sigma(inv.samp, xi)
  
  ### vectorize beta0
  beta0 <- beta0*ones
  ## calculate sigma.t square
  sigma2.t <- t((t.val-beta0)) %*%inv.sigma %*% (t.val-beta0)/n
  
  return(sigma2.t)
}




solve.equations <- function(t.val, samp.rho)
{
  n <- length(t.val)  # number of observations
  ones <- rep(1, n)  #  create a column of 1s
  ## initial value for beta0 (vector)
  beta0 <- mean(t.val)
  # starting value for xi
  ini.xi <- 0.5
  
  ## eigen decomposition of sample correlation matrix
  ###------- Calculated ONLY ONCE, used for all the steps that follows ---------
  inv.samp <- eigen(samp.rho)

  ### update beta0 first
  beta0.update <- update.beta0(xi = ini.xi)
  ### followed by sigma^2.t
  sigma2.t <- update.sigma2.t(beta0 = beta0.update, ini.xi)
  
  ## use optim() to get xi.values
  xi.update <- optim(ini.xi, obj.xi, method="Brent", lower=0, upper=1)

  return(list(beta0 = beta0.update, sigma2.t = sigma2.t, xi = xi.update$par))
}


solve.equations(t.val, samp.rho)
