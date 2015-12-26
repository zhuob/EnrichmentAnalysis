## solve the equation systems ##




#' @title Get the inverse of Sigma
#' 
#' @param sample.rho   sample correlation of genes
#' @param xi     what is the xi
#' @return the updated inverse sigma


#  t.val:  t value to be used in estimation
#  samp.rho: sample correlation matrix
#  start: starting value of optimize()
#  end:   end value of optimize()
#  n.or.nminus1: using biased or unbiased estimater of sigma

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
  #print(paste("beta=", beta.temp))
  
  ## caculate matrix C, defined in manuscript
  CC.temp <- BB - beta.temp*AA
  sigma2.t <- drop(CC.temp %*% inv.gamma %*% t(CC.temp)/n)      # update sigma2^2_T
  #print(sigma2.t)
  
  ## optimize objective function of xi
  a1 <- CC.temp %*% inv.gamma %*% DD %*% inv.gamma %*% t(CC.temp)       # term on the left part of equation
  a2 <- sum( diag(DD %*% inv.gamma ) )          # the trace, right part of the equation
  (a1/sigma2.t - a2)^2                          # return the objective function 
}



################  Original estimating equations (USE THIS ONE)---------------------------------
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
    #print(paste("beta=", beta.temp))
    
    ## caculate matrix C, defined in manuscript
    CC.temp <- BB - beta.temp*AA
    sigma2.t <- drop(CC.temp %*% inv.gamma %*% t(CC.temp)/n)      # update sigma2^2_T
    #print(sigma2.t)
    
    ## optimize objective function of xi
    a1 <- CC.temp %*% inv.gamma %*% DD %*% inv.gamma %*% t(CC.temp)       # term on the left part of equation
    a2 <- sum( diag(DD %*% inv.gamma ) )          # the trace, right part of the equation
    (a1/sigma2.t - a2)^2                          # return the objective function 
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




###################   Modified estimating equations ----------------------------


solve2.equations <- function(t.val, samp.rho, start=0, end=1, n.or.nminus1=n)
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


##################### Enrichment Score test ----------------------------


EnrichTest <- function(t.val, solve_eq,  samp.rho, go.term) # solve_eq is a object from solve.equations() function
{
  ones <- rep(1, length(t.val))
  xi <- solve_eq$xi
  sigma.t <- solve_eq$sigma2.t
  beta0 <- solve_eq$beta0
  
  inv.samp <- eigen(samp.rho)
  eig.value <- inv.samp$values                              # the original eigen values
  u <- inv.samp$vectors                                     # the original eigen vectors
  
  gamma <- diag(xi*eig.value +  1- xi)
  sigma <- u%*%gamma%*%t(u)
  sigma_inv=u%*%(diag(1/(xi*eig.value +  1- xi)))%*%t(u)    # the inverse of sigma
  
  p1 <- 1/sigma.t* (t(go.term) %*% sigma_inv %*% (t.val - beta0*ones))^2  # numerator
  p2 <- t(go.term) %*% sigma_inv %*% go.term
  p3 <- (t(go.term) %*% sigma_inv %*% ones)^2 / (t(ones) %*% sigma_inv %*% ones)
  
  test.stat <- p1/(p2- p3)                                  # the test statistic
  
  pval <-  1 - pchisq(test.stat, 1)
  c(test.stat, pval)
  return(list(beta=beta0, sigma.t= sigma.t, xi=xi,
              statistic = test.stat, p = pval))
}







EnrichmentTest <- function(microarray, trt, index)
{
  size <- ncol(microarray)
  go_term <- index
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  
  samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
  t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment
  
  
  ######## begin the test procedure ---------------------------------------
  solve_eq <- solve.equations(t_val, samp_rho)             # estimate the xi, beta and sigma
  obj0 <- EnrichTest(t_val, solve_eq, samp_rho, go_term)  # use the sample correlation matrix
  #  pval3 <- summary(lm(t_val ~ go_term))$coe[2, 4]             # simple linear regression
  ## the reason why this does not work is because the t_val is not normal any more. 
  return(obj0)
}

