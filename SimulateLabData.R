#'   run poisson regression with random effect , and keep track of warnings 
#'   or errors of regression for a single gene
#' 
#' 
#' @title Simulate normally distributed expression data 
#'
#' @param size  number of biological samples
#' @param N     number of genes
#' @param prob  DE probability
#' @param rho   correlation coefficient for exchangeable covariance
#' @param mu1   a 1 by N vector for mean of trt 1
#' @param mu2   a 1 by N vector for mean of trt 2
#' @return a list
#' \item{data} simulated data with column for sample, row for genes
#' \item{mu1} group mean from which the treatment 1 is simulated.
#' \item{mu2} group mean from which the treatment 2 is simulated
#'


library(MASS)

simu.microarray <- function(size, N=50, prob=0.1, rho=0.4, mu1, mu2, sigma)
{
  set.seed(100)
  # mu1 <- rep(0, N) # mean for treatment group 1
  DE <- N*prob
  # make the first DE genes shift by 5 and -5 alternatively
  # and others are not DE
  
  # mean for treatment group 2
  # mu2 <- c(rep(c(5,-5), DE)[1:DE], rep(0, N-DE))
  
  microarray <- matrix(nrow=N, ncol=size)
  n <- size/2
  
  cor.struct <- matrix(rho, N, N);diag(cor.struct) <- 1
 
 #  sigma <- 1.5*cor.struct # + 0.5*diag(N)
  
  # each row represents one sample from MVN, need to transpose
  set.seed(50)
  microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
  set.seed(51)
  microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))

  ls <- list(data=microarray, mu1=mu1, mu2=mu2)
  return(ls)

}






##  first remove the treatment effects and get the correlation matrix
calculate.cor <- function(data)
{
  size <- dim(data)[2]
  n <- size/2
  trt <- rep(c(1, 2), each=n)
  
  demean <- matrix(NA, nrow = dim(data)[1], ncol=dim(data)[2])
  mu1 <- apply(data[, 1:n], 1, mean)
  mu2 <- apply(data[, -(1:n)], 1, mean)
  mu1.new <- matrix(rep(mu1, each=n), ncol= n, byrow=T)
  mu2.new <- matrix(rep(mu2, each=n), ncol= n, byrow=T)
  
  mu.mat <- cbind(mu1.new, mu2.new)
  demean <- data - mu.mat
  
  return(cor(t(demean)))
}




###  calculate T statisitcs 
t.result <- function(data, quantity= "stat")
{
  n <- length(data)
  t1 <- data[1:n/2]
  t2 <- data[-(1:n/2)]
  t.test(t1, t2)$stat
}



### permute function for the exression matrix

permute.dat <- function(data)
{
  numcol <- dim(data)[2]
  data.new <- matrix(NA, nrow = dim(data)[1], ncol=dim(data)[2])
  size <- numcol/2
  id1 <- sample(1:numcol, size)
  data.new[, 1:size] <- data[,id1]
  data.new[, -(1:size)] <- data[, -id1]
  return(data.new)
}





# a function to calculate group mean
group.mean <- function(data, group)
{
  
  id1 <- which(group == 1)
  nrep <- length(id1)
  mean1 <- apply(data[, id1], 1, mean)
  mean1 <- matrix(rep(mean1, nrep), ncol=nrep, byrow=F) 
  mean2 <- apply(data[, -id1], 1, mean)
  mean2 <- matrix(rep(mean2, nrep), ncol=nrep, byrow=F) 
  
  data.mean <- data.frame(mean1, mean2)
  colnames(data.mean) <- colnames(data)
  rownames(data.mean) <- rownames(data)
  data.mean
}
