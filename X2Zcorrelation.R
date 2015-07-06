
## plot z.cor, sample.cor and true cor
plot.cor.z.test <- function(mu1, mu2, sig1, sig2, sig3, sig4, sigqua){
  
  
  true.cor  <- seq(-0.9, 0.9, 0.01) 
  z.cor <- samp.cor <- c()
  
  for ( j in 1: length(true.cor))
  {
    permu.B <- 100
    rho <- true.cor[j]
    
    # mean 
    #mu1 <- c(20, 30)
    #mu2 <- c(0, 32)
    
    # covariance 
    sigma1 <- matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2)
    sigma2 <- matrix(c(sig3^2, rho*sig3*sig4, rho*sig3*sig4, sig4^2), 2, 2)
    nsim1 <- 20
    nsim2 <- 30
    
    z.stat <- matrix(NA, permu.B, 2)
    s.cor <- c()
    for ( i in 1:permu.B)
    {
      
      x1 <- mvrnorm(nsim1, mu = mu1, sigma1)
      x2 <- mvrnorm(nsim2, mu = mu2, sigma2)
      
      s1 <- sqrt(sigma1[1,1]/nsim1 + sigma2[1, 1]/nsim2)
      s2 <- sqrt(sigma1[2,2]/nsim1 + sigma2[2, 2]/nsim2)
      
      t1 <- (mean(x1[, 1]) - mean(x2[, 1]))/s1
      t2 <- (mean(x1[, 2]) - mean(x2[, 2]))/s2
      
      
      
      rm1 <- x1[, 1] - mean(x1[, 1])
      rm2 <- x1[, 2] - mean(x1[, 2])
      
      rm3 <- x2[, 1] - mean(x2[, 1])
      rm4 <- x2[, 2] - mean(x2[, 2])
      
      ## I use mean of sample correlation as the sample correlation.
      rho1 <- cor(x1)[1,2]
      rho2 <- cor(x2)[1,2]
      s.cor[i] <- mean(c(rho1, rho2))
     # s.cor[i] <- cor(c(rm1, rm3), c(rm2, rm4)) # sample corr
      
      z.stat[i, ] <- c(t1, t2)
      
    }
    
    z.cor[j] <- cor(z.stat)[1, 2]
    samp.cor[j] <- mean(s.cor)
    
    
  }
  
  plot(true.cor, samp.cor, pch=3, cex=0.4, main=sigqua)
  abline(0, 1)
  points(true.cor, z.cor, pch=20, cex=0.5, col='red')
  legend("topleft", legend=c("sample", "test stat"), pch=c(3, 20), col=c("black", "red"))
}

mu1 <- c(20, 60)
mu2 <- c(5, 30)

plot.cor.z.test(mu1, mu2, sig1=11, sig2= 3,sig3=0.2, sig4=0.2, sigqua="DE")

mu1 <- c(20, 60)
mu2 <- c(200, 120)
plot.cor.z.test(mu1, mu2, sig1=100, sig2= 5,sig3=10, sig4=0.5, sigqua="rho constant")
plot.cor.z.test(mu1, mu2, sig1=10, sig2= 10, sig3=10, sig4=10, sigqua="Null") 
# for constant rho,  sig1*sig4 must equal to sig2*sig3
plot.cor.z.test(mu1, mu2, sig1=0.2, sig2= 0.3,sig3=0.2, sig4=0.3, sigqua="DE")
plot.cor.z.test(mu1, mu2, sig1=50, sig2= 0.3,sig3=0.03, sig4=10, sigqua="DE")

