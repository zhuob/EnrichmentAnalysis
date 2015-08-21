

#####################  a simulation study  -----------------------------
source("solve.R")

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




go.term <- ifelse(t.val < -4 | t.val > 2, 1, 0)

mean(t.val[go.term==1])
mean(t.val[go.term==0])

EnrichTest(t.val, r1, go.term)



# ## convergence == 0 does not mean it's converged. 
# fy <- function(y){ (y^2- 3*y + 2)^2}
# optim(0.6, fy,  method="Brent", lower=0, upper=0.8)
# optimize(fy, c(0, 1))

BetaSigmaXi <- read.table("BetaSigmaXi4.txt", header=T)
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

sigma <- seq(0.05, 0.95, by = 0.01)
plot(sigma, BetaSigmaXi$xi, pch=20, cex=0.5)
abline(0 ,1)
