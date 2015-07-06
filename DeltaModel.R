## simulate Duo's model

p <- 1 # proportion of DE
m <- 500 # number of genes
n <- 20 # biological sample in each group


## generate correlation matrix
rho <- 0.4
cor.struct <- matrix(rho, N, N);diag(cor.struct) <- 1
sigma <- 1.5*cor.struct
dat <- mvrnorm(1000, mu=rep(0, m), sigma)

sigma1 <- cor(dat1)
sigma2 <- 1.5*sigma1
dat <- mvrnorm(1000, mu=rep(0, m), sigma2)

sim.sigma <- cor(dat)


reps <- 100
tstat <- matrix(NA, nrow=reps, ncol=m)
for ( k in 1: reps)
{
  # whether a gene is DE or NOT ###
  z <- rbinom(m, size=1, prob=p)
  intens <- 1
  delta <- rexp(m, rate=intens)
  Delta <- z*delta
  
  x <- mvrnorm(n, mu=rep(0, m), sim.sigma)
  y <- mvrnorm(n, mu=Delta, sim.sigma)
  
  tstat[k,] <- apply(y, 2, mean)-apply(x, 2, mean) 
}


## if delta ~ exp(intens) and z~bernolli(p) are independent
## E(Delta) = intens*p, var(Delta)= intens^2*(2*p-p^2)
var.Delta <- intens^2*(2*p-p^2)


cov.tstat <- cov(tstat)

theoretical.cov <- sim.sigma*(1/n + 1/n)
diag(theoretical.cov) <- 1/n + 1/n + var.Delta

## see the first few elements of cov
e1 <- as.vector(cov.tstat)[1:10000]
e2 <- as.vector(theoretical.cov)[1:10000]
ids <- which(e2>0.2)
plot(e1[-ids], e2[-ids], pch=20)

plot(e1, e2, log="xy", pch=20)

range(diag(theoretical.cov))
diag.cov.t <- diag(cov.tstat)
hist(diag.cov.t)
