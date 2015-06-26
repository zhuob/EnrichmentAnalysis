#  poisson regression statistic correlation VS sample correlation
#

N <- 100
mu <- rep(c(5,10), each=N/2)
x <- rep(c(1,0), each=N/2)

y <- rpois(N, mu)

pois.mod <- glm(y~x, family = poisson(link = "log"))
beta0 <- summary(pois.mod)$coe[1, 1]
beta1 <- summary(pois.mod)$coe[2, 1]

### calculate the hessian matrix 
var.beta1 <- sum(exp(beta0 + beta1*x)*x^2)
var.beta0 <- sum(exp(beta0 + beta1*x))
cov.beta0beta1 <- sum(exp(beta0+ beta1*x)*x)

hessian <- matrix(c(var.beta0, cov.beta0beta1, cov.beta0beta1, var.beta1), 2, 2)
## calculate the standard deviation for each parameter. 
## It's the same as the output of glm function.
sqrt(diag(solve(hessian)))
  
cor.pois <- function(n, lambda1, r)  ## generate correlated poisson random variables
{
  
lambda2 <- (1-r^2)/(0.1+r^2)*lambda1 # this makes lambda3 always positive
lambda3 <- lambda1^2/(r^2*(lambda1+lambda2))-lambda1

cor <- lambda1/(sqrt((lambda1+lambda3)*(lambda1+lambda2)))

  y1 <- rpois(n, lambda1)
  y2 <- rpois(n, lambda2)
  y3 <- rpois(n, lambda3)
  
  x1 <- y1  + y2
  x2 <- y1  + y3
  
  # the correlation between x1 and x2 
  # rho <- mu1 /sqrt( (mu1 + mu3)*(mu1 + mu2))
  #  print(round(rho, 4))
  return(cbind(x1, x2))
#  return(c(lambda1, lambda2, lambda3, r))
}



num <- 20  # how many samples in each group
# in this case, gene 1 is not DE, but gene 2 is DE
# for treatment 1

source("SimulateLabData.R")

nsim <- 1000
test.correlation <- sample.correlation <- c()
r <- seq(0.1, 0.95, length=50)
trt <- rep(c(0, 1), each=num)


calculate.lambda <- function(lambda1, r)
{
  lambda2 <- (1-r^2)/(0.1+r^2)*lambda1 # this makes lambda3 always positive
  lambda3 <- lambda1^2/(r^2*(lambda1+lambda2))-lambda1
  return(c(lambda1, lambda2, lambda3, r))
}

lambda.20 <- matrix(NA, ncol=4, nrow=length(r))
for ( j in 1:length(r))
{
  lambda.20[j, ] <- calculate.lambda(37, r[j])
  
}

for( k in 1: length(r))
{
  
  t1 <- t2 <- pois.cor <-  c()
  
for ( i in 1:nsim)
{
  
  exp1 <- cor.pois(num, 100, r[k]) # pi1 = 4, p2= 5, rho = 0.2236
  # for treatment 2
  exp2 <- cor.pois(num, 37, r[k]) #pi1 = 4, p2= 20 rho = 0.2236
  
  poisdata <- rbind(exp1, exp2)
  pois.cor[i] <- calculate.cor(t(poisdata))[1, 2]
  
  trt <- rep(c(0, 1), each=num)
  pois1 <- glm(poisdata[, 1]~trt, family = poisson(link = "log"))
  #pois1 <- lm(poisdata[, 1]~trt)
  t1[i] <- summary(pois1)$coe[2, 3] # wald statistic
  
  # pois2 <- lm(poisdata[, 2]~trt)
  pois2 <- glm(poisdata[, 2]~trt, family = poisson(link = "log"))
  t2[i] <- summary(pois2)$coe[2, 3] # wald statistic
  
}
  
sample.correlation[k] <- mean(pois.cor)
test.correlation[k] <- cor(t1, t2)
}

plot(sample.correlation, r, pch=20)
points(sample.correlation, test.correlation, col="red", pch=3)

pois.wald.stat.sample.cor <- data.frame(
  wald.stat= test.correlation, samp.cor= sample.correlation, true=r
)

saveRDS(pois.wald.stat.sample.cor, "pois.regression.wald.stat.sample.cor.rds")




############ DOES regression have the power to remove DE effect?
###    NO!  It's the same as two sample t-test
## gene1 not DE, gene2 DE
mu1 <- c(5, 3)  ## gene1 gene2 groupA
mu2 <- c(20, 18) ## gene1 gene2 groupB

## normally correlated random variables

nsim <- 1000
num <- 20
test.correlation <- sample.correlation <- c()
r <- seq(0.1, 0.95, length=50)
trt <- rep(c(0, 1), each=num)

for ( k in 1:length(r))
{
  t1 <- t2 <- pois.cor <-  c()
  
  sigma <- matrix(c(1.5^2, r[k]*1.5*2, r[k]*1.5*2, 2^2), 2, 2)
  for ( i in 1:nsim)
  {
    
    
    exp1 <- mvrnorm(num, mu1, sigma) 
    # for treatment 2
    exp2 <- mvrnorm(num, mu2, sigma) #pi1 = 4, p2= 20 rho = 0.2236
    
    poisdata <- rbind(exp1, exp2)
    pois.cor[i] <- calculate.cor(t(poisdata))[1, 2]
    
    trt <- rep(c(0, 1), each=num)
    #pois1 <- glm(poisdata[, 1]~trt, family = poisson(link = "log"))
    pois1 <- lm(poisdata[, 1]~trt)
    t1[i] <- summary(pois1)$coe[2, 3] # wald statistic
    
    pois2 <- lm(poisdata[, 2]~trt)
    #pois2 <- glm(poisdata[, 2]~trt, family = poisson(link = "log"))
    t2[i] <- summary(pois2)$coe[2, 3] # wald statistic
    
  }
  
  sample.correlation[k] <- mean(pois.cor)
  test.correlation[k] <- cor(t1, t2)
  
}

plot(sample.correlation,r,  pch=20, main="sample.size=100")
points(sample.correlation, test.correlation,, col="red", pch=3)






