setwd("~/Google Drive/Study/Graduation/Correlation/Data")

p.mat <- readRDS("p.NBP.21labs.rds")
arab.rm.mean <- readRDS("arab.rm.mean.by.group.rds")


# 
# sta <- read.csv("structuralMoleculeActivity.csv", header=T)
# sta <- unique(sta[-dim(sta)[1], 2])
# 
# eta <- read.table("ElectronTransport EnergyPathways.txt", header=F)
# # this one might be helpful
#  DNA <- read.table("DNAmetabolism.txt", header=F)
# signal <- read.table("signalTransduction.txt", header=F)
# proteinMetabolism <- read.table("ProteinMetabolism.txt", header=F)
# transport <- read.table("transport.txt", header=F)
# 
#  GOterm <- as.vector(DNA)
ribosome <- read.csv("GOribosome.csv", header=F)
GOterm <- unique(ribosome[-dim(ribosome)[1], 2])

## get the p values for the first experiment, and estimate COV from the same exp.
id.gene <- which(p.mat$Gene %in% GOterm)

## genes in GO terms
id.gene2 <- which(row.names(arab.rm.mean) %in% GOterm)
cor.mat.1 <- cor(t(arab.rm.mean[id.gene2, 1:6]))
p.1 <-  p.mat[id.gene, 2]


library(regress)
model1 <- regress(p.1 ~ 1, ~cor.mat.1)


## logistic regression 
p.all <- p.mat[, 2]
log.p <- -log(p.all)
respon <- rep(0, dim(p.mat)[1])
respon[id.gene] <- 1
summary(lm(p.all~respon))


log.p <- -log(p.all)
log.log.p <- log(1-log(p.all))
# respon <- ifelse(p.mat$Gene %in% GOterm[, 1], 1, 0)


summary(lm(log.p~respon))
model2 <- regress(-log(p.1) ~ 1, ~cor.mat.1)
summary(lm(log.log.p~respon))
model3 <- regress(log(1-log(p.1)) ~ 1, ~cor.mat.1)

c(mean(p.all[respon==0]), mean(p.all[respon==1]))
summary(model1)

c(mean(log.p[respon==0]), mean(log.p[respon==1]))

# fisher's exact test

p.cutoff <- function(alpha)
{
  dichot <- ifelse(p.all < alpha, 1, 0)
  chi.table <- table(data.frame(dichot, respon))
  chi.table
  fisher.test(chi.table)$p.value
}

alpha <- seq(0.01, 0.10, by=0.001)
fisher <- c()
for ( i in 1: length(alpha))
{
  fisher[i]  <- p.cutoff(alpha[i])
  
}

plot(alpha, fisher, type="l", main="cutoff vs p value for fisher exact test")


## some simulation to
nsim <- 100
sigma1 <- 0.048
sigma2 <- 0.056
n <- dim(cor.mat.1)[1]

beta <- s1 <- s2 <- samp.mean <- rep(0, nsim)
for ( i in 1:nsim)
{p.sim <- mvrnorm(n=1, mu = rep(0.426, n), Sigma = sigma1*cor.mat.1) + rnorm(n, 0, sigma2)
model.sim <- regress(p.sim ~ 1, ~cor.mat.1)   
  beta[i] <- as.numeric(model.sim$beta)
  s1[i] <- as.numeric(model.sim$sigma[1])
  s2[i] <- as.numeric(model.sim$sigma[2])
  samp.mean[i] <- mean(p.sim)

}










# 
# 
# Fatal error: cannot open file 'regressP.r': No such file or directory
# Likelihood kernel: K = (Intercept)
# 
# Maximized log likelihood with kernel K is  410.597 
# 
# Linear Coefficients:
#   Estimate Std. Error
# (Intercept)    0.426      0.054
# 
# Variance Coefficients:
#   Estimate Std. Error
# cor.mat.1    0.048      0.036
# In           0.056      0.004
# 
# Likelihood kernel: K = (Intercept)
# 
# Maximized log likelihood with kernel K is  10389.15 
# 
# Linear Coefficients:
#   Estimate Std. Error
# (Intercept)    0.388       0.01
# 
# Variance Coefficients:
#              Estimate Std. Error
# cor.mat.0    0.001      0.001
# In           0.099      0.001



# 
# Call:
#   glm(formula = respon ~ log.p, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.8251  -0.6983  -0.6975  -0.6972   1.7514  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -1.290724   0.020571 -62.746   <2e-16 ***
#   log.p        0.002004   0.003283   0.611    0.542    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 17026  on 16292  degrees of freedom
# Residual deviance: 17025  on 16291  degrees of freedom
# AIC: 17029
# 
# Number of Fisher Scoring iterations: 4
# 
# respon
# dichot     0     1
# 0 10338  2752
# 1  2427   776
# 
# 
# 
# Fisher's Exact Test for Count Data
# 
# data:  dichot.table
# p-value = 9.535e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.094969 1.316652
# sample estimates:
# odds ratio 
#   1.201047 
# 
# Likelihood kernel: K = (Intercept)
# 
# Maximized log likelihood with kernel K is  10637.17 
# 
# Linear Coefficients:
#              Estimate Std. Error
#  (Intercept)    0.393       0.01
# 
# Variance Coefficients:
#              Estimate Std. Error
#      cor.mat    0.001      0.001
#      In         0.100      0.001



##




# sample in the pool of non DE genes
threshold <- 0.1
candidate <- which(apply(p.mat[, -1], 1, min) > threshold)
arab.rm.mean.nonDE <- arab.rm.mean[candidate, ]
p.val.nonDE <- p.mat[candidate, ]
z.val <- qnorm(p.val.nonDE[, 2])

cor.mat.nonDE <- cor(t(arab.rm.mean.nonDE))

library(regress)
model2 <- regress(z.val ~ 1, ~cor.mat.nonDE)
model3 <- regress(p.val.nonDE[, 2] ~ 1, ~cor.mat.nonDE)

model2
model3

hist(cor.mat.nonDE)
hist(z.val)

## with this covariance matrix, I simulate the responses 

library(regress)
library(MASS)
subset <- dim(cor.mat.nonDE)[1]

N <- 200
mu <- c(rep(-.2, N), rep(0, subset-N))
cor.struct <- cor.mat.nonDE
# cor.struct <- matrix(0.7, subset, subset);diag(cor.struct) <- 1

sigma <-  1*cor.struct + diag(1, subset)
y <- mvrnorm(1, mu, sigma) 

x <- c(rep(0, N), rep(1, subset-N))

model.sim <- lm(y~x)
model.sim2 <- regress(y~x, ~cor.struct)
summary(model.sim)
model.sim2


X <- cbind(rep(1, subset), x)
sig.inv <- solve(t(X)%*% solve(sigma) %*% X)
sqrt(sig.inv[1,1])*model.sim2$sigma[1]
sqrt(sig.inv[2,2])*model.sim2$sigma[2]

A <- eigen(cor.struct)


#############################################################
#             log transform of raw counts
#############################################################

library(DESeq2)

arab <- readRDS("arabidopsis.21.rds")
trt <- readRDS("treatment.21.rds")

head(p.mat[, 2])
arab1 <- arab[[1]]
idx <- which(row.names(arab1) %in% p.mat$Gene)
arab1 <- as.matrix(arab1[idx, ])
arab.log <- rlog(arab1, blind=F)


candidate <- sample(1:dim(arab.log)[1], 500)
arab.log2 <- arab.log[candidate, ]
arab.log2[1:3, ] <- arab.log2[1:3, ] - apply(arab.log2[1:3,], 1, mean)
arab.log2[-(1:3), ] <- arab.log2[-(1:3), ] - apply(arab.log2[-(1:3),], 1, mean)

cor.mat2 <- cor(t(arab.log2))


dists <- dist(t(arab.log))
plot(hclust(dists))

DEtest <- function(x){
  x1 <- x[1:3]
  x2 <- x[-(1:3)]
  test.p <- t.test(x1, x2)$p.value
  return(test.p)
}

p.two.sample <- apply(arab.log, 1, DEtest)

cor(p.two.sample, p.mat[, 2])
hist(p.two.sample)
hist(p.mat[, 2])



model4 <- regress(p.two.sample[candidate] ~ 1, ~cor.mat2)
model4

z.t.test <- qnorm(p.two.sample[candidate])

model5 <- regress(z.t.test ~ 1, ~cor.mat2)
model5

cor(z.val, z.t.test)
cor(p.two.sample, p.mat[, 2])






### from z to p, correlation can be preserved if z is standard normal
rho <- seq(-0.9, 0.9, by=0.01)
mu <- c(0, 0.3)

cor.p2z <- matrix(NA, length(rho), 2)

for ( k in 1:length(rho))
{
  r <- rho[k]
  sigma <- matrix(c(1, r, r, 1), 2,2)
  rand1 <- mvrnorm(100, mu, sigma)
  z.cor <- cor(rand1)[1,2]
  p <- pnorm(rand1)
  p.cor <- cor(p)[1, 2]
  cor.p2z[k, ] <- c(z.cor, p.cor)
}

plot(rho, cor.p2z[, 1], ylab="sample", pch=3, cex=0.5)
points(rho, cor.p2z[, 2], pch=20, cex=0.5, col="blue")

z1 <- qnorm(p.two.sample)
c(mean(z1), var(z1))

z2 <- qnorm(p.mat[, 6])
c(mean(z2), var(z2))


