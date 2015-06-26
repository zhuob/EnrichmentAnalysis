##  the basic setup is that I simulate a dataset, and permute the columns to
##  nullify the treatment effect. Then for each permutation, I calculate a 
## significant statistic for the coefficient beta1, and see the type I error rate.



## data simulation ####

source("SimulateLabData.R")

size <- 1000
N <- 500
rho <- 0.2
prob <- 0.1
DE <- N*prob

mu1 <- rep(0, N)
mu2 <- c(rep(5, DE)[1:DE], rep(0, N-DE))

dat <- simu.microarray(size = size, N = N,prob = prob, rho=rho, mu1=mu1, mu2=mu2)




### ####  calculate sample correlations  and fix it when permuting columns


new.dat <- permute.dat(dat$data)

mean(cor(t(new.dat)))
mean(cor(t(dat$data)))



## do 1000 enrichment tests

nsim <- 100
enrichment.sim <- matrix(NA, nrow=nsim, ncol= 7)
library(regress)

for ( i in 1: nsim)
{
  
  # a size of 100 genes are in the GO term
  GO <- rep(0, N)
  samp.GO <- sample(1:N, size)
  GO[samp.GO] <- 1
  
  # permute the data
  new.dat <- permute.dat(dat$data)
  # calculate correlation matrix
  cor.permute <- calculate.cor(new.dat)
  
  # calculate the test stat
  T1 <- apply(new.dat, 1, t.result)
  id <- which(T1< qt(0.05, size-2))
  T1[id] <- -1*T1[id]
  
  print(i)
  # model fitting
  model.sim2 <- regress(T1~GO, ~cor.permute)
  mod.lm <- lm(T1~GO)
  p.beta <- summary(mod.lm)$coef[2, 4]
  value2 <- c(model.sim2$beta[1], model.sim2$beta.se[1], 
              model.sim2$beta[2], model.sim2$beta.se[2],
              model.sim2$sigma[1], model.sim2$sigma[2])
  
  enrichment.sim[i, ] <- c(value2, p.beta)
  
}

p <- pt(enrichment.sim[, 3]/enrichment.sim[, 4], 98)
  
## error: Sigma is not positive definite on contrasts: range(eig)= -2.46726e-16 1.198346

# https://github.com/cran/regress/blob/master/R/regress.R




############## #########  
###  SIMULATE THE NULL CASE EXPRESSION DATA
source("SimulateLabData.R")

lm.pval <- function(rho, prob)
{
  
  size <- 100
  N <- 500
 # rho <- 0.4
# prob <- 0.0
  DE <- N*prob
  
  GO <- rep(0, N)
  samp1 <- sample(1:DE,  0.4*DE)
  samp2 <- sample((DE+1):N, 100-0.4*DE)
  samp.index <- c(samp1, samp2)
  GO[samp.index] <- 1
  

  mu1 <-  rep(0, N)
  mu2 <-  rep(0, N)
  if (DE > 0 ){
    mu2 <- c(rep(5, DE)[1:DE], rep(0, N-DE))
  }
  
  
  cor.struct <- matrix(rho, N, N);diag(cor.struct) <- 1
  sigma <- 1.5*cor.struct
  
  dat <- simu.microarray(size = size, N = N,prob = prob,
                         mu1=mu1, mu2=mu2, sigma=sigma)
 # mean(cor(t(dat$data)))
 #  mean(calculate.cor(dat$data))
  
  
  #  Step 1: sample the GO term randomly
  #  Step 2: calculate p value using regress()
  #  Step 3: calculate p value by LM
  #  Step 4: permuate columns to get estimate for beta1, and therefore p value
  
  nsim <- 1000
  p.val <- c()
  
  for ( i in 1:nsim)
  {
    
    dat1 <- permute.dat(dat$data)
    # GO <- rep(0, N)
    # ids <- sample(1:N, 100)
    # GO[ids] <- 1
    y <- apply(dat1, 1, t.result)
    #   y <- dat$data[, 1]
   # p.t <- pt(y, 98)
  #  hist(p.t, breaks=30)
    
    mod <- summary(lm(y~GO))
    
    p.val[i] <- mod$coefficients[2, 4]
  }
 return(p.val) 
}

p0 <- lm.pval(rho=0, prob=0)
p10 <- lm.pval(rho=0, prob=0.1)
p11 <- lm.pval(rho=0, prob=0.2)

p1 <- lm.pval(rho=0.2, prob=0)
p2 <- lm.pval(rho=0.4, prob=0)
p3 <- lm.pval(rho=0.6, prob=0)

p4 <- lm.pval(rho=0.2, prob=0.1)
p5 <- lm.pval(rho=0.4, prob=0.1)
p6 <- lm.pval(rho=0.6, prob=0.1)

p7 <- lm.pval(rho=0.2, prob=0.2)
p8 <- lm.pval(rho=0.4, prob=0.2)
p9 <- lm.pval(rho=0.6, prob=0.2)



permute.simu.p <- data.frame(p0, p10, p11, p1, p2, p3, p4, p5, p6,
                             p7, p8, p9)
saveRDS(permute.simu.p,"/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/EnrichmentTypeIerrorLinearRegSimulation.rds")

permute.simu.p <- readRDS("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/EnrichmentTypeIerrorLinearRegSimulation.rds")


## p1  rho =0.2, NO DE
## p2  rho = 0.2, 10% DE
## p3  rho = 0.6, 10% DE
## p4  rho = 0.6, 20% DE
## p5  rho = 0.6, 40% DE

## distribution of p is not related to the proportion of DE based on my simulation setup


rho.vec <- c(rep(0, 3), rep(c(0.2, 0.4, 0.6), 3) )
prob.vec <- c(c(0, 0.1, 0.2), rep(c(0, 0.1, 0.2), each=3))

par(mfrow=c(3, 2))
 for ( i in 1:dim(permute.simu.p)[2])
 {
   hist(permute.simu.p[, i], breaks=20, xlab="p",
        main=paste("rho=", rho.vec[i], ", DE = ", prob.vec[i]))
   
 }


dev.off()



p.est.cor <- readRDS("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/enrichmentTypeIerrorSim.est.cov.rds")

p.true.cor <- readRDS("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/enrichmentTypeIerrorSim.true.cov.rds")

p.est <- pnorm(p.est.cor$go/p.est.cor$go.sd)

hist(p.est)

p.true <- pt(p.true.cor$go/p.true.cor$go.sd, 98)
hist(p.true)



pvals0 <- readRDS("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/enrichmentTypeIerror.pvalue.rds")
pvals1 <- readRDS("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/enrichmentTypeIerror.pvalue.4.1.rds")
pvals2 <- readRDS("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/enrichmentTypeIerror.pvalue.4.2.rds")



par(mfrow=c(3, 2))
xaxis <- c("lm", "sample cor", "true cor")
for ( i in 1:3)
{
  hist(pvals0[, i], breaks=20, xlab=xaxis[i],
       main=paste("rho= 0.4", ", DE = 0%"))
}

for ( i in 1:3)
{
  hist(pvals1[, i], breaks=20,  xlab=xaxis[i],
       main=paste("rho= 0.4", ", DE = 10%"))
}

for ( i in 1:3)
{
  hist(pvals2[, i], breaks=20, xlab=xaxis[i],
       main=paste("rho= 0.4", ", DE = 20%"))
}

 






