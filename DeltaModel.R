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



## # use a group of real data to estimate 

GSE38400 <- read.table("/Users/Bin/Google Drive/Study/Thesis/NCBI/Data/GSE38400.Rsubread.txt", header=T)
row.mean <- apply(GSE38400, 1, mean)
select.gene <- which(row.mean > 5)
GSE38400 <- GSE38400[select.gene, ]


##  select genes
n.gene <- dim(GSE38400)[1]
set.seed(500)
sample.gene <- sample(1:n.gene, 500)
sample.data <- GSE38400[sample.gene, ]

# calculate sample correlation
source("SimulateLabData.R")

# calculate group mean
d1 <- group.mean(sample.data[, 1:6], group=c(1,1, 1, 2, 2, 2))
d2 <- group.mean(sample.data[, -(1:6)], group=c(1,1, 1, 2, 2, 2))
d.mean <- cbind(d1, d2)

# sample correlation
data.rm.mean <- sample.data-d.mean
sample.correlation <- cor(t(data.rm.mean))


## sample GO term
set.seed(100)
GO.id <- sample(1:500, 100)
GO <- rep(0, 500)
GO[GO.id] <- 1

## calculate test statistics
t.stat <- apply(sample.data, 1, t.result)
## make the negative t.stat postive
ids <- which(t.stat<0)
t.stat[ids] <- t.stat[ids]*(-1)


## fit by regress function
v1 <- v2 <- matrix(0, 500, 500)
diag(v1) <- GO
diag(v2) <- -1*(GO-1)

library(regress)
model1 <- regress(t.stat~ GO, ~ v1 + v2 + sample.correlation, identity = F)
model2 <- regress(t.stat~ GO, ~sample.correlation, identity = T)


##  instead of randomly selecting the GO term, I choose genes that have large t.stat 

DE.genes <- which(t.stat > 3)     
set.seed(123)
GO.id <- sample(DE.genes, 100)
GO <- rep(0, 500)
GO[GO.id] <- 1

v1 <- v2 <- matrix(0, 500, 500)
diag(v1) <- GO
diag(v2) <- -1*(GO-1)

model3 <- regress(t.stat~ GO, ~ v1 + v2 + sample.correlation, identity = F)
model4 <- regress(t.stat~ GO, ~sample.correlation, identity = T)



##   A simulation study ##
a1 <- model3$sigma
 sigma <- a1[1]*v1 + a1[2]*v2 + a1[3]*sample.correlation

 m <- 500
 p <- 0.2
 z <- rbinom(m, size=1, prob=p)
 intens <- 2
 delta <- rexp(m, rate=intens)
 Delta <- z*delta
 
 n <- 30

 x <- mvrnorm(n, mu=rep(0, m), sigma)
 y <- mvrnorm(n, mu=Delta, sigma)
 simu.data <-  t(rbind(x,y))
simu.tstat <- apply(simu.data, 1, t.result)

set.seed(100)
GO.id <- sample(1:500, 100)
GO <- rep(0, 500)
GO[GO.id] <- 1

v1 <- v2 <- matrix(0, 500, 500)
diag(v1) <- GO
diag(v2) <- -1*(GO-1)

model5 <- regress(simu.tstat~ GO, ~ v1 + v2 + sample.correlation, identity = F)
model6 <- regress(simu.tstat~ GO, ~sample.correlation, identity = T)


