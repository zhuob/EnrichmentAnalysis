

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
    
    answer <- solve.equations(t.val, r1)
    return(unlist(answer))
  }
)



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







###########################  Real Case Study -------------------------------------
arab <- readRDS("~/Google Drive/Study/Thesis/Correlation/Data/arabidopsis.21.rds")
trt <- readRDS("~/Google Drive/Study/Thesis/Correlation/Data/treatment.21.rds")

source("SimulateLabData.R")

arab_1 <- scale(arab[[1]])                            # standardize by column
group_mean <- group.mean(arab_1, trt[[1]])
resid_mat <- arab_1 - group_mean

for (i in 2:length(trt) ) {
  
  arab_i <- scale(arab[[i]])                          # STANDARDIZE BY COLUMN
  group_mean <- group.mean(arab_i, trt[[i]])          # calculate the group means within each data set
  resid_temp <- arab_i - group_mean                   # the residual after the mean is removed
  resid_mat <- merge(resid_mat, resid_temp, by=0)     # merge the residuals
  rname <- resid_mat[, 1]                             # keep the row names
  resid_mat <- resid_mat[, -1]
  row.names(resid_mat) <- rname                       # give the rows names
}




#########    prepare the go terms -------------------------------

AllGoterm <- read.table("ATH_GO_GOSLIM.txt", sep="\t", header=F, fill=T) %>%
  filter(grepl('GO:', V6))

GotermList <- group_by(AllGoterm, V6) %>%              # find the go terms
  summarise(NumberofGene=n()) %>%        # count how many genes in each go term
  filter(NumberofGene >= 100)            # select Go terms that have at least 200 genes

GotermList <- as.matrix(GotermList)




######## consider only the genes present in the go term lists
go_term_genes <- unique(AllGoterm$V1)
ids <- which(row.names(arab[[1]]) %in% go_term_genes) # get 13672 genes in total



#####  get a random sample from the Genes, and match them with a GO term.

n_gene <- 1000                                         # number of genes to be considered
# ids <- which(rowMeans(arab[[1]]) > 3)                 # keep the genes that have at least 3 counts for each sample
set.seed(200)                                          # this gene list can be repeated
get_id <- sample(ids, n_gene)
total_gene <- scale(arab[[1]][get_id, ])                     # sample n_gene from the list
sample_size <- dim(total_gene)[2]

samp_rho <- cor(t(resid_mat[get_id, ]))               # the sample correlations


####  Loop through each of the go terms to find the interesting ones ------------------
choose_Go_term <- matrix(NA, nrow=dim(GotermList)[1], 2)
for ( k in 1:dim(GotermList)[1])
{
  go <- filter(AllGoterm, V6 == GotermList[k, 1])$V1             # the members in the go term
  go.in <- which(noquote(row.names(total_gene)) %in% go)
  choose_Go_term[k, ] <- c(as.numeric(GotermList[k, 2]), length(go.in))
  
}
ids2 <- which(choose_Go_term[, 2] >= n_gene*0.01)

go <- filter(AllGoterm, V6 == GotermList[ids2[1], 1])$V1             # the members in the go term


sample_size <- dim(total_gene)[2]                         # sample size
go.in <- which(noquote(row.names(total_gene)) %in% go)    # which Go term genes are present in this list
go_term <- rep(0, n_gene)
go_term[go.in] <- 1
# go_term <- rbinom(n_gene, 1, 0.2)

x_bar <- apply(total_gene[, 1:sample_size/2], 1, mean)
y_bar <- apply(total_gene[, -(1:sample_size/2)], 1, mean)

t_val <- (x_bar-y_bar)# /apply(total_gene, 1, sd)           # mean of treat - mean of control

t1 <- mean(t_val[go_term==1])                               # mean diff within go term
t2 <- mean(t_val[go_term==0])                               # mean diff outside go term
go_num <- length(go.in)                                     # number of genes from go term


xlable <- 1:n_gene
plot(t_val[-go.in], xlable[-go.in], pch=20, cex=0.5)
points(t_val[go.in], xlable[go.in], pch=3, col="red")


# go_term <- ifelse(t_val<0, 1, 0)

result <- EnrichTest(t_val, samp_rho, go_term)              # enrichment test
result




############## Simulation Study -------------------------------------------
source("SimulateLabData.R")
source("Solve.R")

p.mat <- matrix(NA, nrow=100, ncol=2)

for (p in 1: 100) {
  
size <- 50                                     # number of samples to be considered
n_gene <- 500                                  # number of genes to be considered
de_prop <- 0.2                                 # proportion of DE
de <- n_gene*de_prop                           # number of DE genes
delta <- rexp(n_gene, 0.2)                     # magnitude of DE
z <- rbinom(n_gene, 1, de_prop)                # indicator, DE =1, NON_DE=0 
delta <- delta*z
delta <- rep(0, n_gene)             ## NULL, NO TREATMENT effect

rho <- 0.5                                     # correlation matrix
sigma0 <- 2                                    # variance part 
cor.struct <- matrix(rho, n_gene, n_gene);    
diag(cor.struct) <- 1                          # correlation matrix
sigma <- sigma0*cor.struct                     # covariance matrix for both groups

mu_initial <- 0
mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
mu2 <- delta  + mu_initial                     # mean expression value for treatment


microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
n <- size/2
#set.seed(50)
microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
#set.seed(51)
microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))

trt <- rep(c(0, 1), each= size/2)             # group indicator


# prepare the t_stat, samp_rho and Go term 


group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix

samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment


# create a go term where more de genes are contained
n_go <- 0.2*n_gene                                       # number of genes in the go term
go_term <- rep(0, n_gene)                                # go term indicator
go_prop <- 0.4                                           # proportion of DE genes in the go term 
go_de <- sample(which(z==1), n_go*go_prop)               # de genes in go term 
go_not_de <- sample(which(z==0), n_go*(1-go_prop))       # non DE genes in go term 
go_in <- c(go_de, go_not_de)
go_term[go_in] <- 1                                      # 50% go term genes are de
#go_term <- rep(0, length(t_val))
#ind <- which(t_val< -2)
#go_term[ind] <- 1
 

mean(t_val[go_term==1])                                   # mean diff within go term
mean(t_val[go_term==0])                                   # mean diff outside go term
length(go_in)

# xlable <- 1:n_gene
# plot(t_val[-go_in], xlable[-go_in], pch=20, cex=0.5)
# points(t_val[go_in], xlable[go_in], pch=3, col="red")


pval1 <-  EnrichTest(t_val, samp_rho, go_term)$p
pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]

p.mat[p, ] <- c(pval1, pval2)

}





#########   EVALUATION OF THE TEST FUNCTION --------------------------------

## rewrite the enrichment test function that allows true parameter to be used ---


EnrichTest <- function(t.val, xi, sigma.t, beta0, samp.rho, go.term)
{
  ones <- rep(1, length(t.val))
  
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


# generate null of t_vals and calculate corresponding p values

p.mat <- matrix(NA, nrow=100, ncol=2)

for (p in 1: 100) {
n_gene <- m <- 500                              # number of genes
library(MASS)

## generate correlation matrix
rho <- 0.4
cor.struct <- matrix(rho, n_gene, n_gene);diag(cor.struct) <- 1
sigma <- 1.5*cor.struct
dat <- mvrnorm(1000, mu=rep(0, n_gene), sigma)
sigma1 <- cor(dat)
sigma2 <- 1.5*sigma1
dat <- mvrnorm(1000, mu=rep(0, n_gene), sigma2)

#r1 <- cor(dat)                                  # True correlation matrix
r1 <- diag(n_gene)                              # make the t_vals independent of each other

xi <- 0.8                                       # true xi
sigma.t <- 4                                    # true sigma2.t
beta0 <- -2                                     # true beta0
ones <- rep(1, n_gene)

SIGMA <- sigma.t * ( (1- xi)*diag(1, n_gene)  + xi* r1)  # the covariance
t_val <- mvrnorm(n=1, mu=beta0*ones, SIGMA)              # generate the t values


go_id <- sample(1:n_gene, 0.2*n_gene)         # random sample of go term 
go_term <- rep(0, n_gene); 
go_term[go_id] <- 1
samp_rho <- r1

pval1 <- EnrichTest(t_val, xi, sigma.t, beta0, samp.rho, go_term)$p
pval1 <-  EnrichTest(t_val, samp_rho, go_term)$p
pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]

p.mat[p, ] <- c(pval1, pval2)

}




p_mat_uncorrelated <- read.table("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/Enrichment_testStatisticIndepT_0830.txt")
hist(p_mat_uncorrelated[, 1], main="our test")                           # this test seems to be anti-conservative
hist(p_mat_uncorrelated[, 2], main="linear regression")                            # looks uniform
plot(p_mat_uncorrelated[, 2],p_mat_uncorrelated[, 1], cex=0.5, pch=20, xlab="lm", ylab="enrich_test")

## the two test p values are identical


p_mat_correlated <- read.table("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/Enrichment_testStatistic_correlatedT_0830.txt")
hist(p_mat_correlated[, 1], main="our test")
hist(p_mat_correlated[, 2], main="linear regression")
# our test is more conservative


p_mat_uncorrelated_estimate <- read.table("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/Enrichment_testStatistic0830_estimatedCor_Indep.txt")
hist(p_mat_uncorrelated_estimate[, 1], main="our test")
hist(p_mat_uncorrelated_estimate[, 2], main="linear regression")
plot(p_mat_uncorrelated_estimate[, 2],p_mat_uncorrelated_estimate[, 1], cex=0.5, pch=20, xlab="lm", ylab="enrich_test")
## if t_vals are independent, and we estiamte the sample correlation, the p values look OK
  
# if there is DE, our test is more conservative  
p_mat_de <- read.table("/Users/Bin/Google Drive/Study/Thesis/Correlation/Data/Enrichment_testStatistic_DE_0830.txt")
hist(p_mat_de[, 1], main="our test")
hist(p_mat_de[, 2], main="linear regression")



