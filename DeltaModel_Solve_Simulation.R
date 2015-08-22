

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

group_mean <- group.mean(arab[[1]], trt[[1]])
resid_mat <- arab[[1]] - group_mean

for (i in 2:length(trt) ) {
  group_mean <- group.mean(arab[[i]], trt[[i]])       # calculate the group means within each data set
  resid_temp <- arab[[i]] - group_mean                # the residual after the mean is removed
  resid_mat <- merge(resid_mat, resid_temp, by=0)     # merge the residuals 
  rname <- resid_mat[, 1]                             # keep the row names
  resid_mat <- resid_mat[, -1]                        
  row.names(resid_mat) <- rname                       # give the rows names
}


#####  get a random sample from the Genes, and match them with a GO term. 

n_gene <- 500                                         # number of genes to be considered
ids <- which(rowMeans(arab[[1]]) > 3)                 # keep the genes that have at least 3 counts for each sample
#set.seed(20)                                         # this gene list can be repeated
get_id <- sample(ids, n_gene)
total_gene <- arab[[1]][get_id, ]                     # sample 2000 genes from the list
sample_size <- dim(total_gene)[2]

samp_rho <- cor(t(resid_mat[get_id, ]))

##    the GO term genes 
# go <- read.table("~/Google Drive/Study/Thesis/Correlation/Data/DNAmetabolism.txt", header=F)
go <- read.table("~/Google Drive/Study/Thesis/Correlation/Data/ProteinMetabolism.txt", header=F)
go <- go[, 1]
go.in <- which(noquote(row.names(total_gene)) %in% go)    # which Go term genes are present in this list
go_term <- rep(0, n_gene)
go_term[go.in] <- 1
# go_term <- rbinom(n_gene, 1, 0.2)

x_bar <- apply(total_gene[, 1:sample_size/2], 1, mean)
y_bar <- apply(total_gene[, -(1:sample_size/2)], 1, mean) 

t_val <- (x_bar-y_bar)# /apply(total_gene, 1, sd)           # mean of treat - mean of control

mean(t_val[go_term==1])                                   # mean diff within go term
mean(t_val[go_term==0])                                   # mean diff outside go term
length(go.in)

xlable <- 1:n_gene
plot(t_val[-go.in], xlable[-go.in], pch=20, cex=0.5)
points(t_val[go.in], xlable[go.in], pch=3, col="red")


# go_term <- ifelse(t_val<0, 1, 0)
result <- EnrichTest(t_val, samp_rho, go_term)
result




############## Simulation Study -------------------------------------------
size <- 50                                     # number of samples to be considered
n_gene <- 500                                  # number of genes to be considered
de_prop <- 0.2                                 # proportion of DE
de <- n_gene*de_prop                           # number of DE genes
delta <- rexp(n_gene, 0.2)                     # magnitude of DE
z <- rbinom(n_gene, 1, de_prop)                # indicator, DE =1, NON_DE=0 
delta <- delta*z

rho <- 0.5                                     # correlation matrix
sigma0 <- 26                                    # variance part 
cor.struct <- matrix(rho, n_gene, n_gene);    
diag(cor.struct) <- 1                          # correlation matrix
sigma <- sigma0*cor.struct                     # covariance matrix for both groups

mu_initial <- 0
mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
mu2 <- delta  + mu_initial                     # mean expression value for treatment


microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
n <- size/2
set.seed(50)
microarray[, 1:n] <- t(mvrnorm(n, mu1, sigma)) 
set.seed(51)
microarray[, -(1:n)]  <- t(mvrnorm(n, mu2, sigma))

trt <- rep(c(0, 1), each= size/2)             # group indicator


# prepare the t_stat, samp_rho and Go term 

source("SimulateLabData.R")

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



mean(t_val[go_term==1])                                   # mean diff within go term
mean(t_val[go_term==0])                                   # mean diff outside go term
length(go.in)

xlable <- 1:n_gene
plot(t_val[-go_in], xlable[-go_in], pch=20, cex=0.5)
points(t_val[go_in], xlable[go_in], pch=3, col="red")


EnrichTest(t_val, samp_rho, go_term)




