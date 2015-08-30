## use method of moment to estimate var(delta) and then conduct the
## enrichment analysis...

estimate.sigma <- function(microarray, trt){
  
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  samp_rho <- cor(t(resid_mat))                 
  
  factors <- unique(trt)
  m1 <- sum(trt == factors[1])
  m2 <- sum(trt == factors[2]) 
  sigma2 <- 1/m1 + 1/m2
  
  var_diag <- apply(resid_mat, 1, var)
  sigma0 <- var_diag - sigma2
  sigma0[sigma0 < 0] <- 0
  
  
  # what if I don't use samp_rho but some random correlation matrix
  #n_gene <- nrow(microarray); size <- length(trt)
  #zm <- rnorm(n_gene*size)
  #dim(zm) <- c(n_gene, size)
  #samp_rho <- cor(t(zm))
  
  # What if I use the true correlation?
 # samp_rho <- dat$sigma                                             # this has been passed by compare_test()
 # print(samp_rho[1:10, 1:10])
  ndim <- nrow(microarray)
  sigma <-  sigma0*diag(ndim) + sigma2*samp_rho

  size <- length(trt)
  idex1 <- which(trt == factors[1])[1]
  idex2 <- which(trt == factors[2])[1]
  t_val <- group_mean[, idex2] - group_mean[, idex1]            # control - treatment
  
  
  return(list(sigma = sigma, t_val = t_val))

}
  

#####  our test
MOM_test <- function(microarray, trt, go_term){
  
   est_sigma <- estimate.sigma(microarray, trt)
   sigma <- est_sigma$sigma
   t_val <- est_sigma$t_val
  
   beta0 <- mean(t_val)
   ones <- rep(1, length(t_val))
   
   ## calculate the test statistics
   inv_sigma <- solve(sigma)
   numer <- (crossprod(go_term, (t_val - beta0*ones)))^2
   denom <- t(go_term) %*% sigma %*% (go_term) - (crossprod(go_term, ones))^2/(t(ones) %*% inv_sigma %*% ones)
   
   go_bar <- rep(mean(go_term), length(go_term))
   denom2 <- t(go_term-go_bar) %*% sigma %*% (go_term- go_bar)
   print(c(denom, denom2))
   
   test_stat <- drop(numer/denom2)
   
   pval <- 1 - pchisq(test_stat, 1)
   
   return(list(stat = test_stat, p = pval))
}
  

compare_test <- function(dat){
  
  microarray <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term
  sigma <- dat$sigma
  
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  
  samp_rho <- cor(t(resid_mat))                            # sample correlation matrix
  t_val <- group_mean[, size] - group_mean[, 1]            # control - treatment
  
  
  #### our test
  pval1 <- MOM_test(microarray, trt, go_term)$p
  pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression
  
  library(limma)
  design <- model.matrix(~trt) 
  
  fit <- lmFit(microarray, design)
  fit <- eBayes(fit)								      # Emperical Bayes t test	
  stat <- fit$t[, 2]            							      #	
  alter <- "either"                                              			      # undirectional test	 
  index1 <- which(go_term==1)                                                           # which genes are in GOTERM
  
  # up: positive statistics
  # down: negative stats
  # either: the set is either up or down-regulated as a whole
  # mixed: tests whether the genes in the set tend to be DE without regard for direction.
  #        In this case, the test will be significant if the set contains mostly large statistics, negative or positive
  
  tes1 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= F)         # moderated t GeneSet Rank
  tes2 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= T)         # moderated t GeneSet
  tes3 <- camera(microarray, index1,  design)$PValue                                    # camera proedure
  tes4 <- camera(microarray, index1,  design, use.ranks= T)$PValue       
  
  return(c(pval1, pval2, tes1, tes2, tes3, tes4))
  
}



#########  simulation
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/SimulateLabData.R")

nsim <- 1000
size <- 50
n_gene <- 500
go_prop <- 0.0
no_go_prop <- 0.0
delta_mat <- 2
case <- "a"

pval <- data.frame(matrix(NA, nsim, 6))

for ( k in 1:nsim)
{
  dat <- simulate.microarry(size, n_gene, go_prop, no_go_prop, delta_mat, case = case)
  pval[k, ] <- 
    compare_test(dat)

}

colnames(pval)[1:6] <- c("OurTest", "lm",  "ModeratedT", "Rank", "Camera", "CameraRank")

write.table(pval, "/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationMOMEnrichment/cor_f_50.txt")




source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/SimulateLabData.R")
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationMOMEnrichment/")
create.hist2(read.table("DE_a_6.txt"))

create.hist2(read.table("DE_a_50_sample.txt"))
create.hist2(read.table("DE_a_50_random.txt"))
create.hist2(read.table("DE_a_50_true.txt"))

create.hist2(read.table("DE_a_50.txt"))
create.hist2(read.table("DE_b_50.txt"))
create.hist2(read.table("DE_c_50.txt"))
create.hist2(read.table("DE_d_50.txt"))
create.hist2(read.table("DE_e_50.txt"))
create.hist2(read.table("DE_f_50.txt"))


create.hist2(read.table("NO_DE_a_50.txt"))
create.hist2(read.table("NO_DE_c_50.txt"))
create.hist2(read.table("NO_DE_e_50.txt"))
create.hist2(read.table("NO_DE_f_50.txt"))
 

# power: go 20% DE, no_go, 0% DE
create.hist2(read.table("DE_a_50_power.txt"))
create.hist2(read.table("DE_b_50_power.txt"))
create.hist2(read.table("DE_c_50_power.txt"))
create.hist2(read.table("DE_d_50_power.txt"))
create.hist2(read.table("DE_e_50_power.txt"))
create.hist2(read.table("DE_f_50_power.txt"))





source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/SimulateLabData.R")
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationRealCor20151207/")

##  reproduce 
create.hist2(read.table("Reproduce/DE_a0_50.txt"))
create.hist2(read.table("Reproduce/DE_a_50.txt"))
create.hist2(read.table("Reproduce/DE_b_50.txt"))
create.hist2(read.table("Reproduce/DE_c_50.txt"))
create.hist2(read.table("Reproduce/DE_d_50.txt"))
create.hist2(read.table("Reproduce/DE_e_50.txt"))
create.hist2(read.table("Reproduce/DE_f_50.txt"))
create.hist2(read.table("Reproduce/DE_g_50.txt"))
create.hist2(read.table("DE_g_50.txt"))
create.hist2(read.table("DE_g_501.txt"))

resid <- readRDS("arab_residual.rds")
samp <- sample(1:nrow(resid), 1000)
x1 <- resid[samp, ]
sigma <- cor(t(x1))

zm <- rnorm(n_gene*100)
dim(zm) <- c(n_gene, 100)
sigma2 <- cor(t(zm))
mean(sigma2)

##  standardization 
create.hist2(read.table("DE_a0_50.txt"))
create.hist2(read.table("DE_a_50.txt"))
create.hist2(read.table("DE_b_50.txt"))
create.hist2(read.table("DE_c_50.txt"))
create.hist2(read.table("DE_d_50.txt"))
create.hist2(read.table("DE_e_50.txt"))
create.hist2(read.table("DE_f_50.txt"))
create.hist2(read.table("DE_g_50.txt"))
create.hist2(read.table("DE_h_50.txt"))

create.hist2(read.table("NO_DE_c_50.txt"))


##  standardization
create.hist2(read.table("DE_a0_50_New.txt"))
create.hist2(read.table("DE_a_50_new.txt"))
create.hist2(read.table("DE_b_50_new.txt"))
create.hist2(read.table("DE_c_50_new.txt"))
create.hist2(read.table("DE_d_50_new.txt"))
create.hist2(read.table("DE_e_50_new.txt"))
create.hist2(read.table("DE_f_50_new.txt"))
create.hist2(read.table("DE_g_50_New.txt"))
create.hist2(read.table("DE_h_50_new.txt"))


create.hist2(read.table("NO_DE_a_50_New.txt"))
create.hist2(read.table("NO_DE_c_50_New.txt"))
create.hist2(read.table("NO_DE_e_50_New.txt"))
create.hist2(read.table("NO_DE_f_50_New.txt"))
create.hist2(read.table("NO_DE_g_50_New.txt"))




create.hist2(read.table("Power_a_50_New.txt"))
create.hist2(read.table("Power_c_50_New.txt"))
create.hist2(read.table("Power_e_50_New.txt"))
create.hist2(read.table("Power_f_50_New.txt"))
create.hist2(read.table("Power_g_50_New.txt"))


n <- 1000; nu <- 4;tau <- 0.25
x1 <- 1/rchisq(n, nu)                                   # inverse chi square distribution
s2 <- tau^2*nu*x1                                       # scaled inverse chi-square distribution
stdevs <- sqrt(s2)        

