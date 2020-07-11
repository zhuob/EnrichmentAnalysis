
size <- 50; n_gene <- 500; n_test <- 100; prop <- c(0.1, 0.1)
de_mu <- 2; de_sd <- 1; seed <- 123; rho1 <- 0.1; rho2 <- 0.05; rho3 <- -0.05
data_gen_method = "chol" 
file_to_source = "~/Desktop/Paper/EnrichmentAnalysis/meaca/R/GSEA.1.0.R"
package_used = c("MASS", "qusage")
verbose_show <- FALSE
meaca_only <- FALSE
ncore <- 2
nsim <- 2

fti <- run_sim_meaca(seed = seed, nsim = nsim, ncore = ncore,  
                     package_used = package_used, verbose_show = verbose_show, 
                     file_to_source = file_to_source, meaca_only = meaca_only, 
                     n_gene = n_gene, n_test = n_test, prop = prop, 
                     rho1 = rho1, rho2 = rho2, rho3 = rho3, 
                     size = size, de_mu = de_mu, de_sd = de_sd, 
                     data_gen_method = data_gen_method)


n <- 100000
mu1 <- c(0, 0, 0)
sigma <- matrix(c(1, 0.1, 0.05, 0.1, 1, -0.05, 0.05, -0.05, 1), 3, 3)
s1 <- simu_normal(n, mu1, sigma, method = "chol")
cor(s1); colMeans(s1)

p1 <- Sys.time()
expression_data[, 1:n] <- t(simu_normal(n, mu1, sigma, method="MAS"))
expression_data[, (n+1):size]  <- t(simu_normal(n = size - n, mu2, sigma, method="MAS"))
trt <- rep(c(0, 1), c(n, size - n))              # group indicator
Sys.time() - p1


p1 <- Sys.time()
expression_data[, 1:n] <- t(simu_normal(n, mu1, sigma, method="chol"))
expression_data[, (n+1):size]  <- t(simu_normal(n = size - n, mu2, sigma, method="chol"))
trt <- rep(c(0, 1), c(n, size - n))              # group indicator
Sys.time() - p1



########################
###
case <- "a"
files <- "TypeIerror_e_0PCT.txt"
destination <-  paste("/home/stats/zhuob/data/computing/", files, sep = "")

nsim <- 100
size <- 50           # number of samples to be simulated

# rho <- c(0.7, 0.5, -0.5)   # extreme case
rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation
num_gene <- c(500, 100)
prop <- c(0.1, 0.1)

run_sim_meaca(nsim = nsim, size = size, rho = rho, num_gene = num_gene, 
              prop = prop, n_gene = num_gene[1], case = case, ncore = 2, 
              )
