library(tidyverse)
source("R-code-paper/sim-wrapper.R")

use_r_file <- c("R-code-paper/compare-methods.R")

case_a <- c(0,      0,     0)
case_b <- c(0.1,  0.1,   0.1)
case_c <- c(0.1,    0,     0)
case_d <- c(0.1, 0.05,     0)
case_e <- c(0.1, 0.05, -0.05)

run_case <- case_b

# additional simulations: small sample size 
t1 <- data_simu(nsim = 100,
                ncore = 3,   
                package_used = c("MASS", "qusage"), 
                verbose_show = TRUE, 
                meaca_only = FALSE,
                file_to_source = use_r_file,
                dest = "simulation-results/gene-500/n10vs10/", 
                n_gene = 500,
                n_test = 100,
                prop = c(0.0, 0.0),
                rho1 = run_case[1], 
                rho2 = run_case[2], 
                rho3 = run_case[3],
                size = 50,
                de_mu = 2,
                de_sd = 1,
                data_gen_method = "chol",
                seed = 123)
