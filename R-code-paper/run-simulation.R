library(tidyverse)
parent_folder <- "/home/stats/zhuob/EnrichmentAnalysis/"
source_r_file <- "R-code-paper/sim-wrapper.R"
dest <- "simulation-results/gene-500/n10vs10/"

source_r_file <- paste0(parent_folder, source_r_file)


source(source_r_file)

use_r_file <- c("R-code-paper/compare-methods.R",
                "R-code-paper/other-methods.R", 
                "R-code-paper/GSEA.1.0.R")

use_r_file <- paste0(parent_folder, use_r_file)
dest <- paste0(parent_folder, dest)

case_a <- c(0,      0,     0)
case_b <- c(0.1,  0.1,   0.1)
case_c <- c(0.1,    0,     0)
case_d <- c(0.1, 0.05,     0)
case_e <- c(0.1, 0.05, -0.05)

run_case <- case_c

# additional simulations: small sample size 
t1 <- data_simu(nsim = 10000,
                ncore = 40,   
                package_used = c("MASS", "qusage"), 
                verbose_show = FALSE, 
                meaca_only = FALSE,
                file_to_source = use_r_file,
                dest = dest, 
                n_gene = 500,
                n_test = 100,
                prop = c(0.0, 0.0),
                rho1 = run_case[1], 
                rho2 = run_case[2], 
                rho3 = run_case[3],
                size = 20,
                de_mu = 2,
                de_sd = 1,
                data_gen_method = "chol",
                seed = 123)
