library(tidyverse)
parent_folder <- "/home/stats/zhuob/EnrichmentAnalysis/"
# parent_folder <- ""
source_r_file <- "R-code-paper/sim-wrapper.R"
dest <- "simulation-results/gene-500/n25vs25/DEmean1sd1"

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


for(i in 1:5){
  
  if(i == 1){run_case <- case_a}
  else if(i == 2) {run_case <- case_b}
  else if(i == 3) {run_case <- case_c}
  else if(i == 4) {run_case <- case_d}
  else if(i == 5) {run_case <- case_e}
  
  # additional simulations: small sample size 
  t1 <- data_simu(nsim = 5,
                  ncore = 2,   
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
                  size = 50,
                  de_mu = 1,
                  de_sd = 1,
                  data_gen_method = "chol",
                  seed = 123)
}
