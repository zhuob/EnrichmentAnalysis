run_sim_meaca <- function(nsim, ncore = 6, package_used = c("MASS", "qusage"), 
                          verbose_show = FALSE, meaca_only = FALSE,
                          file_to_source = "/home/stats/zhuob/Rcode/Enrichment/GSEA.1.0.R",
                          seed = 123, ...){
  
  # generate seed for reproducing results
  set.seed(seed)
  subseed <- sample(1:1e7, nsim)
  
  platform <- .Platform$OS.type
  if(platform == "windows"){
    cl <- parallel::makeCluster(ncore)                           
  } else{
    cl <- parallel::makeCluster(ncore, type = "FORK")                           
  }
  doParallel::registerDoParallel(cl)                           
  
  library(foreach)
  fti <- foreach(i = 1:nsim, .combine = rbind, .packages = package_used, 
                 .verbose = verbose_show) %dopar% {
    
    if(!is.null(file_to_source)){
      for(kk in 1:length(file_to_source)){
        source(file_to_source[kk])
      }
    }

    # data simulation
    dat <- simulate_expression_data(seed = subseed[i], ...)
    
    # our test, NO standardization for the simulation
    if(meaca_only){
      MEQ <- meaca_single(expression_data = dat$data, trt = dat$trt, 
                          go_term = dat$go_term, standardize = F)             
      pvals <- MEQ$p1							# chi-square test
      names(pvals) <- "meaca"
    } else{
      pvals <- compare_test(dat)
    }
    return(pvals)
  }
  parallel::stopCluster(cl)
  
  return(fti)
}


get_sim_case <- function(rho1, rho2, rho3){
  
  if(rho1 == rho2 & rho1 == rho3 & rho1 == 0){
    case <- "a"
  } else if (rho1 == rho2 & rho1 == rho3 & rho1 == 0.1){
    case <- "b"
  } else if (rho1 == 0.1 & rho2 == rho3 & rho2 == 0){
    case <- "c"
  } else if (rho1 == 0.1 & rho2 == 0.05 & rho3 == 0){
    case <- "d"
  } else if (rho1 == 0.1 & rho2 == 0.05 & rho3 == -0.05){
    case <- "e"
  }
  
  return(case)
}


#' Produce p value matrix for simulation discussed in the paper. 
#' 
#' 
#' @title Compare meaca to existing methods
#'
#' @param n_gene total number of genes to be simulated
#' @param n_test number of genes in the test set.
#' @param prop   a vector of length 2, proportion of DE genes within go term and
#'   outside go_term, corresponding to $p_t$ and $p_b$.
#' @param rho1  a scalar, correlation between two test genes (i.e., \eqn{\rho_1} in the paper)
#' @param rho2  a scalar, correlation between two background genes (i.e., \eqn{\rho_2} in the paper)
#' @param rho3  correlation between a test gene and a background gene (i.e., \eqn{\rho_3} in the paper)
#' @param size number of samples to be simulated
#' @param seed the seed used for simulation (for reproducibility purpose)
#' @param package_used the packages to be used in the simulation
#' @param verbose_show for debug purpose, set to `FALSE` if not in debug mode
#' @param meaca_only Should all the methods to be compared? If `TRUE`, produce
#'   Figure 1; otherwise Figure 2
#' @param file_to_source the R files containing functions to be sourced
#' @param nsim number of simulation to run
#' @param de_mu,de_sd if the gene is DE, delta ~ N(de_mu, de_sd)
#' @param data_gen_method data generation method; if `data_gen_method = MASS`,
#'   then `MASS::mvrnorm` is used, otherwise see function \link[mvtnorm]{rmvnorm}
#' @param seed the seed used for simulation (for reproducibility purpose)
#' @param dest where to store the results
#' @param ncore number of CPUs to be used in the parallel simulation 
#'
#' @return a text file containing the p value matrix
#' @export

data_simu <- function(nsim = 1000,
                      ncore = 6,   
                      package_used = c("MASS", "qusage"), 
                      verbose_show = FALSE, 
                      meaca_only = FALSE,
                      file_to_source = "/home/stats/zhuob/Rcode/Enrichment/GSEA.1.0.R",
                      dest = "/home/stats/zhuob/data/computing/", 
                      n_gene = 500,
                      n_test = 100,
                      prop = c(0.1, 0.1),
                      rho1 = 0.1, 
                      rho2 = 0.05, 
                      rho3 = -0.05,
                      case = "e",
                      size = 50,
                      de_mu = 2,
                      de_sd = 1,
                      data_gen_method = "chol",
                      seed = 123){
  
    sim_case <- get_sim_case(rho1 = rho1, rho2 = rho2, rho3 = rho3)  
  
    if (prop[1] == prop[2]){                            
      # this corresponds to type I error simulation
      files <- paste("TypeIerror_", sim_case, "_", as.integer(prop[1]*100), "PCT.txt", sep= "")
    } else {
      files <- paste("Power_", sim_case, "_", as.integer(prop[1]*100), "VS", 
                     as.integer(prop[2]*100), "PCT.txt", sep = "")
    }
    
    destination <-  paste(dest, "/", files, sep = "")
    
    fti <- run_sim_meaca(seed = seed, nsim = nsim, ncore = ncore,  
                         package_used = package_used, verbose_show = verbose_show, 
                         file_to_source = file_to_source, meaca_only = meaca_only, 
                         n_gene = n_gene, n_test = n_test, prop = prop, 
                         rho1 = rho1, rho2 = rho2, rho3 = rho3, 
                         size = size, de_mu = de_mu, de_sd = de_sd, 
                         data_gen_method = data_gen_method)
    
    write.table(fti, destination, row.names = F)
  
}
