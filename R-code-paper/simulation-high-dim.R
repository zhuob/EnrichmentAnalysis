# this method is specific to high dimensional data simulation under correlation 
# structures (a) - (e)
simulate_expression_data_high_dim <- function(size, n_gene, n_test, prop, de_mu, de_sd,
                                              rho1, rho2, rho3, seed = 123, data_gen_method = NULL){
  set.seed(seed)
  # simulate go term, and delta vector
  obj <- meaca:::prepare_simulation(rho1 = rho1, rho2 = rho2, rho3 = rho3, n_gene = n_gene, 
                            n_test = n_test, prop = prop, de_mu = de_mu, 
                            de_sd = de_sd, seed = seed)
  
  go_term <- obj$go_term
  delta <- obj$delta
  
  # the following row is not used, but for consistency of function input purpose
  dat_gen_method <- data_gen_method 
  # which correlation structure it is simulating
  # sim_case <- get_sim_case(rho1 = rho1, rho2 = rho2, rho3 = rho3)  
  # initialize the expression matrix
  expression_data <- matrix(NA, nrow = n_gene, ncol = size)
  n <- size/2
  go_term_genes <- which(go_term == 1)
  
  ## Simulate expression data according to correlation structure --------------- 
  # xi = sqrt(rho1-abs(rho3))*e1 + sqrt(1-rho1)*ei + sqrt(abs(rho3))*e3
  # xj = sqrt(rho2-abs(rho3))*e2 + sqrt(1-rho2)*ej + sqrt(abs(rho3))*e3
  # make sure that cov(x[p], x[q]) = rho1 if (p, q) belong to test genes
  #                cov(x[p], x[q]) = rho2 if (p, q) belong to background genes 
  # and            cov(x[p], x[q]) = rho3 if (p, q) belong to (test, background) genes
  
    for(i in 1:size){
      
      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)
      mu1 <- rnorm(n_gene); mu2 <- rnorm(n_gene)
      y1_temp <- sqrt(rho1-abs(rho3))*e1 + sqrt(1-rho1)*mu1 + sqrt(abs(rho3))*e3
      y2_temp <- sqrt(rho2-abs(rho3))*e2 + sqrt(1-rho2)*mu2 + sign(rho3)*sqrt(abs(rho3))*e3
      y1 <- rep(NA, n_gene)
      y1[ go_term_genes] <- y1_temp[ go_term_genes]
      y1[-go_term_genes] <- y2_temp[-go_term_genes]
      
      if(i <= n){ # simulate control
        expression_data[, i] <- y1
      } else{ # simulate treatment
        expression_data[, i] <- y1 + delta
      }
    }
  
  ## end of simulation data ----------------------------------------------------
  
  trt <- rep(c(0, 1), c(n, size - n))              # group indicator
  
  rownames(expression_data) <- paste("Gene", 1:n_gene, sep="")
  colnames(expression_data) <- paste("trt", trt, sep = "")
  ls <- list(data=expression_data, trt = trt, go_term = go_term, 
             sigma = obj$sigma)
  
  return(ls)
  
}
  
# dat <- simulate_expression_data_high_dim(size = 50, n_gene = 500, n_test = 100, 
#                                          prop = c(0.1, 0.1), de_mu = 2, de_sd = 1,
#                                          rho1 = 0.1, rho2 = 0.05, rho3 = -0.05, seed = 123)

