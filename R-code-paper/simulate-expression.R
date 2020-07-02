#' Simulate test set gene labels
#' 
#' @keywords internal
#' @title Simulate gene labels.
#' @param n_gene total number of genes to be simulated
#' @param n_test number of genes in the test set.
#' @param prop   a vector of length 2, proportion of DE genes within go term and
#'   outside go_term, corresponding to $p_t$ and $p_b$.
#' @param de_mu,de_sd if the gene is DE, delta ~ N(de_mu, de_sd)
#' @return a tibble with two columns
#' \item{go_term}{ DE indicator for each simulated gene. 1 is for genes in the
#' test set, and 0 othwerwise.}
#' \item{delta}{DE effect for each simulated gene.}
# #' @export
# #' @examples

prepare_go_term <- function(n_gene, n_test, prop, de_mu, de_sd, seed = 123){
  
  ## generate go term and non-go term 
  set.seed(seed)
  # n_gene <- num_gene[1]; n_test <- num_gene[2]
  go_prop <- prop[1]; no_go_prop <- prop[2]
  delta <- rnorm(n_gene, mean = de_mu, sd = de_sd)
    
  # sample the go genes
  go_term <- sample(1:n_gene, n_test)                        
  # randomly assign DE to genes
  go_de <-go_term[which(rbinom(n_test, 1, go_prop)==1)]      
  # the rest are non go term genes
  no_go_term <- (1:n_gene)[-go_term]                          
  # which genes are DE in non-go term
  no_go_de <- no_go_term[which(rbinom(n_gene - n_test, 1, no_go_prop)==1)]; 
  
  if (length(go_de) == 0 & length(no_go_de) == 0) {
    # if no DE in both go_term and no_go_term
    delta <- rep(0, n_gene)                                  
  }
  delta[-c(go_de, no_go_de)] <- 0   
  
  z1 <- rep(0, n_gene)
  z1[go_term] <- 1
  # prepare the GO term
  go_term <- z1                                           
  
  obj <-  tibble::tibble(go_term = go_term, delta = delta)
  return(obj)
}


#' Prepare simulation parameters. 
#' @keywords internal
#' 
#' @title Simulate expression data
#' @param rho1  a scalar, correlation between two test genes (i.e., \eqn{\rho_1} in the paper)
#' @param rho2  a scalar, correlation between two background genes (i.e., \eqn{\rho_2} in the paper)
#' @param rho3  correlation between a test gene and a background gene (i.e., \eqn{\rho_3} in the paper)
#' @param ... other parameters inherited from \code{prepare_go_term}
#' @return a list 
#' \item{go_term}{ DE indicator for each simulated gene. 1 is for genes in the
#' test set, and 0 othwerwise}
#' \item{delta}{DE effect for each simulated gene}
#' \item{sigma}{the true correlation, used to simulate  the expression data}
# #' @export
# #' @examples

################ prepare simulation --------------
prepare_simulation <- function(rho1, rho2, rho3, ...){
  
  obj <- prepare_go_term(...)
  n_gene <- nrow(obj); 
  go_term <- obj$go_term
  delta <- obj$delta
  
  sigma <- diag(n_gene)
  ids <- which(go_term ==1)
  sigma[ids, ids] <- rho1
  sigma[-ids, -ids] <- rho2
  sigma[ids, -ids] <- rho3   
  sigma[-ids, ids] <- rho3
  diag(sigma) <- 1
  
  return(list(go_term = go_term, delta = delta, sigma = sigma))
}


simu_normal <- function(n, mu, sigma, method = "chol"){

  #### generate correlated multivariate normal data, could use package "mvtnorm"
  #### or "MASS"
  # method: "chol" (recommended for high dimensional, much faster)
  #         "eigen"  or "svd",     or "MAS" [method from MASS]
  #     mu: mean vectors
  #  sigma: covariance structure
  #      n: number of samples
  
  if (method == "MASS"){
    data <- MASS::mvrnorm(n, mu, sigma)
  }
  else {
    data <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma, method=method)
  }
  
  return(data)
}




#' simulate normally distributed expression data with desired DE probabilities
#' for genes in the test set and for those not in the test set..
#'
#' @title Simulate expression data. 
#'
#' @param n_gene total number of genes to be simulated
#' @param n_test number of genes in the test set.
#' @param prop   a vector of length 2, proportion of DE genes within go term and
#'   outside go_term, corresponding to $p_t$ and $p_b$.
#' @param rho1  a scalar, correlation between two test genes (i.e., \eqn{\rho_1} in the paper)
#' @param rho2  a scalar, correlation between two background genes (i.e., \eqn{\rho_2} in the paper)
#' @param rho3  correlation between a test gene and a background gene (i.e., \eqn{\rho_3} in the paper)
#' @param size number of samples to be simulated
#' @param de_mu,de_sd if the gene is DE, delta ~ N(de_mu, de_sd)
#' @param data_gen_method data generation method; if `data_gen_method = MASS`,
#'   then \link[MASS]{mvrnorm} is used, otherwise see function \link[mvtnorm]{rmvnorm}
#' @param seed the seed used for simulation (for reproducibility purpose)
#'
#' @return a list 
#' \item{data}{a expression matrix of \eqn{m\times n} where m is the number of
#' genes and n is the number of samples. }
#' \item{trt}{sample labels of length \code{n}, 1 for treatment and 0 for control.}
#' \item{go_term}{gene labels of length \code{m}, 1 for go_term genes and 0 otherwise.}
#' \item{sigma}{true covariance matrix upon which data is simulated.}
#' @export
#' @examples
#' t1 <- simulate_expression_data(size = 50, n_gene = 500, n_test = 100, 
#'                                prop = c(0.1, 0.1), de_mu = 2, de_sd = 1, 
#'                                rho1 = 0.1, rho2 = 0.05, rho3 = -0.05, 
#'                                data_gen_method = "chol", seed = 123)


simulate_expression_data <- function(size, n_gene, n_test, prop, de_mu, de_sd,
                                     rho1, rho2, rho3, data_gen_method = "chol", 
                                     seed = 123){
  # obj: object returned by prepare.simulation()  or  real.go.term()
  #  size: number of samples to be simulated
  # case: which covariance structure
  
  obj <- prepare_simulation(rho1 = rho1, rho2 = rho2, rho3 = rho3, n_gene = n_gene, 
                            n_test = n_test, prop = prop, de_mu = de_mu, 
                            de_sd = de_sd, seed = seed)
  go_term <- obj$go_term
  delta <- obj$delta
  cor_sigma <- obj$sigma
  # the parameters are set to be the same as Smyth's 2012 paper
  # cor2cov <- CovSimu(4, 0.25, cor_sigma)            
  #  sigma <-  cor2cov$sigma 
  sigma <- cor_sigma 				  # for now I just use the correlation as Cov
  
  ############ generate expression matrix  ------------------------
  
  mu_initial <- 0
  # mean expression value for control
  mu1 <- rep(mu_initial, n_gene)                 
  # mean expression value for treatment
  mu2 <- delta  + mu_initial                     
  
  expression_data <- matrix(nrow=n_gene, ncol=size)   # expression matrix
  n <- floor(size/2)
  
  expression_data[, 1:n] <- t(simu_normal(n, mu1, sigma, method = data_gen_method))
  expression_data[, (n+1):size]  <- t(simu_normal(n = size - n, mu2, sigma, 
                                                  method = data_gen_method))
  trt <- rep(c(0, 1), c(n, size - n))              # group indicator
  
  rownames(expression_data) <- paste("Gene", 1:n_gene, sep="")
  colnames(expression_data) <- paste("trt", trt, sep = "")
  ls <- list(data=expression_data, trt = trt, go_term = go_term, sigma = sigma)
  return(ls)
  
}
