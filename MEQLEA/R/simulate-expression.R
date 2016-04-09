#' Simulate test set gene labels
#' 
#' @title Simulate gene labels.
#' @param num_gene  a vector of length 2, total number of genes to be simulated and number of genes in the test set.
#' @param prop   a vector of length 2, proportion of DE genes within go term and outside go_term, corresponding to $p_t$ and $p_b$.
#' @param delta  a vector of length num_gene, representing DE effect for each gene
#' @return a list
#' \item{go_term}{ DE indicator for each simulated gene. 1 is for genes in the test set, and 0 othwerwise.}
#' \item{delta}{DE effect for each simulated gene.}
#' \item{n_gene}{ number of genes to be simulated.}
#' @export
#' @examples




prepare_go_term <- function(num_gene, prop, delta){
  # num_gene:  a vector of length 2, number of genes to be simulated and number of go_term genes
  # prop :  a vector of length 2, proportion of DE genes within go term and outside go_term
  # delta: magnitidue of DE, length = num_gene[1]
  
  
  ## generate go term and non-go term 
  n_gene <- num_gene[1]; go_gene <- num_gene[2]
  go_prop <- prop[1]; no_go_prop <- prop[2]
  
  go_term <- sample(1:n_gene, go_gene)                        # sample the go genes
  go_de <-go_term[which(rbinom(go_gene, 1, go_prop)==1)]      # randomly assign DE to genes
  no_go_term <- (1:n_gene)[-go_term]                          # the rest are non go term genes
  no_go_de <- no_go_term[which(rbinom(n_gene-go_gene, 1, no_go_prop)==1)]; # which genes are DE in non-go term
  
  if (length(go_de) == 0 & length(no_go_de) == 0) {
    delta <- rep(0, n_gene)                                  # if no DE in both go_term and no_go_term
  }
  delta[-c(go_de, no_go_de)] <- 0   
  
  z1 <- rep(0, n_gene)
  z1[go_term] <- 1
  go_term <- z1                                           # prepare the GO term
  
  obj <-  list(go_term = go_term, delta = delta, n_gene = n_gene)
  return(obj)
}




#' Prepare simulation parameters. 
#' 
#' 
#' @title Simulate data
#' @param num_gene  number of genes to be simulated.
#' @param prop   a vector of length 2, corresponding to $p_t$ and $p_b$. 
#' @param delta  a vector of length num_gene, representing DE effect for each gene
#' @param case   which case to simulate (see simulation setup in the paper)
#' @param rho   a vector of length 3. Corresponds \eqn{\rho_1}, \eqn{\rho_2} and \eqn{\rho_3}. 
#' @return a list 
#' \item{go_term}{ DE indicator for each simulated gene. 1 is for genes in the test set, and 0 othwerwise.}
#' \item{delta}{DE effect for each simulated gene. }
#' \item{sigma}{the true correlation, used to simulate  the expression data}
#' @export
#' @examples






################ prepare simulation --------------
prepare_simulation <- function(num_gene, prop, delta, case, rho){
  ## write a function to generate correlation matrix, go_term, and size of DE
  ## rho, vector for case e & f to work
  # cor(go_term) = rho[1], cor(no_go_term) = rho[2], cor(go_term, no_go_term) = rho[3]
  
  ######################## NOTE ######################
  ## reshuffled the labels as follows
  # a0 -> a
  # a -> c
  # b -> f
  # c -> b
  # d -> g
  # e -> d 
  # f -> e
  # h ->   deleted 
  rho1 <- rho[1]; rho2 <- rho[2]; rho3 <- rho[3]
  obj <- prepare_go_term(num_gene, prop, delta=delta)
  n_gene <- obj$n_gene; 
  go_term <- obj$go_term
  delta <- obj$delta
  
  #sigma <- diag(n_gene)  # initialize the covariance matrix
  #   print(obj)	
  if (case == "a"){
    sigma <- diag(n_gene)
  }
  
  else if (case == "c"){
    
    sigma <- diag(n_gene)
    ids <- which(go_term ==1)
    sigma[ids, ids] <- rho1
    diag(sigma) <- 1
  }
  
  else if (case == "b") {
    
    sigma <- matrix(rho1, n_gene, n_gene)
    diag(sigma) <- 1
    
  }
  
  else if (case == "d"){
    
    sigma <- diag(n_gene)
    ids <- which(go_term ==1)
    sigma[ids, ids] <- rho1
    sigma[-ids, -ids] <- rho2
    diag(sigma) <- 1
  }
  
  else if (case == "e"){
    
    sigma <- diag(n_gene)
    ids <- which(go_term ==1)
    sigma[ids, ids] <- rho1
    sigma[-ids, -ids] <- rho2
    sigma[ids, -ids] <- rho3   ## negative correlation between go term and non-go term
    sigma[-ids, ids] <- rho3   ##
    diag(sigma) <- 1
  } 
  
  else if (case == "f"){
    
    set.seed(100)
    sigma <- diag(n_gene)
    ids <- which(go_term ==1)
    zm <- rnorm(length(ids)*1000)
    dim(zm) <- c(length(ids), 1000)
    sigma.set <- cor((t(zm)))
    
    sigma[ids, ids] <- sigma.set
    diag(sigma) <- 1
    
    set.seed(NULL)
  }
  
  
  else if (case == "g"){
    
    set.seed(100)
    zm <- rnorm(n_gene*1000)
    dim(zm) <- c(n_gene, 1000)
    sigma <- cor(t(zm))
    set.seed(NULL)
  }
 
  
  return(list(go_term = go_term, delta = delta, sigma=sigma))
}





simu_normal <- function(n, mu, sigma, method = "chol"){

    #### generate correlated multivariate normal data, could use package "mvtnorm"  or "MASS"
  # method: "chol" (recommended for high dimensional, much faster)
  #         "eigen"  or "svd",     or "MAS" [method from MASS]
  #     mu: mean vectors
  #  sigma: covariance structure
  #      n: number of samples
  
  if (method == "MAS"){
    data <- mvrnorm(n, mu, sigma)
  }
  else {
    data <- rmvnorm(n, mean = mu, sigma = sigma, method=method)
  }
  
  return(data)
}




#' simulate normally distributed expression data with desired DE probabilities for genes in the test set and for those not in the test set..
#'
#' @title Simulate expression data. 
#' @param obj  object returned by \code{prepare_simulation}.
#' @param size number of samples to be simulated
#' @return a list 
#' \item{data}{a expression matrix of \eqn{m\times n} where m is the number of genes and n is the number of samples. }
#' \item{trt}{sample labels of length \code{n}, 1 for treatment and 0 for control.}
#' \item{go_term}{gene labels of length \code{m}, 1 for go_term genes and 0 otherwise.}
#' \item{sigma}{true covariance matrix upon which data is simulated.}
#' @export
#' @examples




simulate_expression_data <- function(size, obj){
  # obj: object returned by prepare.simulation()  or  real.go.term()
  #  size: number of samples to be simulated
  # case: which covariance structure
  
  
  go_term <- obj$go_term
  delta <- obj$delta
  cor_sigma <- obj$sigma
  n_gene <- length(delta)
# cor2cov <- CovSimu(4, 0.25, cor_sigma)            # the parameters are set to be the same as Smyth's 2012 paper
#  sigma <-  cor2cov$sigma 
  sigma <- cor_sigma 				  # for now I just use the correlation as Cov
  
  ############ generate expression matrix  ------------------------
  
  mu_initial <- 0
  mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
  mu2 <- delta  + mu_initial                     # mean expression value for treatment
  
  expression_data <- matrix(nrow=n_gene, ncol=size)   # expression matrix
  n <- size/2
  expression_data[, 1:n] <- t(simu_normal(n, mu1, sigma, method="MAS"))
  expression_data[, -(1:n)]  <- t(simu_normal(n, mu2, sigma, method="MAS"))
  trt <- rep(c(0, 1), each= size/2)              # group indicator
  
  rownames(expression_data) <- paste("Gene", 1:n_gene, sep="")
  colnames(expression_data) <- paste("trt", trt, sep = "")
  ls <- list(data=expression_data, trt = trt, go_term = go_term, sigma = sigma)
  return(ls)
  
}

