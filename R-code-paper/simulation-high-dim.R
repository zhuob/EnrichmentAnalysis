
#'         Produce p value matrix for simulation discussed in the paper (NO parellel computing).
#'
#'
#' @title Compare meaca to existing methods (NO parellel computing)
#' @param prop A vector of two, specifying the proportion of DE genes for the test set and the background set
#' @param de A vector of two, specifying the (normally distributed) DE effect size and its std
#' @param nsim number of simulations to run
#' @param size number of samples to be simulated
#' @param rho  A vector of three, for correlation coefficients corresponding to rho1, rho2, rho3 in the paper
#' @param dest where to store the results
#' @param num_gene A vector of two, specifying the number of genes to be simulated in total and in the test set
#' @param post_txt formated name of output files
#' @param save  should the results be saved in a file?
#' @return a text file containing the p value matrix
#' @export



data_simu3 <- function(case, 
                       block = 1, 
                       prop = c(0.1, 0.1),             # the proportion of DE genes
                       de = c(2, 1),                   # DE effect size and its std
                       nsim = 1000,                   	# number of simulations to run. 
                       size = 20,                     	# number of samples to be simulated
                       rho = c(0.1, 0.05, -0.05),       # correlations, corresponding to rho1, rho2, rho3 in the paper
                       dest = "/home/stats/zhuob/data/computing/meaca_GB/", # where to store the results
                       num_gene = c(10000, 100),         # number of genes to be simulated
                       post_txt = ".txt",           # format of output files
                       save = TRUE                    # whether save the result
){
  blockid <- sprintf("%05d", block)
  
  if (prop[1] == prop[2]){                            # this corresponds to type I error simulation
    files <- paste("TypeIerror_", case, "_", as.integer(prop[1]*100), "PCT_",  blockid,  post_txt, sep= "")
  } else {
    files <- paste("Power_", case, "_", as.integer(prop[1]*100), "VS", as.integer(prop[2]*100), "PCT_", blockid, post_txt, sep = "")
  }
  
  destination <- paste(dest, "/", files, sep = "")  
  
  delta <- rnorm(num_gene[1], de[1], de[2])

  
  fti <- foreach(i = 1:nsim, .combine = rbind, .packages = c("meaca")) %dopar% 
  {
    ## keep the true de effect delta fixed for each simulation.
    obj <- prepare_simulation3(num_gene, prop = prop, delta)
    dat <- simulate_expression_data3(obj, case = case, size = size, rho = rho)
    pvals <- compare_test(dat)
    return(pvals)
    # rm(list=ls())
  }
  if (save) {
    write.table(fti, destination, row.names = F)
  }
  else {
    return(fti)
  }
  
}


#' Simulate test set gene labels
#' 
#' @title Simulate gene labels.
#' @param num_gene  a vector of length 2, total number of genes to be simulated and number of genes in the test set.
#' @param prop   a vector of length 2, proportion of DE genes within go term and outside go_term, corresponding to $p_t$ and $p_b$.
#' @param delta  a vector of length num_gene, representing DE effect for each gene
#' @return a list
#' \item{go_term}{ DE indicator for each simulated gene. 1 is for genes in the test set, and 0 othwerwise. For simplicity, the GO term genes are set to be gene 1 to gene m1, where m1 is the size of the GO term}
#' \item{go_gene}{DE labels: 1 for DE , 0 for non DE}
#' @export


prepare_simulation3 <- function(num_gene, prop, delta){
  # num_gene:  a vector of length 2, number of genes to be simulated and number of go_term genes
  # prop :  a vector of length 2, proportion of DE genes within go term and outside go_term
  # delta: magnitidue of DE, length = num_gene[1]


  ## generate go term and non-go term
  n_gene <- num_gene[1]; go_gene <- num_gene[2]
  go_prop <- prop[1]; no_go_prop <- prop[2]
  de_gene <- rep(0, n_gene)       # initiate the DE labels: 1 for DE , 0 for non DE

  # sample the go_term
  go_term <- sample(1:n_gene, go_gene)
  go_de_id <- rbinom(go_gene, 1, go_prop)
  go_de <-go_term[go_de_id ==1]      # randomly assign DE to genes
  no_go_term <- setdiff((1:n_gene), go_term)                          # the rest are non go term genes
  no_go_de_id <- rbinom(n_gene - go_gene, 1, go_prop)
  no_go_de <- no_go_term[no_go_de_id==1]; # which genes are DE in non-go term

  k1 <- length(go_de); k2 <- length(no_go_de)

  if (k1*k2 != 0) {  # if DE genes exist in both go_term and no_go_term
    de_gene[c(go_de, no_go_de)] <- 1
  }

  else if (k1 == 0 & k2 != 0){   # no de genes only in the GO term
    de_gene[no_go_de] <- 1
  }
  else if (k1 != 0 & k2 == 0){   # no DE genes only in the background set
    de_gene[go_de] <- 1
  }
  #  the last case, k1 = 0 and k2 = 0, means NO DE genes in the entire genome.
  #  de_gene remains to be all 0

  z1 <- rep(0, n_gene)
  z1[go_term] <- 1
  
  delta[de_gene == 0] <- 0

  obj <-  list(go_term = z1, delta = delta)
  return(obj)
}





# the independent case 
#' simulate expression data 
#' 
#' @title Simulate normally distributed data
#' @param obj   the object returned by \code{\link{prepare_simulation2}}
#' @param size how many sample should be simulated 
#' @param case which corelation structure should be used to simulate data?
#' @param rho  a vector of 3, rho1, rho2 and rho3, respectively.
#' @return a list 
#' \item{data}{a expression matrix of \eqn{m\times n} where m is the number of genes and n is the number of samples. }
#' \item{trt}{sample labels of length \code{n}, 1 for treatment and 0 for control.}
#' \item{go_term}{gene labels of length \code{m}, 1 for go_term genes and 0 otherwise.}
#' @export 


simulate_expression_data3 <- function(obj, 
                                      case = "a",
                                      size,
                                      rho = c(0.1, 0.05, -0.05)){


  go_term <- obj$go_term   # GO term status: 1 for GO term genes and 0 otherwise
  delta <- obj$delta
  n_gene <- length(delta)     # number of genes in the genome
  rho1 <- rho[1]; rho2 <- rho[2]; rho3 <- rho[3]

  go_gene <- which(go_term == 1)    # the GO term genes
  no_go_gene <- which(go_term != 1) # the NON GO term genes
  test_size <- length(go_gene)      # the size of the test set

  expression_data <- matrix(nrow=n_gene, ncol=size)   # expression matrix


  for (i in 1:(size/2)){


    y1 <- rep(NA, n_gene)  # initialize sample 1 and sample 2
    y2 <- rep(NA, n_gene)

    if (case == "a"){

      y1 <- rnorm(n_gene);
      y2_temp <- rnorm(n_gene)
    }

    else if (case == "b"){
      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # the error terms used to generate correlated data
      y1 <- sqrt(rho1)*rep(e1, n_gene) + sqrt(1-rho1)*rnorm(n_gene)
      
      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # run this twice to make y1 and y2 independent
      y2_temp <- sqrt(rho1)*rep(e1, n_gene) + sqrt(1-rho1)*rnorm(n_gene)
    }

    else if (case == "c"){
      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # the error terms used to generate correlated data
      y1[go_gene] <- sqrt(rho1)*rep(e1, test_size) + sqrt(1-rho1)*rnorm(test_size)
      y1[no_go_gene] <- rnorm(n_gene - test_size)

      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # run this twice to make y1 and y2 independent
      y2_temp <- rep(NA, n_gene)
      y2_temp[go_gene] <- sqrt(rho1)*rep(e1, test_size) + sqrt(1-rho1)*rnorm(test_size)
      y2_temp[no_go_gene] <- rnorm(n_gene - test_size)

    }

    else if (case == "d" ){
      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # the error terms used to generate correlated data
      # generate sample 1 for control
      y1[go_gene] <- sqrt(rho1)*rep(e1, test_size) + sqrt(1-rho1)*rnorm(test_size)
      y1[no_go_gene] <- sqrt(rho2)*rep(e2, n_gene - test_size) + sqrt(1-rho2)*rnorm(n_gene-test_size)

      # generate sample 1 for case
      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # run this twice to make y1 and y2 independent
      y2_temp <- rep(NA, n_gene)
      y2_temp[go_gene] <- sqrt(rho1)*rep(e1, test_size) + sqrt(1-rho1)*rnorm(test_size)
      y2_temp[no_go_gene] <- sqrt(rho2)*rep(e2, n_gene - test_size) + sqrt(1-rho2)*rnorm(n_gene-test_size)

    }

    else if (case == "e"){
      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # the error terms used to generate correlated data
      y1[go_gene]<- sqrt(rho1-abs(rho3))*rep(e1, test_size) + sqrt(1-rho1)*rnorm(test_size) + sqrt(abs(rho3))*rep(e3, test_size)
      y1[no_go_gene] <- sqrt(rho2-abs(rho3))*rep(e2, n_gene-test_size) + sqrt(1-rho2)* rnorm(n_gene-test_size) - sqrt(abs(rho3))*rep(e3, n_gene-test_size)

      e1 <- rnorm(1); e2 <- rnorm(1); e3 <- rnorm(1)  # run this twice to make y1 and y2 independent
      y2_temp <- rep(NA, n_gene)
      y2_temp[go_gene] <- sqrt(rho1-abs(rho3))*rep(e1, test_size) + sqrt(1-rho1)*rnorm(test_size) + sqrt(abs(rho3))*rep(e3, test_size)
      y2_temp[no_go_gene] <- sqrt(rho2-abs(rho3))*rep(e2, n_gene-test_size) + sqrt(1-rho2)* rnorm(n_gene-test_size) - sqrt(abs(rho3))*rep(e3, n_gene-test_size)
    }

    ### modificantion made on May 25, 2017
    # de_effect <- rnorm(n_gene, delta, rep(de[2], n_gene))  # fix the DE effect for each gene
    # de_effect[delta == 0] <- 0      # set the DE effect to be 0 for NON-DE genes
    y2 <- y2_temp + delta

    expression_data[, i] <- y1
    expression_data[, i + size/2] <- y2

  }

  trt <- rep(c(0, 1), each = size/2)

  rownames(expression_data) <- paste("Gene", 1:n_gene, sep="")
  colnames(expression_data) <- paste("trt", trt, sep = "")

  ls <- list(data=expression_data, trt = trt, go_term = go_term)

  return(ls)

}

