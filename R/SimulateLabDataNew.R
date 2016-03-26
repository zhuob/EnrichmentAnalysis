
library(MASS)
library(dplyr)


CovSimu <- function(nu, tau, corMat){
# convert correlation matrix to covariance matrix
# nu: degree of freedom for inversed-chi-square distribution
# tau: scaling parameter for inversed-chi-square distribution
# see https://en.wikipedia.org/wiki/Scaled_inverse_chi-squared_distribution

	n <- dim(corMat)[1]
  	x1 <- 1/rchisq(n, nu)                                   # inverse chi square distribution
  	s2 <- tau^2*nu*x1                                       # scaled inverse chi-square distribution
  	stdevs <- sqrt(s2)                                      # get the standard errors
  	b <- stdevs %*% t(stdevs)
  	covariance <- b * corMat
  
	return(list(sigma = covariance, sd = stdevs))
}



SimuNormal <- function(n, mu, sigma, method = "chol"){
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





go.term <- function(num_gene, prop, delta){
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
  
 	obj <-  list(go_term = go_term, delta = delta, n_gene1 = n_gene)
    	return(obj)
}
  


################ prepare simulation --------------

prepare.simulation <- function(num_gene, prop, delta, case, rho){
## write a function to generate correlation matrix, go_term, and size of DE
## rho, vector for case e & f to work
# cor(go_term) = rho[1], cor(no_go_term) = rho[2], cor(go_term, no_go_term) = rho[3]
  
  
    rho1 <- rho[1]; rho2 <- rho[2]; rho3 <- rho[3]
    ## a major mistake here !!!! the DE effect is set to be N(2, 1)
    
    obj <- go.term(num_gene, prop, delta= delta)
    n_gene <- obj$n_gene1; 
    go_term <- obj$go_term
    delta <- obj$delta

#   print(obj)	
    if (case == "a0"){
      sigma <- diag(n_gene)
    }
    
   else if (case == "a"){
      sigma <- diag(n_gene)
      ids <- which(go_term ==1)
      sigma[ids, ids] <- rho1
      diag(sigma) <- 1
  }
  
  else if (case == "b"){
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
  
  else if (case == "c") {
      sigma <- matrix(rho1, n_gene, n_gene)
      diag(sigma) <- 1
    
  }
  
  else if (case == "d"){
      set.seed(100)
      zm <- rnorm(n_gene*1000)
      dim(zm) <- c(n_gene, 1000)
      sigma <- cor(t(zm))
      set.seed(NULL)
  }
  
  else if (case == "e"){
      
      sigma <- diag(n_gene)
      ids <- which(go_term ==1)
      sigma[ids, ids] <- rho1
      sigma[-ids, -ids] <- rho2
      diag(sigma) <- 1
  }
  
  else if (case == "f"){
    
      sigma <- diag(n_gene)
      ids <- which(go_term ==1)
      sigma[ids, ids] <- rho1
      sigma[-ids, -ids] <- rho2
      sigma[ids, -ids] <- rho3   ## negative correlation between go term and non-go term
      sigma[-ids, ids] <- rho3   ##
      diag(sigma) <- 1
  }
	
	# similar to d	
   else if (case == "h"){			
	arab_res <- readRDS("/home/stats/zhuob/data/computing/arab_residual.rds")
	#set.seed(100)
	rows_samp <- sample(1:nrow(arab_res), n_gene)
	resid <- arab_res[rows_samp, ]
	sigma <- cor(t(resid))	
	}
  
  return(list(go_term = go_term, delta = delta, sigma=sigma))
}





real.go.term <- function(nth_largest, num_gene, prop, delta){
## select the n-th largest go term gene
# this function uses the residual matrix and go_term_database as input, and 
# returns the correlation matrix and the go term OF Desired number of genes
# prop :  a vector of length 2, proportion of DE genes within go term and outside go_term
# delta: magnitidue of DE, length = num_gene[1]
  
    
    	#arab_res <- readRDS("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/arab_residual.rds")
    	arab_res <- readRDS("/home/stats/zhuob/data/computing/arab_residual.rds")
	AllGoterm <- read.table("/home/stats/zhuob/data/computing/ATH_GO_GOSLIM.txt", sep="\t", quote = "", header=F, fill=T) %>%	
	# AllGoterm <- read.table("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/ATH_GO_GOSLIM.txt", sep="\t", quote = "", header=F, fill=T) %>%
	filter(grepl('GO:', V6))              # nrow(AllGoterm)  # 251297
  
    	## remove duplicated obs of combination:  GOterm *geneID
    	AllGoterm <- AllGoterm[!duplicated(AllGoterm[c("V1", "V6")]),]
    	nrow(AllGoterm)  # 146108
  
    	GotermList <- group_by(AllGoterm, V6) %>%              # find the go terms
   	summarise(NumberofGene=n()) %>%        # count how many genes in each go term
    	filter(NumberofGene >= 100)            # select Go terms that have at least 100 genes
   
    	idx <- sort(GotermList$NumberofGene,  decreasing = T)[nth_largest]
    	extract_go_term <- as.matrix(GotermList[GotermList$NumberofGene==idx, 1]) # the 4th largest go term
    	gene <- AllGoterm[ AllGoterm$V6 == drop(noquote(extract_go_term)), 1]    # the size of this go term
    
   	gene <- unlist(lapply(gene, as.character))				 # this is the list of go term genes, 9749 in total
    	# which of the go term genes are in the data set
    	#  id_gene <- gene %in% rownames(arab_res)
    	id_gene <- rownames(arab_res) %in% gene					# total number of go_term genes in arab
    	go_term <- rep(0, dim(arab_res)[1])					# create go_term indicator
    	go_term[id_gene] <- 1							# sum(go_term==1) = 5061 
#    	res <- cor(t(arab_res))
 
  
### Updated Dec 16th sample GO_term and non_go_term with fixed numbers    ###
    go_prop <- prop[1]; no_go_prop <- prop[2];
    # subset the arab to two parts: the go_term and non_go_term
    subset1 <- which(go_term ==1)
    length(subset1)  # 5061,  		DEPENDING ON THE RESIDUAL MATRIX, THIS MAY VARY
    subset2 <- which(go_term ==0)
    length(subset2)  # 15906            
    
    ## sample genes from go term and non_go term separately
    s1 <- num_gene[2]      			# genes in the go term
    s2 <- num_gene[1]-num_gene[2]		# number of genes not in go term 
    
    g1 <- sample(subset1, s1)			# sampled go term genes row number
    g2 <- sample(subset2, s2)			# sampled non go term genes row number
    
    gene_select <- c(g1, g2)
    
### end of update ############    
    
   
#    	gene_select <- sample(1:dim(arab_res)[1], n_gene)
    	go_term <- go_term[gene_select]
    
    	sigma <- cor(t(arab_res[gene_select, ]))	# sigma used to simulate correlation structure
    	go_prop <- prop[1]; no_go_prop <- prop[2]
    	go_gene <- sum(go_term == 1)
    
    	go <- which(go_term == 1)           # the go term 
    	no_go <- which(go_term == 0)        # the non-go term
    	go_de <- go[which(rbinom(go_gene, 1, go_prop)==1)]      # randomly assign DE to genes
    	no_go_de <- no_go[which(rbinom(n_gene-go_gene, 1, no_go_prop)==1)]; # which genes are DE in non-go term
    	# length(go_de)/length(go)
    	# length(no_go_de)/length(no_go)
    	if (length(go_de) == 0 & length(no_go_de) == 0) {
      	delta <- rep(0, n_gene)                                  # if no DE in both go_term and no_go_term
   	 }
    	delta[-c(go_de, no_go_de)] <- 0   
    
    	return(list(go_term = go_term, delta = delta, sigma = sigma))
    
}




simulate.microarry <- function(size, obj){
  # obj: object returned by prepare.simulation()  or  real.go.term()
  #  size: number of samples to be simulated
  # case: which covariance structure
  
  
	go_term <- obj$go_term
  	delta <- obj$delta
  	cor_sigma <- obj$sigma
  	n_gene <- length(delta)
  	cor2cov <- CovSimu(4, 0.25, cor_sigma)            # the parameters are set to be the same as Smyth's 2012 paper
  	sigma <-  cor2cov$sigma 
	sigma <- cor_sigma 				  # for now I just use the correlation as Cov

 ############ generate expression matrix  ------------------------
  
  	mu_initial <- 0
 	mu1 <- rep(mu_initial, n_gene)                 # mean expression value for control
  	mu2 <- delta  + mu_initial                     # mean expression value for treatment
  
  	microarray <- matrix(nrow=n_gene, ncol=size)   # expression matrix
  	n <- size/2
  	microarray[, 1:n] <- t(SimuNormal(n, mu1, sigma, method="MAS"))
  	microarray[, -(1:n)]  <- t(SimuNormal(n, mu2, sigma, method="MAS"))
  	trt <- rep(c(0, 1), each= size/2)              # group indicator
  
	rownames(microarray) <- paste("Gene", 1:n_gene, sep="")
	colnames(microarray) <- paste("trt", trt, sep = "")
  	ls <- list(data=microarray, trt = trt, go_term = go_term, sigma = sigma)
  	return(ls)
  
}




get.residual.matrix <- function(set){

	library(NBPSeq)  
	data <- set$count
  	count <- as.matrix(data[rowSums(data)>=3*dim(data)[2], ])   # filter lowly expressed genes
  
 	cpm <- prepare.nb.data(count)$rel.frequencies*1e6           # count perMillion
  	trt <- as.numeric(noquote(paste(set$trt, set$lab, sep="")))
  
  	residua  <- matrix(NA, dim(cpm)[1], dim(cpm)[2])
    
    	sep.mean <- function(y, group){                           # this function is used to calculate the residuals
    		me<- aggregate(y, by = list(trt), FUN = mean)             # and standardize them 
    		names(me) <- c("group", "average")
    		da1 <- data.frame(y, group = trt)
    		newd <- merge(da1, me)
    		resd <- newd$y-newd$average
    		s_resd <- scale(resd, scale=T, center = F)
  	return(s_resd)
	}
  
 	 for (i in 1:nrow(cpm)){
    	residua[i, ] <- sep.mean(cpm[i, ], trt)
 	 }
  
	rownames(residua) <- rownames(cpm)
  	colnames(residua) <- colnames(cpm)

	  return(residua)
}





