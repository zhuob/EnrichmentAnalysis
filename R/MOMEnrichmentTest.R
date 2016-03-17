## use method of moment to estimate var(delta) and then conduct the
## enrichment analysis...

library(limma)         # needed for camera() and geneSetTest()


readGeneSet <- function(msigdb){
 ## read the gene sets from MsigDB   
  gene_set <- readLines(msigdb)                                        # the gene sets
  temp <- gene_set
  gs.line <- list()
  
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric") 
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line[[i]] <- noquote(unlist(strsplit(temp[[i]], "\t")))
  }
  return(list(NumofSets = length(temp.size.G), size = temp.size.G, geneSet = gs.line))
}




group.mean <- function(data, group){
  # a function to calculate group mean
  
  id1 <- which(group == 1)
  nrep1 <- length(id1)
  nrep2 <- ncol(data)-nrep1
  mean1 <- apply(data[, id1], 1, mean)
  mean1 <- matrix(rep(mean1, nrep1), ncol=nrep1, byrow=F) 
  mean2 <- apply(data[, -id1], 1, mean)
  mean2 <- matrix(rep(mean2, nrep2), ncol=nrep2, byrow=F) 
  
  data.mean <- data.frame(matrix(NA, nrow(data), ncol(data)))
  data.mean[, id1] <- mean1
  data.mean[, -id1] <- mean2
  colnames(data.mean) <- colnames(data)
  rownames(data.mean) <- rownames(data)
  data.mean
}




standardize.microarray <- function(microarray, trt){
  ## standardize the expression data, to make it have variance 1 in each of treatmeent/control group
  #  first remove the means within each group to get sd, and then original data by sd
  
  group1 <- trt == 1
  group2 <- trt == 0
  col_index <- 1:ncol(microarray)
  
  set1 <- microarray[, group1]
  set2 <- microarray[, group2]
  
  s1 <- apply(set1, 1, sd)
  s2 <- apply(set2, 1, sd)

### modified on Feb 18, 2016, use pooled variacne to do standardization....

  n1 <- length(group1)
  n2 <- length(group2)
  
  poolsd <-  sqrt( ( s1^2*(n1-1) + s2^2*(n2-1) ) / (n1 + n2-2) )  # calculate the pooled standard deviation  
  poolsd[poolsd ==0] <- 1           ### if the sd is 0, replace with 1. 
  
  new_data <- data.frame(matrix(NA, nrow(microarray), ncol(microarray)))
 
  new_data <- microarray/poolsd
 
  rownames(new_data) <- rownames(microarray); 
  colnames(new_data) <- colnames(microarray)
  
  return(new_data)
  }


estimate.sigma <- function(microarray, trt){
  
  microarray <- as.matrix(microarray)
  group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  samp_rho <- cor(t(resid_mat))                 

# two lines correction for inflated type I error.	  
#   nsamp <- ncol(microarray)
#  samp_rho <- cov(t(resid_mat)) * (nsamp-1)/(nsamp-2)

  factors <- unique(trt)
  m1 <- sum(trt == factors[1])
  m2 <- sum(trt == factors[2]) 
  sigma2 <- 1/m1 + 1/m2

  ndim <- nrow(microarray)
  
  size <- length(trt)
  idex1 <- which(trt == factors[1])[1]
  idex2 <- which(trt == factors[2])[1]
  t_val <- group_mean[, idex2] - group_mean[, idex1]            # control - treatment

  sigma0 <- var(t_val) - sigma2
  sigma0 <- max(0, sigma0)                                     # in case of negative sigma0
  sigma <-  sigma0*diag(ndim) + sigma2*samp_rho
  
  
  return(list(sigma = sigma, t_val = t_val))

}
  

#####  our test
MOM_test <- function(microarray, trt, go_term, standardize=T){
##  we need to standardize it here.

  if (standardize == T){             # do the standardization 
    	microarray <- standardize.microarray(microarray, trt)
  }

   	est_sigma <- estimate.sigma(microarray, trt)
   	sigma <- est_sigma$sigma
   	t_val <- est_sigma$t_val

   	beta0 <- mean(t_val)
   	ones <- rep(1, length(t_val))
   
   	numer <- (crossprod(go_term, (t_val - beta0*ones)))^2
   	go_bar <- rep(mean(go_term), length(go_term))
   	denom2 <- t(go_term-go_bar) %*% sigma %*% (go_term- go_bar)

   	test_stat <- drop(numer/denom2)
   	pval <- 1 - pchisq(test_stat, 1)	

	###  March 16:  one sided test modification
	t_temp <- mean(t_val[go_term == 1]) - mean(t_val[go_term ==0])
	
	if (t_temp > 0 ){
	test_stat2 <- sqrt(test_stat)	
	pval2 <- 1- pnorm(test_stat2)
	status <- "up"
	}

	else {
	test_stat2 <- sqrt(test_stat)*(-1)
	pval2 <- pnorm(test_stat2)
	status  <- "down"
	}  

   
   return(list(stat = test_stat, p.two = pval, status = status, p.one = pval2))
}
  



  


  

BwtGeneCorr <- function(expression, trt, geneset, standardize = T, minSetSize = 5){
## calculate the mean correlation for testset genes, the background set genes and the inter-set correlations
  
# for GSE64810
# 	microarray <- expression[, -(1:2)]
#  	all_genes <-  expression[, 2]
  
# for other typical data, where row names are genes
	microarray <- expression
	all_genes <- rownames(microarray)
  
  if (standardize == T){             # do the standardization
    	microarray <- standardize.microarray(microarray, trt)
  }
  
  #	## calcuate the sample correlations
  	microarray <- as.matrix(microarray)
  	group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
  	resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
  	samp_rho <- cor(t(resid_mat))
  
  	keep_term <- which(geneset$size >= minSetSize)
  
  	set_name <- set_size <- c()
  	testSet_cor <- backSet_cor <- interCor <- sumTestCor <- rep(NA, length(keep_term))
  
  for ( i in 1:length(keep_term)){
    	gset1 <- geneset$geneSet[[keep_term[i]]][-(1:2)]
    	set_name[i] <- geneset$geneSet[[keep_term[i]]][1]
    	go_term <- ifelse(all_genes %in% gset1, 1, 0)
    	set_size[i] <- sum(go_term)
    
    	go_ind <- which(go_term==1)
    	sumTestCor[i] <- sum(samp_rho[go_ind, go_ind])
    	testSet_cor[i] <- (sumTestCor[i] - set_size[i])/( set_size[i]*(set_size[i]-1) )
    	interCor[i] <- mean(samp_rho[-go_ind, go_ind])
    
  	 # print(c(set_size[i], sumTestCor[i], interCor[i], testSet_cor[i])
  }

  # the correlations for the background genes
  	sumAll <- sum(samp_rho); ndim <- nrow(samp_rho)
  	n2 <- ndim - set_size
  	backSet_cor <- ( sumAll - (interCor*n2*2 + sumTestCor) - n2 ) / (n2^2 - n2)
  
  	cor.matrix <- data.frame(set.name = set_name, testSetCor = testSet_cor, interCor= interCor, backSetCor = backSet_cor)
   
  return(cor.matrix)
}






MOM_test_multiple <- function(expression, trt, geneset, standardize = T, minSetSize = 5){
## conduct enrichment analysis for a battery of gene sets and return the raw p values 
## as well as adjusted p values
##  input: 
##      microarray: the raw expression matrix
##             trt: the treatment lables
##         geneset: object from readGeneSet function
##          

############### for GSE64810 USE the following two lines. Otherwise comment it
# for GSE64810
#	microarray <- expression[, -(1:2)]
#	all_genes <-  expression[, 2]

# for other typical data, where row names are genes	
	microarray <- expression
	all_genes <- rownames(microarray)	

   if (standardize == T){             # do the standardization 
    	microarray <- standardize.microarray(microarray, trt)
    }
  
    	est_sigma <- estimate.sigma(microarray, trt)
    	sigma <- est_sigma$sigma
    	t_val <- est_sigma$t_val


    	beta0 <- mean(t_val)
    	ones <- rep(1, length(t_val))
      
    	keep_term <- which(geneset$size >= minSetSize)
    	print(paste("number of sets to be tested:", length(keep_term)))
    
#    	pval <- c()
	pval <- pval2 <- status <- c() 
   	set_size <- c()      
    	set_name <- c()

	cat("\n this might take a long time...  \n")    
	cat("\n testing gene set: \n")

   for ( i in 1:length(keep_term)){
	gset1 <- geneset$geneSet[[keep_term[i]]][-(1:2)]
	cat('\r', i)      
	set_name[i] <- geneset$geneSet[[keep_term[i]]][1]
	go_term <- ifelse(all_genes %in% gset1, 1, 0)
	set_size[i] <- sum(go_term) 
     		
      	numer <- (crossprod(go_term, (t_val - beta0*ones)))^2
      	go_bar <- rep(mean(go_term), length(go_term))
      	denom2 <- t(go_term-go_bar) %*% sigma %*% (go_term- go_bar)
      
      	test_stat <- drop(numer/denom2)
      	pval[i] <- 1 - pchisq(test_stat, 1)
	
	
	  ###  March 16:  one sided test modification
        t_temp <- mean(t_val[go_term == 1]) - mean(t_val[go_term ==0])

	if(! is.finite(t_temp)) { 
		pval2[i] <- NA
		status[i] <- "Unknown"
	}  else {
		   if (t_temp > 0 ){
                test_stat2 <- sqrt(test_stat)
                pval2[i] <- 1- pnorm(test_stat2)
                status[i] <- "up"
       		} else {

                test_stat2 <- sqrt(test_stat)*(-1)
                pval2[i] <- pnorm(test_stat2)
                status[i]  <- "down"
           	}

	}


    }  

	cat("\n")
	p_fdr <- p.adjust(pval, method = "BH")

	results <- data.frame(set.name = set_name, set.size = set_size, p.two=pval,
			p.one = pval2, status = status, p.fdr = p_fdr)
	
	## calculate the between gene correlations
	btw <- BwtGeneCorr(expression, trt, geneset, minSetSize=2, standardize = T)   
 	result_comb <- merge(results, btw, by = "set.name")	
	
	return(result_comb)
}


Camera_multiple <- function(expression, trt, geneset, use.rank = F){
##  Use Camera procedure to do a battery of gene set test. 

# for GSE64810 data
	microarray <- expression[, -(1:2)]
	all_genes <-  expression[, 2]

# for other typical data, where row names are genes
#        microarray <- expression
#        all_genes <- rownames(microarray)

	gset1 <- list()
	set.name <- c()
	for ( i in 1:length(geneset$size)){
                gset1[[i]] <- geneset$geneSet[[i]][-(1:2)]
	  	set.name[i] <- geneset$geneSet[[i]][1]
	}
	names(gset1) <- set.name
	c2.indices <- ids2indices(gset1, all_genes)      # it contains multiple lists

	design <- model.matrix(~trt)
  	Results <- camera(microarray, c2.indices,  design, use.ranks = use.rank, sort = F)
	return(Results)

}


MRSGE_multiple <- function(expression, trt, geneset, use.rank = T){


# for GSE64810 data
#        microarray <- expression[, -(1:2)]
#        all_genes <-  expression[, 2]

# for other typical data, where row names are genes
        microarray <- expression
        all_genes <- rownames(microarray)

        gset1 <- list()
        set.name <- c()
        for ( i in 1:length(geneset$size)){
                gset1[[i]] <- geneset$geneSet[[i]][-(1:2)]
                set.name[i] <- geneset$geneSet[[i]][1]
        }
	names(gset1) <- set.name
        c2.indices <- ids2indices(gset1, all_genes)	 # it contains multiple lists 
	
	design <- model.matrix(~trt) 
  	fit <- lmFit(microarray, design)
  	fit <- eBayes(fit)								                # Emperical Bayes t test	
  	stat <- fit$t[, 2]            							      #	use the moderated t statistics to do enrichment
  	alter <- "mixed"                     

	n_genes <- nrow(microarray)
	set.name <- names(c2.indices)	
	set.size <- p.MRGSE <- rep(0, length(set.name))
	for ( i in 1:length(c2.indices)) {
		index1 <- c2.indices[[i]]
	# print(i)
	 	set.size[i] <- length(index1)
		p.MRGSE[i] <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= use.rank)  
	}

	result <- data.frame(set.name, size = set.size, p.MRGSE)
	return(result)
}

compare_test <- function(dat){
########  a function to incorporate all test procedures

  library(limma)

   	microarray <- dat$data
  	trt <- dat$trt
  	go_term <- dat$go_term

 #  	pval1 <- MOM_test(microarray, trt, go_term, standardize=F)$p             # our test, NO standardization for the simulation
  	
## modified on March 16  
  	MEQ <- MOM_test(microarray, trt, go_term, standardize=F)             # our test, NO standardization for the simulation
#	print(MEQ)
  	pval1 <- MEQ$p.two							# two sided test
	pval1.2 <- MEQ$p.one							# one sided test
	
  	est_sigma <- estimate.sigma(microarray, trt)
  	t_val <- est_sigma$t_val
  	pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression
  	

  	design <- model.matrix(~trt) 
  	fit <- lmFit(microarray, design)
  	fit <- eBayes(fit)								                # Emperical Bayes t test	
  	stat <- fit$t[, 2]            							      #	use the moderated t statistics to do enrichment

  	alter <- "mixed"                                  # This is the default option, see below for explanation.
  	index1 <- which(go_term==1)                       # which genes are in GOTERM
  
  # up: positive statistics
  # down: negative stats
  # either: the set is either up or down-regulated as a whole
  # mixed: tests whether the genes in the set tend to be DE without regard for direction.
  #        In this case, the test will be significant if the set contains mostly large statistics, negative or positive
  
  	tes1 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= F)         # moderated t GeneSet
  	tes2 <- geneSetTest(index = index1, stat, alternative = alter, ranks.only= T)         # moderated t GeneSet Rank
  	tes3 <- camera(microarray, index1,  design)$PValue                                    # camera proedure
  	tes4 <- camera(microarray, index1,  design, use.ranks= T)$PValue                      # camera rank 

 	tes5 <- GSEA.SingleSet(dat$data, dat$trt, dat$go_term, nperm=1000)$p.vals		# GSEA


  ## qusage <Yaari, 2013>
  	geneSets <- list()
  	geneSets[[paste("Set",1)]] <- which(go_term == 1)
  	labels <- rep(NA, length(trt))
  	labels[trt == 1] <- "B"; labels[trt==0] <- "A"
  	qsarray <- qusage(microarray, labels, contrast = "B-A" , geneSets)  			# calculate the probability for all genes
  	tes6 <- pdf.pVal(qsarray, alternative = "two.sided", selfContained = F) 		# competitive test
  
  
  
#  	return(c(pval1, pval2, tes1, tes2, tes3, tes4, tes5, tes6))

	return(c(pval1, pval1.2, pval2, tes1, tes2, tes3, tes4, tes5, tes6))
  
}






