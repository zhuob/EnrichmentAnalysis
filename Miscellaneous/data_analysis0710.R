

trt_mean <- function(data, trt){
  # a function to calculate treatment means
  # return a matrix
  id1 <- which(trt == 1)
  nrep1 <- length(id1)
  nrep2 <- ncol(data)-nrep1
  mean1 <- apply(data[, id1], 1, mean)
  mean1 <- matrix(rep(mean1, nrep1), ncol=nrep1, byrow=F) 
  mean2 <- apply(data[, -id1], 1, mean)
  mean2 <- matrix(rep(mean2, nrep2), ncol=nrep2, byrow=F) 
  
  data_mean <- data.frame(matrix(NA, nrow(data), ncol(data)))
  data_mean[, id1] <- mean1
  data_mean[, -id1] <- mean2
  colnames(data_mean) <- colnames(data)
  rownames(data_mean) <- rownames(data)
  data_mean
}



btw_gene_corr1 <- function(expression_data, trt, go_term, standardize = T){
  
  if (standardize == T){             # do the standardization
    expression_data <- standardize_expression_data(expression_data, trt)
  }
  
  #	## calcuate the sample correlations
  expression_data <- as.matrix(expression_data)
  # calculate correlation matrix
  group_mean <- as.matrix(trt_mean(expression_data, trt))     
  # the trt effects are removed from matrix
  resid_mat <- expression_data - group_mean                     
  
  samp_rho <- cov(t(resid_mat))
  set_name <- set_size <- c()
  testSet_cor <- backSet_cor <- interCor <- sumTestCor <- NA
  
  set_size <- sum(go_term)
  
  go_ind <- which(go_term==1)    # which are go_terms
  sumTestCor <- sum(samp_rho[go_ind, go_ind])
  testSet_cor <- (sumTestCor - set_size)/( set_size*(set_size-1) )
  interCor <- mean(samp_rho[-go_ind, go_ind])
  
  # the correlations for the background genes
  sumAll <- sum(samp_rho); ndim <- nrow(samp_rho)
  n2 <- ndim - set_size
  backSet_cor <- ( sumAll - (interCor*n2*2 + sumTestCor) - n2 ) / (n2^2 - n2)
  
  cor_matrix <- data.frame(testSetCor = testSet_cor, interCor= interCor, backSetCor = backSet_cor)
  
  return(cor_matrix)
}


#' The PLAGE method 
#' 
#' @param dat  the simulated data (returned from \code{link{simulate_expression_data}}) consisting of three lists:
#'   \item{data}{the expression matrix}
#'   \item{trt}{what treatment each column of \code{data} belongs to}
#'   \item{go_term}{whether each row of \code{data} belongs to this go term}
#' @param nperm   Number of permutation to run (to get the p value by permuted results t statistics [see equation 4 of the paper]  )
#' @return pval 
#' @export
#'   

plage <- function(dat, nperm = 999){
  
  expression_data <- dat$data
  rownames(expression_data) <- paste("Gene", 1:nrow(expression_data), sep = "")
  go_term <- dat$go_term
  trt <- dat$trt
  gs <- list( set1 =  paste("Gene", which(go_term==1), sep = ""))
  
  ## calculate the activity level for each sample
  es_plage <- plage_score(expression_data, gs)
  # get the t-statistic for each treatment
  es_a <- es_plage[trt== 1]; 
  es_b <- es_plage[trt==0] 
  na <- length(es_a); nb <- length(es_b); n_total <- na + nb
  t_obs <- abs(mean(es_a)- mean(es_b))/sqrt(var(es_a)/na + var(es_b)/nb)
  
  t_perm <- rep(NA, nperm)
  # do the permutations
  for (k in 1:nperm){
    na_perm <- sample(1:n_total, na)
    nb_perm <- setdiff(1:n_total, na_perm)
    expression_data_perm <- expression_data[,c(na_perm, nb_perm)]
    es_plage_perm <- plage_score(expression_data_perm, gs)
    es_a_perm <- es_plage_perm[trt== 1]; 
    es_b_perm <- es_plage_perm[trt==0] 
    t_perm[k] <- abs(mean(es_a_perm)- mean(es_b_perm))/sqrt(var(es_a_perm)/na + var(es_b_perm)/nb)
  }
  
  # the pvalue is (b + 1)/(nperm + 1)
  b <- sum(t_perm >= t_obs)
  pval <- (b + 1)/(nperm + 1)
  return(pval)
}


###    this code is extracted from GSVA package 
plage_score <- function(X, gs){
  geneSets <- gs$set1
  p <- nrow(X)
  n <- ncol(X)
  
  Z <- t(apply(X, 1, function(x) (x-mean(x))/sd(x)))
  
  rightsingularsvdvectorgset <- function(gSetIdx, Z) {
    s <- svd(Z[gSetIdx, ])
    s$v[, 1]
  }
  
  es<-  rightsingularsvdvectorgset(geneSets, Z)
  names(es) <- colnames(X)
  
  return(es)
}






#' the over-representation method (see http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079217 ). 
#' 
#' @param  fit a object returned from \code{eBayes} or \code{lmFit} in \code{limma} package.
#' @param method The method used to adjust for nominal p values if necessary. \code{BH} by default
#' @param thresh the threshold described in the paper, 200 by default.
#' @return pval  the enrichment p value.
#' 
ora <- function(fit, go_term, method = "BH", thresh = 0.01){
  # the implementation of ORA method is described in this paper
  # http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0079217
  
  fit_result <- topTable(fit, number = length(fit$p.value), sort.by = "none", adjust.method = method)
  fold_change <- abs(fit_result$logFC)  ## the absolute fold change
  pvals <- fit_result$P.Value
  adjust_p <- fit_result$adj.P.Val
  
  thresh1 <- thresh
  if (thresh < 1 ) {
    thresh1 <- ceiling(thresh*length(adjust_p))
  }
  
  # procedure 1
  if (sum(adjust_p < 0.1)> thresh1) {
    DE <- adjust_p < 0.1
    cont <- ftable(data.frame(go_term, DE))
  }
  else {
    if (sum(pvals < 0.05 & fold_change > 1.5) > thresh1){ ## procedure 2
      DE <- pvals < 0.05 & fold_change > 1.5
      cont <- ftable(data.frame(go_term, DE))
    } else {    # procedure 3
      cutoff <- quantile(pvals, 0.01)
      DE <- pvals <= cutoff
      cont <- ftable(data.frame(go_term, DE))
    }
  }
  
  k <- colSums(cont)[2]  # the number of DE genes
  ros <- rowSums(cont)
  m <- ros[2]  # the number of GO term genes
  n <- ros[1]  # the number of background genes
  q <- cont[2, 2] # the number of DE genes from the GO term 
  
  # the probability of getting more than q DE genes
  # p_enrich <- phyper(q = q-1, m = m, n = n, k = k, lower.tail =  F) 
  ## BZ made change on June 11, 2017
  p_enrich <- phyper(q = q, m = m, n = n, k = k, lower.tail =  F)  # P(X > q)
  
  return(p_enrich)
}


compare_test1 <- function(dat){
  ########  a function to incorporate all test procedures
  
  library(limma)   # for CAMERA methods
  library(qusage) # for QuSAGE method
  #library(GSVA)
  
  expression_data <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term
  disease <- dat$disease
  n_normal <- dat$n_normal
  n_disease <- dat$n_disease
  
  # our test, NO standardization for the simulation
  #  	pval1 <- MOM_test(microarray, trt, go_term, standardize=F)$p             
  
  ## modified on March 16  
  # our test, NO standardization for the simulation
  MEQ <- meaca::meaca_single(expression_data, trt, go_term, standardize=T)             
  #	print(MEQ)
  p_meaca <- MEQ$p1							# chi-square test
  pval1_2 <- MEQ$p2							# normal reference distribution
  
  # calculate the three rho's
  cor <- btw_gene_corr1(expression_data, trt, go_term, standardize = T)
  
  # est_sigma <- estimate_sigma(expression_data, trt)
  # t_val <- est_sigma$t_val
  # pval2 <- summary(lm(t_val~ go_term))$coe[2, 4]            # linear regression
  
  
  design <- model.matrix(~trt) 
  fit <- limma::lmFit(expression_data, design)
  #stat <- fit$coefficients[, 2]
  # Emperical Bayes t test	
  fit <- limma::eBayes(fit)								                
  #	use the moderated t statistics to do enrichment
  stat <- fit$t[, 2]            							      
  
  
  # modified on May 14, 2017
  # the over-representation method
  p_ora <- ora(fit, go_term, method = "BH", thresh = 0.01) 
  
  # modified on March 18, 2016
  # This is the default option, see below for explanation.
  #  alter <- "mixed"                                  
  # this is comparable to what we are testing
  alternative <- "either"                            
  # which genes are in GOTERM
  index1 <- which(go_term==1)                       
  
  # up: positive statistics
  # down: negative stats
  # either: the set is either up or down-regulated as a whole
  # mixed: tests whether the genes in the set tend to be DE without regard for direction.
  #        In this case, the test will be significant if the set contains mostly large statistics, negative or positive
  
  # MRGSE
  p_mrgse <- limma::geneSetTest(index = index1, stat, alternative = alternative, ranks.only= T)         
  # sigPathway methods
  p_sigpath <- meaca::sig_path(index = index1, stat, nsim = 999)                                                
  # camera proedure
  p_camera <- limma::camera(expression_data, index1,  design)$PValue                                    
  # camera rank 
  p_camera_R <- limma::camera(expression_data, index1,  design, use.ranks= T)$PValue                      
  
  # GSEA
  p_gsea <- meaca:::GSEA.SingleSet(dat$data, dat$trt, dat$go_term, nperm=1000)$p.vals		
  
  ## qusage <Yaari, 2013>
  geneSets <- list()
  geneSets[[paste("Set",1)]] <- which(go_term == 1)
  labels <- rep(NA, length(trt))
  labels[trt == 1] <- "B"; labels[trt==0] <- "A"
  # calculate the probability for all genes
  qsarray <- qusage::qusage(expression_data, labels, contrast = "B-A" , geneSets)  			
  # competitive test
  p_qusage <- qusage::pdf.pVal(qsarray, alternative = "two.sided", selfContained = F) 		
  
  
  # Modified on May 14, 2017
  # the PLAGE method is implemented in GSVA, and the way to calculate p value is based on examples in
  # https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
  p_plage <- plage(dat = dat, nperm = 999)
  
  ## Modified on May 17, 2017
  ## the GSVA method
  # gs <- list(set1 = rownames(expression_data)[geneSets$`Set 1`])
  # es_gsva <- gsva(expression_data, gs, method = "gsva",
  #                     verbose=FALSE, parallel.sz=1)$es.obs
  # p_gsva <- t.test(es_gsva[trt==1], es_gsva[trt==0])$p.value
  
  p_vec <- c(p_meaca, pval1_2, p_plage, p_mrgse, p_sigpath, p_camera, p_camera_R,
             p_gsea, p_qusage, p_ora, cor)
  names(p_vec)[1:10] <- c("meaca", "meaca_n", "PLAGE", "MRGSE", "SigPathway", 
                    "CAMERA_ModT", "CAMERA_Rank", "GSEA", "QuSAGE", "ORA")
  p_vec <- data.frame(p_vec)
  p_vec$disease <- disease
  p_vec$n_normal <- n_normal
  p_vec$n_disease <- n_disease
  return(p_vec)
  
  
}






## making use of the data set and KEGG package



prep_data <- function(set){
  
  # the database to load
  library(KEGGdzPathwaysGEO)
  library(KEGGREST)
  library(hgu133plus2.db)
  
  data(list = set)
  x <- get(set)
  exp <- experimentData(x)
  dataset <- exp@name
  disease <- notes(exp)$disease
  dat_m <- exprs(x)  # expression matrix
  ano <- pData(x)
  design <- notes(exp)$design
  stopifnot(trimws(toupper(design)) == "NOT PAIRED")
  annotation <- paste(x@annotation, ".dp", sep = "")
  # the name of the gene set
  targetGeneSets <- notes(exp)$targetGeneSets  
  
  # the members of the gene set
  pks <- KEGGREST::keggGet(paste("hsa", targetGeneSets, sep = ""))
  genes <- pks[[1]]$GENE
  gene1 <- genes[seq(1, length(genes), by = 2)]
  
  n_normal <- sum(ano$Group == "c")
  n_disease <- sum(ano$Group == "d")
  
  anpack <- "hgu133plus2ENTREZID"
  t1 <- keys(get(anpack))
  ## annotate ids to ENTREZID
  ENTREZID <- unlist(as.list(get(anpack)))
  entire_list <- data.frame(probeID = names(ENTREZID), ENTREZID)
  entire_list <- data.frame(lapply(entire_list, as.character), stringsAsFactors=FALSE)
  entire_list <- entire_list[!is.na(entire_list$ENTREZID), ]
  entire_list <- entire_list[!duplicated(entire_list$ENTREZID), ]
  
  ## match the gene ids to those in the data set
  entire_list$go_term <- ifelse(entire_list$ENTREZID %in% gene1, 1, 0)
  
  trt <- ifelse(ano$Group == "c", 0, 1)

  #filtered expression matrix
  dat_m <- as.data.frame(dat_m)
  dat_m$probeID = rownames(dat_m)
  # merge the probeID 
  data <- dplyr::right_join(dat_m %>% arrange(probeID), 
                            entire_list %>% arrange(probeID), by = "probeID")
  
  dat_out <-  data[, 1:(length(ano$Group))]
  rownames(dat_out) <- data$probeID
  
  result <- list(data = dat_out, 
                 trt = trt, 
                 go_term = data$go_term,
                 disease = disease, n_normal = n_normal, n_disease = n_disease)
  
  return(result)
}


prep_data2 <- function(set){
  data(list=set)
  x=get(set)
  #Extract from the dataset the required info
  exp=experimentData(x);
  dataset= exp@name
  dat_m=exprs(x)
  disease <- notes(exp)$disease
  ano=pData(x)
  design= notes(exp)$design
  annotation= paste(x@annotation,".db",sep="")
  targetGeneSets <- notes(exp)$targetGeneSets 
  
  #get rid of duplicates in the same way as is done for PADOG and assign probesets to ENTREZ IDS
  #get rid of duplicates by choosing the probe(set) with lowest p-value; get ENTREZIDs for probes
  aT2=filteranot(esetm=dat_m,group=ano$Group,paired=(design=="Paired"),block=ano$Block,annotation) 
  

 
  pks <- KEGGREST::keggGet(paste("hsa", targetGeneSets, sep = ""))
  genes <- pks[[1]]$GENE
  gene1 <- genes[seq(1, length(genes), by = 2)]
  # print(dim(aT2))
  aT2$go_term  <- ifelse(aT2$ENTREZID %in% gene1, 1, 0)
  n_normal <- sum(ano$Group == "c")
  n_disease <- sum(ano$Group == "d")
  
  trt <- ifelse(ano$Group == "c", 0, 1)
  #filtered expression matrix
  dat_m <- as.data.frame(dat_m)
  dat_m$ID = rownames(dat_m)
  # merge the probeID 
  data <- dplyr::right_join(dat_m %>% arrange(ID), 
                     aT2 %>% arrange(ID), by = "ID")
  
  dat_out <-  data[, 1:(length(ano$Group))]
  rownames(dat_out) <- data$ID
  
  result <- list(data = dat_out, 
                 trt = trt, 
                 go_term = data$go_term,
                 disease = disease, n_normal = n_normal, n_disease = n_disease)
  
}





library(KEGGdzPathwaysGEO)
# Not found: GSE5281, GSE6956
# PAIRED: GSE8671, GSE15471, GSE16515, GSE3467, GSE3678, GSE18842
setlist1 <- c("GSE20153", "GSE20291", "GSE8762", "GSE4107", "GSE14762",
             "GSE1297", "GSE9348", "GSE781",  
             "GSE19728", "GSE21354",  'GSE9476', 
             "GSE19188", "GSE3585")




library(KEGGandMetacoreDzPathwaysGEO)
setlist2 <- c("GSE1145", "GSE11906", "GSE14924_CD4", "GSE14924_CD8", 
              "GSE16759", 'GSE19420', 'GSE20164', "GSE23878", "GSE24739_G0", 
              "GSE24739_G1", "GSE30153", "GSE32676", "GSE38666_epithelia", 
              "GSE38666_stroma", "GSE4183", "GSE42057", "GSE7305", "GSE22780")
# Bad Request (HTTP 400): 2, 6, 7, 11, 13, 14, 16
# Paird, 8, 17, 18
set <- "GSE4183"
library(PADOG)
source('~/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/meaca/meaca_for_paper/meaca/R/other-methods.R')
library(tidyverse)
result <- prep_data2(set)
r5 <- compare_test1(result)

setlist3 <- setlist2[-c(2, 6:8, 11, 13, 14, 16:18)]
setlist <- c(setlist1, setlist3)

enrich_result <- data.frame(matrix(NA, ncol = 17, nrow = length(setlist)))
names(enrich_result) <- c("set", names(r5))
enrich_result$set[length(setlist)] <- set
enrich_result[length(setlist), -1] <- r5

for(i in 1:length(setlist)){
  cat(setlist[i])
  result <- prep_data2(setlist[i])
 r1 <- compare_test1(result)
 enrich_result$set[i] <- setlist[i]
 enrich_result[i, -1] <- r1

#  print(r1)
}

write.csv(enrich_result, "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/GB_revision/real_data20_new_v1.csv", row.names = F)


### get all the pathways
library(KEGG.db)
pw2id <- as.list(KEGGPATHNAME2ID)
pw2id = as.list(KEGGPATHID2EXTID)
organism <- "hsa"
gslist = pw2id[grep(organism, names(pw2id))]


set <- "GSE4183"
r1 <- prep_data(set)
r2 <- prep_data2(set)
c1 <- btw_gene_corr1(r1$data, trt = r1$trt, go_term = r1$go_term)
c2 <- btw_gene_corr1(r2$data, trt = r2$trt, go_term = r2$go_term)

