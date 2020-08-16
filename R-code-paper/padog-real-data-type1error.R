#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## this file is to run for real-data based type 1 error simulation
##   Method:  use real data, extract test statistics used in different methods
##            use real data to derive between-gene correlation
##            but the gene level statistics for the background gene set is sampled
##            from the test set statitics (of course, with replacement)
##            Therefore, we create a hypothetical null, based on which we 
##            evaluate the model performance on type 1 error
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## 

## making use of the data set and KEGG package
prep_padog_data <- function(set){
  
  data(list = set)
  x = get(set)
  #Extract from the dataset the required info
  exp = Biobase::experimentData(x);
  dataset = exp@name
  dat_m = Biobase::exprs(x)
  disease <- Biobase::notes(exp)$disease
  ano = Biobase::pData(x)
  design = Biobase::notes(exp)$design
  annotation = paste(x@annotation,".db",sep="")
  targetGeneSets <- Biobase::notes(exp)$targetGeneSets 
  
  #get rid of duplicates in the same way as is done for PADOG and assign probesets to ENTREZ IDS
  # anpack <- "hgu133plus2ENTREZID"
  # t1 <- keys(get(anpack))
  # ## annotate ids to ENTREZID
  # ENTREZID <- unlist(as.list(get(anpack)))
  # entire_list <- data.frame(probeID = names(ENTREZID), ENTREZID)
  # entire_list <- data.frame(lapply(entire_list, as.character), stringsAsFactors=FALSE)
  # entire_list <- entire_list[!is.na(entire_list$ENTREZID), ]
  # aT2 <- entire_list[!duplicated(entire_list$ENTREZID), ]
  
  #get rid of duplicates by choosing the probe(set) with lowest p-value; get ENTREZIDs for probes
  aT2 <- PADOG::filteranot(esetm = dat_m, group = ano$Group, 
                           paired = (design == "Paired"), block = ano$Block,annotation) 
  
  pks <- KEGGREST::keggGet(paste("hsa", targetGeneSets, sep = ""))
  genes <- pks[[1]]$GENE
  gene1 <- genes[seq(1, length(genes), by = 2)]
  n_normal <- sum(ano$Group == "c")
  n_disease <- sum(ano$Group == "d")
  
  trt <- ifelse(ano$Group == "c", 0, 1)
  #filtered expression matrix
  dat_m <- as.data.frame(dat_m)
  dat_m$ID <- rownames(dat_m)
  # merge the probeID 
  data <- dplyr::right_join(dat_m %>% arrange(ID), 
                            aT2 %>% arrange(ID), by = "ID")
  # print(dim(aT2))
  go_term  <- ifelse(data$ENTREZID %in% gene1, 1, 0)
  
  dat_out <-  data[, 1:(length(ano$Group))]
  rownames(dat_out) <- data$ID
  
  result <- list(data = dat_out, 
                 trt = trt, 
                 go_term = go_term,
                 disease = disease, 
                 n_normal = n_normal, 
                 n_disease = n_disease,
                 target_geneset = targetGeneSets,
                 n_gene = sum(go_term))
  
}

create_bootstrap_data <- function(expression_data, go_term, trt, seed = 123, 
                                  raw_data = FALSE, no_de = FALSE){
  
  expression_data <- as.matrix(expression_data)
  group_mean <- as.matrix(meaca:::trt_mean(expression_data, trt))
  resid_mat <- expression_data - group_mean
  
  if(no_de){
    n0 <- dim(expression_data)
    group_mean <- matrix(0, nrow = n0[1], ncol = n0[2])
  }
  
  test_genes <- which(go_term == 1)
  back_genes <- which(go_term == 0)
  # t_test_set <- t_val0[test_genes]
  n_back_set <- sum(go_term == 0)
  n_test_set <- sum(go_term == 1)
  
  ## bootstrap beta0, beta1
  
  set.seed(seed)
  resample_test_genes <- base::sample(1:length(go_term), size = n_test_set, replace = TRUE)
  resample_back_genes <- base::sample(1:length(go_term), size = n_back_set, replace = TRUE)
  group_mean_hat <- group_mean # initiate the matrix
  if(!raw_data){
    group_mean_hat[test_genes, ] <- group_mean[resample_test_genes, ]
    group_mean_hat[back_genes, ] <- group_mean[resample_back_genes, ]
  }

  # bootstrap residual independently from beta0 and beta1
  set.seed(seed + 1e5)
  resample_id <- base::sample(1:length(trt), size = length(trt), replace = TRUE)
  resid_mat_hat <- resid_mat
  if(!raw_data){
    resid_mat_hat <- resid_mat[, resample_id]
  }
  expression_data_hat <- group_mean_hat + resid_mat_hat
  
  return(expression_data_hat)
  
}



compare_test_new <- function(dat, seed, nsim = 1e3, ncore = 4, no_de = FALSE, 
                             package_used = c("tidyverse"), verbose_show = FALSE){
  
  ########  a function to incorporate all test procedures
  expression_data <- dat$data
  trt <- dat$trt
  go_term <- dat$go_term

  library(foreach)
  cl <- parallel::makeCluster(ncore, type = "FORK")
  doParallel::registerDoParallel(cl)
  set.seed(seed)
  subseed <- sample(1:1e6, nsim)
  ## add a standardization step 
  expression_data <- meaca::standardize_expression_data(expression_data, trt)
  
  pvals <- foreach(kk = 1:nsim, .combine = bind_rows, 
                   .packages = package_used, 
                   .verbose = verbose_show) %dopar% {
    
    expression_data_hat <- create_bootstrap_data(expression_data = expression_data, 
                                                 go_term = go_term, trt = trt, 
                                                 seed = subseed[kk], raw_data = FALSE, 
                                                 no_de = no_de)
    
    expression_data_hat <- meaca::standardize_expression_data(expression_data_hat, trt)
    
    temp <- meaca::meaca_single(expression_data = expression_data_hat, 
                                trt = trt, go_term = go_term, 
                                standardize = FALSE)
    p_meaca <- temp$p1
    
    design <- model.matrix(~trt) 
    fit <- limma::lmFit(expression_data_hat, design)
    # Emperical Bayes t test	
    fit <- limma::eBayes(fit)								                
    #	use the moderated t statistics to do enrichment
    stat <- fit$t[, 2]            							      

    alternative <- "either"                            
    # which genes are in GOTERM
    index1 <- which(go_term==1)                       
    
    # MRGSE
    p_mrgse <- limma::geneSetTest(index = index1, stat, alternative = alternative, ranks.only= T)         
    # sigPathway methods
    p_sigpath <- sig_path(index = index1, stat, nsim = 999)                                                
    # camera proedure
    p_camera <- limma::camera(y = expression_data_hat, index = index1,  
                              design = design, 
                              allow.neg.cor = TRUE, inter.gene.cor = NA)$PValue                                    
    # camera rank 
    p_camera_R <- limma::camera(y = expression_data_hat, index = index1,  
                                design = design, use.ranks= T, 
                                allow.neg.cor = TRUE, inter.gene.cor = NA)$PValue                      
  
    # GSEA
    p_gsea <- GSEA.SingleSet(data = expression_data, trt = trt, go_term = go_term, nperm = 1000)$p.vals		
     
    ## qusage <Yaari, 2013>
    geneSets <- list()
    geneSets[[paste("Set",1)]] <- which(go_term == 1)
    labels <- rep(NA, length(trt))
    labels[trt == 1] <- "B"; labels[trt==0] <- "A"
    # calculate the probability for all genes
    qsarray <- qusage::qusage(expression_data_hat, labels, contrast = "B-A" , geneSets)  			
    # competitive test
    p_qusage <- qusage::pdf.pVal(qsarray, alternative = "two.sided", selfContained = F) 		
    
    p_plage <- plage(expression_data = expression_data_hat, trt = trt, 
                     go_term = go_term, seed = kk, nperm = 999)
    # modified on May 14, 2017
    # the over-representation method
    p_ora <- ora(expression_data = expression_data_hat, trt = trt, 
                 go_term = go_term, method = "BH", thresh = 0.01) 
    
    temp1 <- tibble::tibble(p_meaca = p_meaca, 
                            p_mrgse = p_mrgse, 
                            p_sigpath = p_sigpath,
                            p_camera = p_camera, 
                            p_camera_R = p_camera_R,
                            p_gsea  = p_gsea,
                            p_qusage = p_qusage, 
                            p_plage = p_plage, 
                            p_ora = p_ora)
    
      return(temp1)
    }
  parallel::stopCluster(cl)
  
  return(pvals)
}

plot_type1error <- function(dat_type = "test", var_col = 1:4){
  dat_name <- paste0("real-data/padog-package/padog-real-data-type1error-simulation-", 
                     dat_type,".csv")
  t1 <- read_csv(dat_name)
  print(dat_name)
  sample_test_genes <- pivot_longer(t1, cols = var_col, names_to = "method")
  
  ggplot(data = sample_test_genes, aes(sample = value)) + 
    stat_qq(distribution = qunif) + 
    geom_abline(intercept = 0, slope = 1) + 
    facet_wrap(~method)
  
}


#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                    RUN Simulation
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
parent_folder <- "/home/stats/zhuob/EnrichmentAnalysis/"
parent_folder <- ""

library(tidyverse)
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)

r_file <- c("R-code-paper/other-methods.R", "R-code-paper/GSEA.1.0.R")
r_file <- paste0(parent_folder, r_file)
for(k in 1:length(r_file)){source(r_file[k])}

df1 <- prep_padog_data("GSE8762")
result <- compare_test_new(dat = df1, seed = 1234, nsim = 1000, ncore = 1, no_de = TRUE)
save_data <- paste0(parent_folder, "real-data/padog-package/padog-real-data-type1error-simulation-all-genes-data-corr-test-bootstrap-separately-de0.csv")
write_csv(result, save_data)

pdf(paste0(parent_folder, "real-data/padog-package/padog-type1error.pdf"), width = 10, height = 10)
print(plot_type1error(dat_type = "all-genes-data-corr-test-bootstrap-separately-de0", var_col = 1:9))
dev.off()

