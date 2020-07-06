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


###########  perform data analysis 
compare_padog_data <- function(dat, seed = 20200702){
  
  disease <- dat$disease
  n_normal <- dat$n_normal
  n_disease <- dat$n_disease
  target_geneset <- dat$target_geneset
  n_gene <- sum(dat$go_term)
  
  # perform enrichment analysis
  p_vec <- compare_test(dat, seed = seed)

  # Modified on May 14, 2017
  # the PLAGE method is implemented in GSVA, and the way to calculate p value is based on examples in
  # https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.pdf
  p_plage <- plage(expression_data = dat$data, trt = dat$trt, go_term = dat$go_term,  
                   seed = seed, nperm = 999)
  # modified on May 14, 2017
  # the over-representation method
  p_ora <- ora(expression_data = dat$data, trt = dat$trt, go_term = dat$go_term, 
               method = "BH", thresh = 0.01) 
  
  p_vec <- p_vec %>% mutate(PLAGE = p_plage, 
                            ORA = p_ora, 
                            disease = disease,
                            n_normal = n_normal, 
                            n_disease = n_disease,
                            target_geneset = target_geneset,
                            gene_setsize = n_gene)
  
  return(p_vec)  
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
setlist3 <- setlist2[-c(2, 6:8, 11, 13, 14, 16:18)]
setlist <- c(setlist1, setlist3)
# set <- "GSE20153"

source('R-code-paper/other-methods.R')
source('R-code-paper/compare-methods.R')
source('R-code-paper/GSEA.1.0.R')

library(tidyverse)
# result <- prep_padog_data(set)
# r5 <- compare_padog_data(result)

analysis_padog_data <- NULL
for(i in 1:length(setlist)){
  cat(setlist[i])
  padog_dat <- prep_padog_data(setlist[i])
  r1 <- compare_padog_data(dat = padog_dat, seed = i + 1000)
  r1 <- r1 %>% mutate(dataset = setlist[i]) %>% 
    select(dataset, meaca, meaca_n, MRGSE, SigPathway, CAMERA_ModT, CAMERA_Rank,
           GSEA, QuSAGE, PLAGE, ORA, testSetCor, interCor, backSetCor, disease, 
           n_normal, n_disease, target_geneset, gene_setsize, status)
  analysis_padog_data <- dplyr::bind_rows(analysis_padog_data, r1)
}

write.csv(analysis_padog_data, "real-data/padog-data-analysis-2.csv", row.names = F)



###########################  BELOW NOT RUN #####################################

## translate it into probe test
prep_padog_data_all_genesets <- function(set, gslist){
  
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
  # targetGeneSets <- Biobase::notes(exp)$targetGeneSets 
  
  #get rid of duplicates by choosing the probe(set) with lowest p-value; get ENTREZIDs for probes
  aT2 <- PADOG::filteranot(esetm = dat_m, group = ano$Group, 
                           paired = (design == "Paired"), block = ano$Block,annotation) 
  
  n_normal <- sum(ano$Group == "c")
  n_disease <- sum(ano$Group == "d")
  trt <- ifelse(ano$Group == "c", 0, 1)
  
  #filtered expression matrix
  dat_m <- as.data.frame(dat_m)
  dat_m$ID <- rownames(dat_m)
  # merge the probeID 
  data <- dplyr::right_join(dat_m %>% arrange(ID), 
                            aT2 %>% arrange(ID), by = "ID")
  dat_out <-  data[, 1:(length(ano$Group))]
  rownames(dat_out) <- data$ID
  
  gs_names <- names(gslist)
  go_terms <- data.frame(matrix(NA, nrow = nrow(data), ncol = length(gs_names)))
  names(go_terms) <- gs_names
  gene_seq <- 1:length(gs_names)
  gene_seq <- gene_seq[-c(104, 105, 199, 200, 204, 205)]
  for (i in gene_seq){
    pks <- KEGGREST::keggGet(gs_names[i])
    genes <- pks[[1]]$GENE
    # print(i)
    if(!is.null(genes)){
      gene1 <- genes[seq(1, length(genes), by = 2)]
      go_terms[, i]  <- ifelse(data$ENTREZID %in% gene1, 1, 0)
    }
  }
  # print(dim(aT2))
  # g1 <- lapply(gslist, function(x) length(x)) %>% unlist()
  
  ## keep only those sets that have at least 5 genes
  go_term_set <- colSums(go_terms)
  go_term <- go_terms[, !is.na(go_term_set) & go_term_set >= 5]
  
  
  result <- list(data = dat_out, 
                 trt = trt, 
                 go_term = go_term,
                 disease = disease, 
                 n_normal = n_normal, 
                 n_disease = n_disease)
  
}

### get all the pathways
library(KEGG.db)
pw2id <- as.list(KEGGPATHNAME2ID)
pw2id = as.list(KEGGPATHID2EXTID)
organism <- "hsa"
gslist = pw2id[grep(organism, names(pw2id))]

set <- "GSE20153"
r1 <- prep_padog_data_all_genesets(set = set, gslist = gslist)
n_geneset <- r1$go_term

gse20153_all_geneset <- NULL
for(i in 1:ncol(n_geneset)){
  geneset_name <- names(n_geneset)[i]
  cat(geneset_name)
  temp0 <- r1
  temp0$go_term <- n_geneset[, i]
  temp <- compare_test(dat = temp0, seed = i + 2000)
  temp <- temp %>% mutate(target_geneset = geneset_name, 
                      gene_setsize = sum(n_geneset[, i]), 
                      n_normal = temp0$n_normal, 
                      n_disease = temp0$n_disease) 
  gse20153_all_geneset <- dplyr::bind_rows(gse20153_all_geneset, temp)
}

write.csv(gse20153_all_geneset, "real-data/gse20153_all_geneset-3.csv", row.names = F)

