# test "meaca" on real data

library(meaca)
library(limma)
source("R-code-paper/read-gene-set.R")

# these two functions come from MEQLEA, only used to read data of this particular type
CLS <- GSEA.ReadClsFile(file="real-data/Gender/Gender.cls")  # treatment lables
dataset <- GSEA.Gct2Frame2(filename = "real-data/Gender/Gender.gct") # expression data
geneset <- read_gene_set("real-data/Gender/C1.gmt")


## run meaca analysis
result_meaca <- meaca_multiple(expression_data = dataset, trt = CLS$class.v, 
                          geneset = geneset, min_set_size=2)


## run CAMERA analysis
result_camera <- Camera_multiple(expression_data = dataset, trt = CLS$class.v, 
                                 geneset = geneset, use.rank = F)

result_camera_rank <- Camera_multiple(expression_data = dataset, trt = CLS$class.v, 
                                 geneset = geneset, use.rank = T)

result_mrsge <- MRGSE_multiple(expression_data = dataset, trt = CLS$class.v, 
                               geneset = geneset, use.rank = T)

# result1 <- result_camera %>% mutate(set_name = row.names(result_camera)) %>% 
#   arrange(set_name) %>% left_join(result_meaca %>% arrange(set_name))


## GSEA method
GSEA_results <- GSEA(dataset, CLS, 
                     gs.db = "real-data/Gender/C1.gmt",
                     #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
                     doc.string            = "Gender_C1",     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
                     non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
                     reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
                     nperm                 = 9999,            # Number of random permutations (default: 1000)
                     weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
                     nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
                     fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
                     fdr.q.val.threshold   = -1,  #I changed this from 0.25          # Significance threshold for FDR q-vals for gene sets (default: 0.25)
                     topgs                 = 10000,    #I changed this from 20          # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
                     adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
                     gs.size.threshold.min = 2,     # I changed this         # Minimum size (in genes) for database gene sets to be considered (default: 25)
                     gs.size.threshold.max = 10000, # I changed this            # Maximum size (in genes) for database gene sets to be considered (default: 500)
                     reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
                     preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
                     random.seed           = 141,             # Random number generator seed. (default: 123456)
                     perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
                     fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
                     replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
                     save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
                     OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
                     use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
)

result_gsea <- rbind(GSEA_results $report1, GSEA_results $report2)


