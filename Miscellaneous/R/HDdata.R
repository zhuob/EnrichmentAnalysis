
source("/home/stats/zhuob/Rcode/Enrichment/MOMEnrichmentTest.R")
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabDataNew.R")               # it contains group.mean function
source("/home/stats/zhuob/Rcode/Enrichment/GSEA.1.0.R")

setwd("/home/stats/zhuob/data/computing/DataSet/")

expressionData <- read.table("GSE64810_mlhd_DESeq2_norm_counts_adjust.txt", header=T)   # the expression data
 genes <- rownames(expressionData)
 Ensembl.Gene.ID <- substring(genes, 1, which(strsplit(genes, "")[[1]] == ".")-1)
 expressionData <- data.frame(Ensembl.Gene.ID, expressionData)
 HGNC <- read.table("mart_export.txt", sep = ",", header=T)
 newData <- merge(HGNC, expressionData,  by = "Ensembl.Gene.ID")  # the data set

 geneset <- readGeneSet("c2.cp.v5.1.symbols.gmt")  		# the gene sets

#geneset <- readGeneSet("c3.tft.v5.1.symbols.gmt")

trt <- c(rep(0, 49), rep(1, 20))  # the treatment structure, 0 as control, 1 as treatment

OurMethod_results <- MOM_test_multiple(newData, trt, geneset, minSetSize=2)   		## this runs ourMethod

Camera_results <- Camera_multiple(newData, trt, geneset, use.rank=F)
Camera_results$set.name <- rownames(Camera_results)

MRSGE_results <-  MRSGE_multiple(newData, trt, geneset, use.rank = T)

write.csv(OurMethod_results, "C2.OurMethod.csv", row.names= F)
write.csv(Camera_results, "C2.CAMERA.csv", row.names=F)
write.csv(MRSGE_results, "C2.MRSGE.csv", row.names=F)


############ For GSEA Method

CLS <- list()
CLS$phen = c("T", "C")
CLS$class.v <- trt


GSEA_results <- GSEA(newData, CLS, gs.db = "c2.cp.v5.1.symbols.gmt",
            #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
            doc.string            = "Gender_C1",     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
            non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
            reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
            nperm                 = 1000,            # Number of random permutations (default: 1000)
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

saveRDS(GSEA_results, "C2.GSEA.rds")
result <- rbind(GSEA_results $report1, GSEA_results $report2)
write.csv(result, "C2.GSEA.csv", row.names=F)


