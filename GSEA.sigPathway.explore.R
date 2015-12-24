### comparison of GSEA code downloaded from Broad Institute

source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/GSEA.1.0.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/LinuxCode/SimulateLabDataNew.R")
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/LinuxCode/MOMEnrichmentTest.R")


data <- read.table("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/GSEASoftware/Gender.txt")
sample.names <- names(data)

x <- "PITX3	SPFH1	NEURL	C10orf12	NDUFB8	C10orf95	DNTT	USMG5	CWF19L1	SUFU	OBFC1	PEO1	PIK3AP1	UBTD1	CUTC	SEC31L2	AS3MT	MGEA5	NPM3	PSD	KAZALD1	FAM26C	C10orf61	NKX2-3	NOLC1	SEMA4G	ARL3	MMS19L	CYP17A1	GBF1	ANKRD2	C10orf33	ALDH18A1	WNT8B	ARHGAP19	PPRC1	TLL2	ZFYVE27	FRAT1	FBXL15	BLOC1S2	GOT1	CNNM2	ENTPD1	KCNIP2	NT5C2	CNNM1	SHFM3	TLX1	INA	C10orf76	LDB1	HPS1	C10orf32	C10orf26	HPSE2	ABCC2	FAM26A	CUEDC2	SCD	BLNK	C10orf65	SLC25A28	FRAT2	TRIM8	EXOSC1	CHUK	MARVELD1	SORBS1	POLL	HIF1AN	SFRP5	C10orf6	AVPI1	CPN1	MRPL43	ZDHHC16	CCNJ	ENTPD7	SFXN3	PGAM1	C10orf77	PDZK7	NFKB2	PDCD11	PKD2L1	SFXN2	FGF8	SH3MD1	TMEM10	KIAA0690	PAX2	LOXL4	BTRC	FAM26B	ELOVL3	C10orf83	TAF5"
#x <- "ALDH7A1	IL13	SEPT8	IRF1	ACSL6	IL4	SLC12A2	PPIC	CSF2	SLC22A5	CSNK1G3	DMXL1	P4HA2	ZNF608	LOX	FTMT	ADAMTS19	IL5	IL3	LMNB1	PRDM6	SNCAIP	HSD17B4	SNX24	PDLIM4	CDC42SE2	AP3S1	KIF3A	HINT1	TNFAIP8	FBN2	SEMA6A"
## NOM p-val:    0.686 for line 1 and 0.7282 for line 2
## glob.p.val:  0.36 for line 1 and 0.438 for line 2

list1 <- strsplit(x, " ")
class.labels <- read.table("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/GSEASoftware/class.labels.txt")
class.labels <- as.vector(class.labels$V1)
gene <- unlist(strsplit(list1[[1]], "\t"))
ids <- which(rownames(data) %in% gene)
go_term <- rep(0, nrow(data))
go_term[ids] <- 1

col.index <- order(class.labels, decreasing=F)
class.labels <- class.labels[col.index]
sample.names <- sample.names[col.index]
A <- data
A2 <- A
cols <- ncol(data)
for (j in 1:cols) {
  A[, j] <- data[, col.index[j]]
}
names(A) <- sample.names


result <- GSEA.SingleSet(A, class.labels, go_term, nperm=1000)
result$p.vals

prepareData <- list()
prepareData$data <- log(A + 1e-6); prepareData$trt <- class.labels; prepareData$go_term = go_term;
#prepareData$data<- standardize.microarray(prepareData$data, prepareData$trt)
compare_test (prepareData)

xm <- estimate.sigma(prepareData$data, prepareData$trt)
mean(xm$t_val[go_term==1])
mean(xm$t_val[go_term==0])
hist(xm$t_val[go_term==1])


###############  This example shows if we don't do standardization, then the estimate will
##  be significantly different, especially when two genes are very different.
size <- 10
x1 <- matrix(NA, 2, size)
x1[1, ] <- rpois(size, 6)
x1[2, ] <- rpois(size, 120)
trt <- rep(c(0, 1), each=size/2)
x0 <- standardize.microarray(x1, trt)
estimate.sigma(x0, trt)
estimate.sigma(x1, trt)
group.mean(x1, trt)
sep.mean(x1[1, ], trt)



cat("\nTesting the funciton based on my understanding\n")

phi <- result$phi
pos.phi <- phi[1, phi[1, ]>=0]
sum(pos.phi>=0)

neg.phi <- phi[1, phi[1, ] < 0]
sum(neg.phi < 0)

pvalue <- sum(neg.phi < result$Obs.ES)/length(neg.phi)
pvalue
result$p.vals

case <- "a"
size <- 50           # number of samples to be simulated
rho <- c(0.1, 0.05, -0.05)   # correlation for case a, e, f
num_gene <- c(500, 100)
prop <- c(0.2, 0.2)
n_gene <- num_gene[1]
delta <- rnorm(n_gene, 2 , 1)

obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
dat <- simulate.microarry(size, obj)

dat_sd <- standardize.microarray(dat$data, dat$trt)
dat$data <- dat_sd

data_exp <- dat$data
rownames(data_exp) <- paste("Gene", 1:nrow(data_exp), sep="_")
class.label <- -1*(dat$trt-1)
gene.label <- rownames(data_exp)


GSEA <- GSEA.GeneRanking(dat$data, class.labels = class.label, 
                              gene.labels = gene.label, nperm = 10)
result <- GSEA.SingleSet(dat$data, dat$trt, dat$go_term, nperm=1000)
result$p.vals 
result$Obs.ES
quantile(result$phi)
compare_test(dat)





# GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute 
#
# R script to run GSEA Analysis of the Gender vs C1 example (cut and paste into R console)

#GSEA.program.location <- "/Users/Bin/Dropbox/Zhuo/Research/Software/GSEA-P-R/GSEA.1.0.R"   #  R source program (change pathname to the rigth location in local machine)
#source(GSEA.program.location, verbose=F, max.deparse.length=9999)
debug(GSEA)
GSEA(                                                                      # Input/Output Files :-------------------------------------------
 input.ds =  "/Users/Bin/Dropbox/Zhuo/Research/Software/GSEA-P-R/Datasets/Gender.gct",               # Input gene expression Affy dataset file in RES or GCT format
 input.cls = "/Users/Bin/Dropbox/Zhuo/Research/Software/GSEA-P-R/Datasets/Gender.cls",               # Input class vector (phenotype) file in CLS format
 gs.db =     "/Users/Bin/Dropbox/Zhuo/Research/Software/GSEA-P-R/GeneSetDatabases/C1.gmt",           # Gene set database in GMT format
 output.directory      = "/Users/Bin/Dropbox/Zhuo/Research/Software/GSEA-P-R/StudyGSEA/",            # Directory where to store output and results (default: "")
 #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
 doc.string            = "Gender_C1",     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
 non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
 reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
 nperm                 = 1000,            # Number of random permutations (default: 1000)
 weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
 nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
 fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
 fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
 topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
 adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
 gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
 gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
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
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
  directory           = "d:/CGP2005/GSEA/GSEA-P-R/Gender_C1/",        # Directory where to store output and results (default: "")
  topgs = 20,                                                           # number of top scoring gene sets used for analysis
  height = 16,
  width = 16
)











############  explore sigPathway   #############
library(sigPathway)

data(MuscleExample)
## Prepare the pathways to analyze
probeID <- rownames(tab)
gsList <- selectGeneSets(G, probeID, 100, 500)
nsim <- 1000
ngroups <- 2
verbose <- TRUE
weightType <- "constant"
methodName <- "NGSk"
npath <- 25
allpathways <- FALSE
annotpkg <- "hgu133a.db"
statV <- calcTStatFast(tab, phenotype, ngroups)$tstat
res.NGSk <- calculate.NGSk(statV, gsList, nsim, verbose)
## Summarize top pathways from NGSk
res.pathways.NGSk <-
  rankPathways.NGSk(res.NGSk, G, gsList, methodName, npath=npath)
print(res.pathways.NGSk)
## Get more information about the probe sets' means and other statistics
## for the top pathway in res.pathways.NGSk
gpsList <-
  getPathwayStatistics.NGSk(statV, probeID, G, res.pathways.NGSk$IndexG,
                            FALSE, annotpkg)



########################
# Well, it seems there's no need to check sigPathway at all: mentioned in Smyth's 2012
# paper, it uses the ordinary t-statistic, and permutes gene lables, which is similar to 
# geneSetTest() in limma package. 
# And I read the method part of the original paper sigPathway. 



microarray <- dat$data
trt <- dat$trt
go_term <- dat$go_term

library(sigPathway)


sigPathway.SingleSet <- function(microarray, trt, go_term, nsim = 1000, ngroups = 2, methodName = "NGSk"){
  ## should use NTk-like statistics becaues we are answering Q1 of that paper
  
  phenotype <- paste("PH", trt, sep="_")
  probID <- rownames(microarray)
  ngroups <- 2
  gsList <- list()
  gsList$nprobes <- c(100, 50)
  ids1 <- which(go_term == 1)
  ids2 <- sample(1:500, 50)
  gsList$indexV <- c(ids1, ids2)
  gsList$indGused <- c(1, 2)
  
  
  statV <- calcTStatFast(microarray, phenotype, ngroups)$tstat
  res.NGSk <- calculate.NGSk(statV, gsList, nsim, verbose=T)
  
  
}
  
  
  































