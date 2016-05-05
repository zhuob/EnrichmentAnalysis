## diagnostic 
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/R/MOMEnrichmentTest.R")                                          # the code 
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/R/SimulateLabDataNew.R")                                          # the code 
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/GSEA.1.0.R")
# source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/summarize.results.R")
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation")

t1 <- read.table("SimulationPower20160203/NODE_c_50.txt")
t2 <- read.table("SimulationPower20160203/DE_f_50_10PCT.txt")
t3 <- read.table("SimulationPower20160211/SizePoint1/Power_a_25VS10PCT.txt")

dat <- t2
par(mfrow=c(2, 2))
plot(dat$OurTest, dat$Camera, pch=20); plot(dat$OurTest, dat$CameraRank, pch=20)
plot(dat$OurTest, dat$GSEA, pch=20); plot(dat$Camera, dat$GSEA, pch=20)


gender <- read.csv("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/gender/combined.csv")
GSE64810 <- read.csv("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/C2.CombinedResults.csv")


dat2 <- GSE64810
par(mfrow=c(2, 3))
plot(dat2$p, dat2$Camera, pch=20); plot(dat2$p, dat2$GSEA, pch=20);plot(dat2$Camera, dat2$GSEA, pch=20)
plot(-log(dat2$p), -log(dat2$Camera), pch=20); plot(-log(dat2$p), -log(dat2$GSEA), pch=20);plot(-log(dat2$Camera), -log(dat2$GSEA), pch=20)
names(dat2)
cor(dat2[, c(7, 9, 10)])



## dig into the difference between our methods and GSEA/CAMERA
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/") # the data set
GSE64810 <- read.csv("C2.combined3.csv")
geneset <- readGeneSet("c2.cp.v5.1.symbols.gmt")
expressionData <- read.table("GSE64810_mlhd_DESeq2_norm_counts_adjust.txt", header=T)   # the expression data
genes <- rownames(expressionData)
Ensembl.Gene.ID <- substring(genes, 1, which(strsplit(genes, "")[[1]] == ".")-1)
expressionData <- data.frame(Ensembl.Gene.ID, expressionData)

HGNC <- read.table("mart_export.txt", sep = ",", header=T)
newData <- merge(HGNC, expressionData,  by = "Ensembl.Gene.ID")


library(limma)
trt <- c(rep(1, 49), rep(0, 20))  # the treatment structure

design <- model.matrix(~trt) 
fit <- lmFit(newData[, -c(1, 2)], design)
fit <- eBayes(fit)		
pval <- fit$p.value[, 2]
t_val <- fit$t[, 2]

## this gene set is enriched by CAMERA and GSEA, but not OurMethod

for ( i in 1:geneset$NumofSets)
  {
    if (geneset$geneSet[[i]][1] == "REACTOME_FGFR4_LIGAND_BINDING_AND_ACTIVATION")
      { print(i)}
   }

## similar are 123, 468, 1194, 1284, or 419, 678 for significant OurMethods, but not GSEA or CAMERA

set1 <- geneset$geneSet[[123]]
ids <- which(newData$HGNC.symbol %in% set1[-(1:2)])


mean(t_val[ids])
mean(t_val[-ids])
t.test(t_val[ids], t_val[-ids])





## ## this shows the problem is in standardization
source("/home/stats/zhuob/Rcode/Enrichment/MOMEnrichmentTest.R")
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabDataNew.R")               # it contains group.mean fun$
source("/home/stats/zhuob/Rcode/Enrichment/GSEA.1.0.R")

setwd("/home/stats/zhuob/data/computing/DataSet/")

expressionData <- read.table("GSE64810_mlhd_DESeq2_norm_counts_adjust.txt", header=T)   # the expression data
genes <- rownames(expressionData)
Ensembl.Gene.ID <- substring(genes, 1, which(strsplit(genes, "")[[1]] == ".")-1)
expressionData <- data.frame(Ensembl.Gene.ID, expressionData)
HGNC <- read.table("mart_export.txt", sep = ",", header=T)
newData <- merge(HGNC, expressionData,  by = "Ensembl.Gene.ID")  # the data set

geneset <- readGeneSet("c2.cp.v5.1.symbols.gmt")               # the gene sets

#geneset <- readGeneSet("c3.tft.v5.1.symbols.gmt")


set1 <- geneset$geneSet[[678]]
ids <- which(newData$HGNC.symbol %in% set1[-(1:2)])
trt <- c(rep(1, 49), rep(0, 20))  # the treatment structure

## before standardization

design <- model.matrix(~trt) 
fit <- lmFit(newData[, -(1:2)], design)
fit <- eBayes(fit)		
pval <- fit$p.value[, 2]
tval_limma <- fit$t[, 2]
t.test(tval_limma[ids], tval_limma[-ids]) # p-value = 8.651e-06

## after standardization
microarray <- standardize.microarray(newData[, -(1:2)], trt)
fit2 <- lmFit(microarray, design)
fit2 <- eBayes(fit2)		
pval <- fit2$p.value[, 2]
tval2_limma <- fit2$t[, 2]
t.test(tval2_limma[ids], tval2_limma[-ids]) # p-value = 0.4511

microarray <- standardize.microarray2(newData[, -(1:2)], trt)
fit3 <- lmFit(microarray, design)
fit3 <- eBayes(fit3)		
pval <- fit3$p.value[, 2]
tval3_limma <- fit3$t[, 2]
t.test(tval3_limma[ids], tval3_limma[-ids]) # p-value = 0.002883


## our test 
go_term <- rep(0, nrow(microarray))
go_term[ids] <- 1
MOM_test(newData[, -(1:2)], trt, go_term) # p-value 0.5558113






## what about simulated data 
## there's no problem becuase the standardization issue is not a problem

testStandardization <- function(case = case, prop = c(0.3, 0.1)){

  nsim <- 1000
  size <- 50           # number of samples to be simulated
  
  # rho <- c(0.7, 0.5, -0.5)   # extreme case
  rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation
  num_gene <- c(500, 100)
  prop <- prop
  n_gene <- num_gene[1]
  delta <- rep(0.20, n_gene)
  # delta <- rnorm(n_gene, 0.2 , 0.1)  # make the delta positive
  obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
  dat <- simulate.microarry(size, obj)
  
  design <- model.matrix(~dat$trt) 
  fit <- lmFit(dat$data, design)
  fit <- eBayes(fit)		
  pval <- fit$p.value[, 2]
  tval_limma <- fit$t[, 2]
  t1 <- t.test(tval_limma[dat$go_term==1], tval_limma[dat$go_term==0]) # p-value = 8.651e-06
  
  dat$data2 <- standardize.microarray(dat$data, dat$trt)
  fit <- lmFit(dat$data2, design)
  fit <- eBayes(fit)		
  pval <- fit$p.value[, 2]
  tval_limma <- fit$t[, 2]
  t2 <- t.test(tval_limma[dat$go_term==1], tval_limma[dat$go_term==0]) # p-value = 8.651e-06
  
  return(c(t1$p.value, t2$p.value))
  
}

testStandardization("c")




# should we use the covariance or the correlation as an estimate of sample correlation?
## the conclusion here: the correlation will give slightly larger values than the covariance after standardization


case <- "e"
nsim <- 1000
size <- 50           # number of samples to be simulated

# rho <- c(0.7, 0.5, -0.5)   # extreme case
rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation
num_gene <- c(500, 100)
prop <- c(0.30, 0.10)
n_gene <- num_gene[1]
delta <- rep(0.20, n_gene)
# delta <- rnorm(n_gene, 0.2 , 0.1)  # make the delta positive

obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
dat <- simulate.microarry(size, obj)

microarray <- dat$data
trt <- dat$trt
go_term <- dat$go_term

microarray <- standardize.microarray(as.matrix(microarray), trt)
group_mean <- as.matrix(group.mean(microarray, trt))     # calculate correlation matrix
resid_mat <- microarray - group_mean                     # the trt effects are removed from matrix
samp_rho <- cov(t(resid_mat)) 
ids <- which(go_term ==1)
mean(samp_rho[ids, ids]); mean(samp_rho[-ids, ids]); mean(samp_rho[-ids, -ids])


##  see from a real data

CLS <- GSEA.ReadClsFile(file="Gender/Gender.cls")
gender <- GSEA.Gct2Frame2(filename = "Gender/Gender.gct")
geneset <- readGeneSet("Gender/C1.gmt")
go <- geneset$geneSet[[278]][-(1:2)]
trt <- CLS$class.v

gene_id <- which(rownames(gender) %in% go)
go_term <- rep(0, nrow(gender))
go_term[gene_id] <- 1

gender <- standardize.microarray(as.matrix(gender), trt)
group_mean <- as.matrix(group.mean(gender, trt))     # calculate correlation matrix
resid_mat <- gender - group_mean                     # the trt effects are removed from matrix
samp_rho1 <- cor(t(resid_mat))  
samp_rho2 <- cov(t(resid_mat))
ids <- which(go_term ==1)
mean(samp_rho[ids, ids]); mean(samp_rho[-ids, ids]); mean(samp_rho[-ids, -ids])



## look into the inflated type I error problem

source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/R/SimulateLabDataNew.R")                                          # the code
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/R/MOMEnrichmentTest.R") 


library(gap);
Path <-"/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160211/TypeIerrorDiagnostic/"
setwd(Path)

case <- "a0"
dat <- read.table(paste("diagnostic_", case, "_0PCT.txt", sep=""))

par(mfrow=c(2, 2))
hist(dat$testCor); hist(dat$backCor); hist(dat$interCor);  qqunif(dat$p)
x1 <- apply(dat[, 3:5], 2, mean)* c(10000, 160000, 80000) - c(100, 400, 0)
x1/c(10000-100, 160000-400, 80000)


quantile(dat$p, c(0.05, 0.10, 0.50, 0.80))

mean(1/pool)


par(mfrow=c(2, 2))
library(gap)
dat <- read.table("diagnostic_a0_0PCT_Standardization.txt");qqunif(dat$p)
dat <- read.table("diagnostic_a0_10PCT_Standardization.txt");qqunif(dat$p)
dat <- read.table("diagnostic_a0_0PCT_NOStandardization.txt");qqunif(dat$p)
dat <- read.table("diagnostic_a0_0PCT_NOStandardization.txt");qqunif(dat$p)
dat <- read.table("diagnostic_a0_0PCT_NOStandardization.txt");qqunif(dat$p)

plot_pval <- function(case){
  par(mfrow=c(2, 2))
  library(gap)
  dat1 <- read.table(paste("diagnostic_", case, "_0PCT_NOStandardization.txt", sep=""));qqunif(dat1$p, main=paste(case, "0PCT NOstd"))
  dat1 <- read.table(paste("diagnostic_", case, "_0PCT_Standardization.txt", sep=""));qqunif(dat1$p, main=paste(case, "0PCT std"))
  dat1 <- read.table(paste("diagnostic_", case, "_10PCT_NOStandardization.txt", sep=""));qqunif(dat1$p, main=paste(case, "10PCT NOstd"))
  dat1 <- read.table(paste("diagnostic_", case, "_10PCT_Standardization.txt", sep=""));qqunif(dat1$p, main=paste(case, "10PCT std"))
  }

plot_pval("a0")
plot_pval("a")
plot_pval("c")
plot_pval("e")
plot_pval("f")


## type I error simulation plots

FigurePath <-"/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160211/SizePoint05/"
setwd(FigurePath)

showcol <- c(1, 3:8)
textsize = rep(20,20, 6, 20)

case <- "c"
typeIdata <- paste("TypeIerror_", case, "_10PCT_SizePoint1.txt", sep = "")
figure_num <- paste("TypeIerror_", case, "_10PCT", sep = "")
figure_name <- paste("TypeIerror_", case, "_10PCT.eps", sep = "")

create.hist2(read.table(typeIdata)[, showcol], textsize = textsize, figure.num =figure_num, figname = figure_name)

library(gap)

case <- c("a0", "a", "c", "e", "f")
par(mfrow = c(3, 2))
method <- 1
for ( i in case){
  typeIdata <- paste("TypeIerror_", i, "_10PCT_SizePoint05.txt", sep = "")
  print(typeIdata)
  x <- read.table(typeIdata)
  qqunif(x[, method], main = paste("Q-Q plot", i, "10PCT") )
}







case <- "c"
nsim <- 1000
size <- 50           # number of samples to be simulated

# rho <- c(0.7, 0.5, -0.5)   # extreme case
rho <- c(0.1, 0.05, -0.05) # correlation for case a, e, f, used in simulation
num_gene <- c(500, 100)
prop <- c(0.0, 0.0)
n_gene <- num_gene[1]
delta <- rep(0.10, n_gene)
# delta <- rnorm(n_gene, 0.2 , 0.1)  # make the delta positive

nsim <- 100; 
pool <- matrix(NA, nrow= 500, ncol=nsim)
for( i in 1:nsim){
  
  obj <- prepare.simulation(num_gene, prop, delta, case = case, rho)
  dat <- simulate.microarry(size, obj)
  pool[, i] <- standardize.microarray(dat$data, dat$trt)
  
}


par(mfrow=c(1, 1))
hist(pool)
mean(pool)
quantile(pool)
x1 <- colMeans(pool)
