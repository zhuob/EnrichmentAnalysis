# test "meaca" on real data

library(MEQLEA)
library(meaca)
library(limma)

setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Share/DataSet/Gender/")

# these two functions come from MEQLEA, only used to read data of this particular type
CLS <- GSEA.ReadClsFile(file="Gender.cls")  # treatment lables
dataset <- GSEA.Gct2Frame2(filename = "Gender.gct") # expression data
geneset <- read_gene_set("C1.gmt")


results <- meaca_multiple(dataset, CLS$class.v, geneset, minSetSize=2)
names(results)
head(results)  # the results are exactly the same as we did before.
