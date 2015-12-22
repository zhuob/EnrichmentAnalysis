# setwd("~/Google Drive/Study/Graduation/Paper3/Data")

# scp /Users/Bin/Google\ Drive/Study/Graduation/Paper3/LinuxServer/pValueGOregress.R zhuob@submit.hpc.cosine.oregonstate.edu:/home/stats/zhuob/Rcode
setwd("/home/stats/zhuob/data/computing")
p.mat <- readRDS("p.NBP.rds")
arab.rm.mean <- readRDS("arab.rm.mean.rds")
ribosome <- read.csv("GOribosome.csv", header=F)
GOribosome <- unique(ribosome[-dim(ribosome)[1], 2])

## get the p values for the first experiment, and estimate COV from the same exp.
id.gene <- which(p.mat$Gene %in% GOribosome)

## genes in GO terms
id.gene2 <- which(row.names(arab.rm.mean) %in% GOribosome)
cor.mat.1 <- cor(t(arab.rm.mean[id.gene2, 1:6]))
p.1 <-  p.mat[id.gene, 2]


## genes not in GO terms
cor.mat.0 <- cor(t(arab.rm.mean[-id.gene2, 1:6]))
p.0 <- p.mat[-id.gene, 2]


library(regress)

model1 <- regress(p.1 ~ 1, ~cor.mat.1)
summary(model1)
model2 <- regress(p.0 ~ 1, ~cor.mat.0)
summary(model2)

## category as a predictor
cor.mat<- cor(t(arab.rm.mean[, 1:6]))
p <- p.mat[, 2]
cat <- ifelse(row.names(arab.rm.mean) %in% GOribosome, 1, 0)
model3 <- regress(p ~ cat, ~cor.mat)
summary(model3)

ls <- list(model1, model2, model3)
saveRDS(ls, "/home/stats/zhuob/data/computing/result_ribosome.rds")

## logistic regression 
p.all <- p.mat[, 2]
log.p <- -log(p.all)
respon <- ifelse(p.mat$Gene %in% GOribosome, 1, 0)

mod3 <- glm(respon ~ log.p, family=binomial(link="logit"))

# fisher's exact test
dichot <- ifelse(p.all < 0.05, 1, 0)
fisher.test(table(data.frame(dichot, respon)))









