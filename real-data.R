setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData/HDresidCorrected/") # the data set


### merge the results from GSEA, CAMERA and OurMethod

## the HD data
d1 <- read.csv("C2.MEQLEA.csv")
d2 <- read.csv("C2.GSEA.10000.csv")
d3 <- read.csv("C2.CAMERA.csv")
d4 <- read.csv("C2.MRSGE.csv")
colnames(d1)[1] <- colnames(d2)[1]<- "set.name"
d5 <- merge(merge(merge(d1, d2, by = "set.name"), d3, by = "set.name"), d4, by = "set.name")
d6 <- data.frame(set.name = d5$set.name, set.size = d5$set_size, 
                 testsetCor = d5$testSetCor, backSetCor = d5$backSetCor, 
                 interCor = d5$interCor,  Camera = d5$PValue, GSEA = d5$NOM.p.val, 
                 MEQLEA = d5$p1, MRSGE = d5$p.MRSGE)
write.csv(d6, "C2.combined.csv", row.names = F)


########  draw some conclusions ------------------------------------------------------





CombinedResults <- read.csv("C2.combined.csv", row.names = NULL)
#CombinedResults <- read.csv("Gender/combined.2.csv")
FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/nar-latex2010/Figures/"
threshold <- 1e-5
wd <- 5
ht <- 5
cor(CombinedResults[, 6:9])[3, ]
cor(log(CombinedResults[, 6:9] + threshold))[3,]

library(ggplot2)

plot1 <- ggplot(data= CombinedResults, aes(-log(CombinedResults$MEQLEA + threshold, 10), -log(CombinedResults$Camera + threshold, 10))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color= "black", size =1) +
  labs(x = "-log10 p values (MEACA)", y = "-log10 p values (CAMERA)") + 
  theme(axis.text=element_text(size=20, face = "bold"), 
        axis.title=element_text(size=20,face="bold")) + 
  scale_y_continuous(limits=c(0, 5)) + 
  scale_x_continuous(limits=c(0, 5))


ggsave(paste(FigurePath,"/MEQLEA_Camera.eps", sep =""), plot1, 
       width = 8, height = 5)

plot2 <- ggplot(data = CombinedResults, aes(-log(MEQLEA + threshold, 10), -log(GSEA + threshold, 10))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1)  +
  labs(x = "-log10 p values (MEACA)", y = "-log10 p values (GSEA)") +
  theme(axis.text=element_text(size=20, face = "bold"), 
        axis.title=element_text(size=20,face="bold")) + 
  scale_y_continuous(limits=c(0, 5)) + 
  scale_x_continuous(limits=c(0, 5))

ggsave(paste(FigurePath,"/MEQLEA_GSEA.eps", sep =""), plot2, 
       width = 8, height = 5)

 plot3 <- ggplot(data = CombinedResults, aes(-log(MEQLEA + threshold, 10), -log(MRSGE + threshold, 10))) + 
   geom_point() + 
   geom_abline(intercept = 0, slope = 1, color = "black", size = 1)  +
   labs(x = "-log10 p values (MEACA)", y = "-log10 p values (MRGSE)") +
   theme(axis.text=element_text(size=20, face = "bold"), 
         axis.title=element_text(size=20,face="bold")) + 
   scale_y_continuous(limits=c(0, 5)) + 
   scale_x_continuous(limits=c(0, 5))
 
 ggsave(paste(FigurePath,"/MEQLEA_MRGSE.eps", sep =""), plot3, 
        width = 8, height = 5)
 
 
 ########### create a table similar to CAMERA, listing top enriched gene sets
 adjust.method <-  "BH"# "bonferroni"
 # https://support.bioconductor.org/p/13804/
 CombinedResults$FDR.OurMethod <- p.adjust(CombinedResults$MEQLEA, method = adjust.method)
 CombinedResults$FDR.Camera<- p.adjust(CombinedResults$Camera, method = adjust.method)
 CombinedResults$FDR.GSEA <- p.adjust(CombinedResults$GSEA, method = adjust.method)
 CombinedResults$FDR.MRSGE <- p.adjust(CombinedResults$MRSGE, method = adjust.method)
 print(sum(CombinedResults$FDR.OurMethod < 0.05))
 print(sum(CombinedResults$FDR.Camera < 0.05))
 print(sum(CombinedResults$FDR.GSEA < 0.05))
 print(sum(CombinedResults$FDR.MRSGE < 0.05))
 # how many sets are overlapped.
 length(which(CombinedResults$FDR.OurMethod < 0.05 & CombinedResults$FDR.GSEA < 0.05))
 length(which(CombinedResults$FDR.OurMethod < 0.05 & CombinedResults$FDR.Camera < 0.05))
 length(which(CombinedResults$FDR.OurMethod < 0.05 & CombinedResults$FDR.MRSGE < 0.05))
 
 
cut_off <- sort(CombinedResults$FDR.OurMethod)[30]
topEnrichedSets <- CombinedResults[CombinedResults$FDR.OurMethod <= cut_off, ]
 # since there are only 33 gene sets, we list all of them
# cut_off <- which(CombinedResults$FDR.OurMethod < 0.05)
# topEnrichedSets <- CombinedResults[cut_off, ]
enrichByGSEA <- which(topEnrichedSets$FDR.GSEA < 0.05)
enrichByMRGSE <- which(topEnrichedSets$FDR.MRSGE <0.05)
topEnrichedSets$indicator <- NA
topEnrichedSets$indicator2 <- NA

topEnrichedSets$indicator[enrichByGSEA] = T
topEnrichedSets$indicator2[enrichByMRGSE] = T

reportEnrichedSets <- topEnrichedSets[order(topEnrichedSets$MEQLEA), c(1:5,8, 10, 14)] # sort the p values
# reportEnrichedSets$set.name <- tolower(reportEnrichedSets$set.name)
reportEnrichedSets <- reportEnrichedSets[order(reportEnrichedSets$FDR.OurMethod, decreasing = F),]

 #reportEnrichedSets$set.name <- gsub("_", "\\1", reportEnrichedSets$set.name)
 library(xtable) 
 # the digits, negative number for scientific notation.
tab <-  xtable(reportEnrichedSets, digits = c(0,0,0, 3, 3, 3, -1, -1, 0))
#align(tab) <- "lp{3in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}{0.5in}"
print(tab, include.rownames=FALSE)
 
x1 <- read.csv("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData/HD/C2.MEQLEA.csv", header=T)
x2 <- x1[x1$p1_fdr<= sort(x1$p1_fdr)[30],]
dim(x2)
intersect(x2$set_name, topEnrichedSets$set.name)

###  the gender data ####---------------------------


# the gender data
d1 <- read.csv("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData/Gender/Gender.meqlea.csv")
d2 <- read.csv("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData/Gender/Gender.CAMERA.csv")
d3 <- read.csv("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData/Gender/Gender.GSEA.csv")
d4 <- read.csv("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData/Gender/Gender.MRSGE.csv")
colnames(d3)[1] <- colnames(d1)[1] <- "set.name"
d5 <- merge(merge(merge(d1, d2, by = "set.name"), d3, by = "set.name"), d4, by = "set.name")
d6 <- data.frame(set.name = d5$set.name, set.size = d5$set_size, 
                 testsetCor = d5$testSetCor, backSetCor = d5$backSetCor, 
                 interCor = d5$interCor,  Camera = d5$PValue, GSEA = d5$NOM.p.val, 
                 MEQLEA = d5$p1, MRSGE = d5$p.MRSGE)
write.csv(d6, "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData/Gender/Gender.combined.csv", row.names = F)




adjust.method  <- "BH"

gender <- read.csv("Results/RealData/Gender/Gender.combined.csv")
gender$FDR.MEQLEA <- p.adjust(gender$MEQLEA, method = adjust.method)
gender$FDR.Camera<- p.adjust(gender$Camera, method = adjust.method)
gender$FDR.GSEA <- p.adjust(gender$GSEA, method = adjust.method)
gender$FDR.MRSGE <- p.adjust(gender$MRSGE, method = adjust.method)

idex <- gender$MEQLEA<0.01 | gender$Camera < 0.01 | gender$GSEA < 0.01 | gender$MRSGE < 0.01
fdr <- 0.05
idex <- gender$FDR.MEQLEA < fdr | gender$FDR.Camera < fdr | gender$FDR.GSEA < fdr | gender$FDR.MRSGE < fdr

print(gender[idex,c(1:5, 10:13)])
# idex <- gender$FDR.OurMethod < 0.05
reportGender <- gender[idex, c(1, 2, 6:10)]
reportGender2 <- data.frame(gene.set = reportGender$set.name, size = reportGender$set.size, 
                            p1 = reportGender$MEQLEA, # FDR.p = reportGender$FDR.OurMethod, 
                            p2 = reportGender$GSEA,#  FDR.GSEA = reportGender$FDR.GSEA, 
                            P3 = reportGender$Camera, #FDR.Camera = reportGender$FDR.Camera
                            p4 = reportGender$MRSGE)
reportGender2 <- reportGender2[order(reportGender2$p1), ]
genderTable <- xtable(reportGender2, digits = c(0, 0, 0, -1, -1, -1, -1), label = "table:gender")
print(genderTable, include.rownames = F)


# ## overlapped gene sets with Diamanti 2013 ###
# Diamanti <- read.csv("HuntingtonDisease/Diamanti2013.csv", header= T)
# Diamanti$set.name <- paste(Diamanti$Path, "_", Diamanti$Gene_set_name, sep ="")
# Diamanti_Overlap <- merge(Diamanti, CombinedResults, by = "set.name")
# print(Diamanti_Overlap$set.name[Diamanti_Overlap$FDR.OurMethod < 0.05])








