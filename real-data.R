setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/DataSet/") # the data set


### merge the results from GSEA, CAMERA and OurMethod

## the HD data
d1 <- read.csv("HuntingtonDisease/C2.CAMERA.csv")
d2 <- read.csv("HuntingtonDisease/C2.GSEA.csv")
d3 <- read.csv("HuntingtonDisease/C2.OurMethod.csv")
colnames(d2)[1] <- "set.name"
d4 <- merge(d3, merge(d1, d2, by = "set.name"), by = "set.name")
d5 <- data.frame(set.name = d4$set.name, set.size = d4$set.size, testsetCor = d4$testSetCor,
                 backSetCor = d4$backSetCor, interCor = d4$interCor,  Camera = d4$PValue, GSEA = d4$NOM.p.val, p = d4$p)
write.csv(d5, "HuntingtonDisease/C2.combinedNew.csv", row.names = F)

# the gender data
d1 <- read.csv("Gender/Gender.CAMERA.csv")
d2 <- read.csv("Gender/Gender.GSEA.csv")
d3 <- read.csv("Gender/Gender.OurMethod.csv")
colnames(d2)[1] <- "set.name"
d4 <- merge(d3, merge(d1, d2, by = "set.name"), by = "set.name")
d5 <- data.frame(set.name = d4$set.name, set.size = d4$set.size, testsetCor = d4$testSetCor,
                 backSetCor = d4$backSetCor, interCor = d4$interCor,  Camera = d4$PValue, GSEA = d4$NOM.p.val, p = d4$p)
write.csv(d5, "Gender/Gender.combinedNew.csv", row.names = F)



########  draw some conclusions ------------------------------------------------------



CombinedResults <- read.csv("HuntingtonDisease/C2.Combined.All.csv", row.names = NULL)
#CombinedResults <- read.csv("Gender/combined.2.csv")
FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"
threshold <- 1e-6
wd <- 5
ht <- 5

library(ggplot2)


plot1 <- ggplot(data = CombinedResults, aes((p + threshold), (Camera + threshold))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color= "black", size =1) +
  labs(x = "p values (MEQLEA)", y = "p values (CAMERA)") + 
  theme(axis.text=element_text(size=20, face = "bold"), 
        axis.title=element_text(size=20,face="bold")) +
  scale_y_continuous(limits=c(0, 1))

ggsave(paste(FigurePath,"/P_Camera.eps", sep =""), plot1, 
       width = 8, height = 5)

plot2 <- ggplot(data = CombinedResults, aes((p + threshold), (GSEA + threshold))) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, color = "black", size = 1)  +
  labs(x = "p values (MEQLEA)", y = "p values (GSEA)") +
  theme(axis.text=element_text(size=20, face = "bold"), 
        axis.title=element_text(size=20,face="bold"))

ggsave(paste(FigurePath,"/P_GSEA.eps", sep =""), plot2, 
       width = 8, height = 5)

 plot3 <- ggplot(data = CombinedResults, aes((Camera + threshold), (GSEA + threshold))) + 
   geom_point() + 
   geom_abline(intercept = 0, slope = 1, color = "black", size = 1)  +
   labs(x = "p values (Camera)", y = "p values (GSEA)") +
   theme(axis.text=element_text(size=20, face = "bold"), 
         axis.title=element_text(size=20,face="bold"))
 
 ggsave(paste(FigurePath,"/P_GSEA_CAMERA.eps", sep =""), plot3, 
        width = 8, height = 5)
 
 
 ########### create a table similar to CAMERA, listing top enriched gene sets
 CombinedResults$FDR.OurMethod <- p.adjust(CombinedResults$p, method = "BH")
 CombinedResults$FDR.Camera<- p.adjust(CombinedResults$Camera, method = "BH")
 CombinedResults$FDR.GSEA <- p.adjust(CombinedResults$GSEA, method = "BH")
 CombinedResults$FDR.MRGSE <- p.adjust(CombinedResults$p.MRGSE, method = "BH")
 print(sum(CombinedResults$FDR.OurMethod < 0.05))
 print(sum(CombinedResults$FDR.Camera < 0.05))
 print(sum(CombinedResults$FDR.GSEA < 0.05))
 print(sum(CombinedResults$FDR.MRGSE < 0.05))
 # how many sets are overlapped.
 length(which(CombinedResults$FDR.OurMethod < 0.05 & CombinedResults$FDR.GSEA < 0.05))
 length(which(CombinedResults$FDR.OurMethod < 0.05 & CombinedResults$FDR.Camera < 0.05))
 length(which(CombinedResults$FDR.OurMethod < 0.05 & CombinedResults$FDR.MRGSE < 0.05))
 
 
cut_off <- sort(CombinedResults$FDR.OurMethod)[30]
topEnrichedSets <- CombinedResults[CombinedResults$FDR.OurMethod <= cut_off, ]
enrichByGSEA <- which(topEnrichedSets$FDR.GSEA < 0.05)
enrichByMRGSE <- which(topEnrichedSets$FDR.MRGSE <0.05)
topEnrichedSets$indicator <- NA
topEnrichedSets$indicator2 <- NA

topEnrichedSets$indicator[enrichByGSEA] = T
topEnrichedSets$indicator2[enrichByMRGSE] = T

reportEnrichedSets <- topEnrichedSets[order(topEnrichedSets$p), c(1:5,8, 11,14)] # sort the p values
# reportEnrichedSets$set.name <- tolower(reportEnrichedSets$set.name)
reportEnrichedSets <- reportEnrichedSets[order(reportEnrichedSets$FDR.OurMethod, decreasing = F),]

 #reportEnrichedSets$set.name <- gsub("_", "\\1", reportEnrichedSets$set.name)
 library(xtable) 
 # the digits, negative number for scientific notation.
tab <-  xtable(reportEnrichedSets, digits = c(0,0,0, 3, 3, 3, -1, -1, 0))
align(tab) <- "lp{3in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}p{0.5in}"
print(tab, include.rownames=FALSE)
 



###  the gender data ####---------------------------


gender <- read.csv("Gender/Gender.Combined.All.csv")
gender$FDR.OurMethod <- p.adjust(gender$p, method = "BH")
gender$FDR.Camera<- p.adjust(gender$Camera, method = "BH")
gender$FDR.GSEA <- p.adjust(gender$GSEA, method = "BH")
gender$FDR.MRGSE <- p.adjust(gender$p.MRGSE, method = "BH")

idex <- gender$p<0.01 | gender$Camera < 0.01 | gender$GSEA < 0.01 | gender$p.MRGSE < 0.01
# idex <- gender$FDR.OurMethod < 0.05
reportGender <- gender[idex, c(1, 2, 6:13)]
reportGender2 <- data.frame(gene.set = reportGender$set.name, size = reportGender$set.size, 
                            p1 = reportGender$p, FDR.p = reportGender$FDR.OurMethod, 
                            p2 = reportGender$GSEA, FDR.GSEA = reportGender$FDR.GSEA, 
                            P3 = reportGender$Camera, FDR.Camera = reportGender$FDR.Camera)
reportGender2 <- reportGender2[order(reportGender2$p1), ]
genderTable <- xtable(reportGender2, digits = c(0, 0, 0,  3, 3, 3, 3, 3, 3), label = "table:gender")
print(genderTable, include.rownames = F)


## overlapped gene sets with Diamanti 2013 ###
Diamanti <- read.csv("HuntingtonDisease/Diamanti2013.csv", header= T)
Diamanti$set.name <- paste(Diamanti$Path, "_", Diamanti$Gene_set_name, sep ="")
Diamanti_Overlap <- merge(Diamanti, CombinedResults, by = "set.name")
print(Diamanti_Overlap$set.name[Diamanti_Overlap$FDR.OurMethod < 0.05])








