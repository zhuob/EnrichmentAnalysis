setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/RealData2017/") # the data set


### merge the results from GSEA, CAMERA and OurMethod

## the HD data
d1 <- read.csv("C2.MEACA.csv")
d2 <- read.csv("C2.GSEA.9999.csv")
d3 <- read.csv("C2.CAMERA.csv")
d4 <- read.csv("C2.MRGSE.csv")
colnames(d1)[1] <- colnames(d2)[1]<- "set.name"
d5 <- merge(merge(merge(d1, d2, by = "set.name"), d3, by = "set.name"), d4, by = "set.name")
d6 <- data.frame(set.name = d5$set.name, set.size = d5$set_size, 
                 testsetCor = d5$testSetCor, backSetCor = d5$backSetCor, 
                 interCor = d5$interCor,  Camera = d5$PValue, GSEA = d5$NOM.p.val, 
                 MEACA = d5$p1, MRGSE = d5$p.MRGSE)
write.csv(d6, "HD.combined.csv", row.names = F)


########  draw some conclusions ------------------------------------------------------



CombinedResults <- read.csv("HD.combined.csv", row.names = NULL)
# adjust the GSEA p values by (b + 1)/(K + 1)
CombinedResults$GSEA <- (CombinedResults$GSEA*9999 +1 )/10000

#CombinedResults <- read.csv("Gender/combined.2.csv")
FigurePath <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/nar-latex2010/Figures/"
threshold <- 1e-5
wd <- 5
ht <- 5
limit = c(0, 5)
cor(CombinedResults[, 6:9])[3, ]
cor(log(CombinedResults[, 6:9] + threshold))[3,]

library(ggplot2)
library(cowplot)

plot_fig3 <- function(data, threshold, textsize = c(20, 20), limit = c(0, 5)){
  
  plot1 <- ggplot(data= data, aes(-log(MEACA + threshold, 10), -log(Camera + threshold, 10))) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color= "black", size =1) +
    labs(x = "", y = "-log10 p values (CAMERA)") + 
    theme(axis.text=element_text(size=textsize[1], face = "bold"), 
          axis.title=element_text(size=textsize[2],face="bold")) + 
    scale_y_continuous(limits=ylimit) + 
    scale_x_continuous(limits=limit)
  
  plot2 <- ggplot(data = data, aes(-log(MEACA + threshold, 10), -log(GSEA + threshold, 10))) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color = "black", size = 1)  +
    labs(x = "", y = "-log10 p values (GSEA)") +
    theme(axis.text=element_text(size=textsize[1], face = "bold"), 
          axis.title=element_text(size=textsize[2],face="bold")) + 
    scale_y_continuous(limits=limit) + 
    scale_x_continuous(limits=limit)
  
  
  plot3 <- ggplot(data = data, aes(-log(MEACA + threshold, 10), -log(MRGSE + threshold, 10))) + 
    geom_point() + 
    geom_abline(intercept = 0, slope = 1, color = "black", size = 1)  +
    labs(x = "-log10 p values (MEACA)", y = "-log10 p values (MRGSE)") +
    theme(axis.text=element_text(size=textsize[1], face = "bold"), 
          axis.title=element_text(size=textsize[2],face="bold")) + 
    scale_y_continuous(limits=limit) + 
    scale_x_continuous(limits=limit)
  

  fig <- plot_grid(plot1, plot2, plot3, labels = c("a", "b", "c"), ncol=1, nrow=3)
  return(fig)
}

fig3 <- plot_fig3(data = CombinedResults, threshold = threshold, textsize = c(12, 15))


 
 ggsave("fig3.eps", fig3,  width = 8, height = 15)
 
 
 ########### create a table similar to CAMERA, listing top enriched gene sets
 adjust.method <-  "BH"# "bonferroni"
 # https://support.bioconductor.org/p/13804/
 CombinedResults$FDR.OurMethod <- p.adjust(CombinedResults$MEACA, method = adjust.method)
 CombinedResults$FDR.Camera<- p.adjust(CombinedResults$Camera, method = adjust.method)
 CombinedResults$FDR.GSEA <- p.adjust(CombinedResults$GSEA, method = adjust.method)
 CombinedResults$FDR.MRSGE <- p.adjust(CombinedResults$MRGSE, method = adjust.method)
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

reportEnrichedSets <- topEnrichedSets[order(topEnrichedSets$MEACA), c(1:5,8, 10, 14)] # sort the p values
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
d1 <- read.csv("Gender.MEACA.csv")
d2 <- read.csv("Gender.CAMERA.csv")
d3 <- read.csv("Gender.GSEA.9999.csv")
d4 <- read.csv("Gender.MRGSE.csv")
colnames(d3)[1] <- colnames(d1)[1] <- "set.name"
d5 <- merge(merge(merge(d1, d2, by = "set.name"), d3, by = "set.name"), d4, by = "set.name")
d6 <- data.frame(set.name = d5$set.name, set.size = d5$set_size, 
                 testsetCor = d5$testSetCor, backSetCor = d5$backSetCor, 
                 interCor = d5$interCor,  Camera = d5$PValue, GSEA = d5$NOM.p.val, 
                 MEACA = d5$p1, MRGSE = d5$p.MRGSE)
write.csv(d6, "Gender.combined.csv", row.names = F)




adjust.method  <- "BH"

gender <- read.csv("Gender.combined.csv")
# adjust the GSEA p values by (b + 1)/(K + 1)
gender$GSEA <- (gender$GSEA*9999 +1 )/10000
gender$FDR.MEACA <- p.adjust(gender$MEACA, method = adjust.method)
gender$FDR.Camera<- p.adjust(gender$Camera, method = adjust.method)
gender$FDR.GSEA <- p.adjust(gender$GSEA, method = adjust.method)
gender$FDR.MRGSE <- p.adjust(gender$MRGSE, method = adjust.method)

idex <- gender$MEACA <0.01 | gender$Camera < 0.01 | gender$GSEA < 0.01 | gender$MRGSE < 0.01
fdr <- 0.05
idex <- gender$FDR.MEACA < fdr | gender$FDR.Camera < fdr | gender$FDR.GSEA < fdr | gender$FDR.MRGSE < fdr

print(gender[idex,c(1:5, 10:13)])
# idex <- gender$FDR.OurMethod < 0.05
reportGender <- gender[idex, c(1, 2, 6:10)]
reportGender2 <- data.frame(gene.set = reportGender$set.name, size = reportGender$set.size, 
                            p1 = reportGender$MEACA, # FDR.p = reportGender$FDR.OurMethod, 
                            p2 = reportGender$GSEA,#  FDR.GSEA = reportGender$FDR.GSEA, 
                            P3 = reportGender$Camera, #FDR.Camera = reportGender$FDR.Camera
                            p4 = reportGender$MRGSE)
reportGender2 <- reportGender2[order(reportGender2$p1), ]
genderTable <- xtable(reportGender2, digits = c(0, 0, 0, -1, -1, -1, -1), label = "table:gender")
print(genderTable, include.rownames = F)


# ## overlapped gene sets with Diamanti 2013 ###
# Diamanti <- read.csv("HuntingtonDisease/Diamanti2013.csv", header= T)
# Diamanti$set.name <- paste(Diamanti$Path, "_", Diamanti$Gene_set_name, sep ="")
# Diamanti_Overlap <- merge(Diamanti, CombinedResults, by = "set.name")
# print(Diamanti_Overlap$set.name[Diamanti_Overlap$FDR.OurMethod < 0.05])








