# run our method with gender data from GSEA paper


source("/home/stats/zhuob/Rcode/Enrichment/MOMEnrichmentTest.R")
source("/home/stats/zhuob/Rcode/Enrichment/SimulateLabDataNew.R")               # it contains group.mean functi$
source("/home/stats/zhuob/Rcode/Enrichment/GSEA.1.0.R")

setwd("/home/stats/zhuob/data/computing/DataSet/")

CLS <- GSEA.ReadClsFile(file="Gender.cls")
dataset <- GSEA.Gct2Frame2(filename = "Gender.gct")
geneset <- readGeneSet("C1.gmt")


## for OurMethod
results <- MOM_test_multiple(dataset, CLS$class.v, geneset, minSetSize=2)
write.csv(results, "Gender.OurMethod.csv", row.names =F)

## for CAMERA
Camera_results <- Camera_multiple(dataset, CLS$class.v, geneset, use.rank=F)
Camera_results$set.name <- rownames(Camera_results)
write.csv(Camera_results, "Gender.CAMERA.csv", row.names=F)


## for MRSGE
MRSGE_results <-  MRSGE_multiple(dataset, CLS$class, geneset, use.rank = T)
write.csv(MRSGE_results, "Gender.MRSGE.csv", row.names=F)
