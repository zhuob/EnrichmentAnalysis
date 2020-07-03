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
