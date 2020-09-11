rm(list = ls())
# path <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/Simulation2017/Normal2"
path <- "simulation-results/gene-500/n25vs25/DEmean2sd1"
source("R-code-paper/plot-functions.R")

methods_explored <- c("MEACA", "MRGSE", "SigPathway", "CAMERA_ModT", 
                      "CAMERA_Rank", "GSEA", "QuSAGE", "ORA")

# showcol <- c(1, 3:9)
showcol_body <- c(1, 2, 4, 6, 8)
showcol_appendix <- c(1, 3, 5, 7)
showcol <- showcol_appendix
selected_m1 <- methods_explored[showcol]

color1 <- "#0099FF"
textsize <- c(15, 15,  10, 10)

fig1 <- plot_fig1(path = path, textsize = textsize, selected_method = methods_explored[showcol_body],
                  color1 = color1, case = letters[1:5])
cowplot::ggsave(paste(path,"/fig1-body", ".eps", sep =""), fig1, width = 18, height = 20)
fig2 <- plot_fig1(path = path, textsize = textsize, selected_method = methods_explored[showcol_appendix],
                  color1 = color1, case = letters[1:5])
cowplot::ggsave(paste(path,"/fig1-supp", ".eps", sep =""), fig2, width = 18, height = 20)

fig1 <- plot_fig1(path = path, textsize = textsize, selected_method = methods_explored,
                  color1 = color1, case = letters[1:3])
cowplot::ggsave(paste(path,"/corr-a-c", ".pdf", sep =""), fig1, width = 18, height = 20)
fig1 <- plot_fig1(path = path, textsize = textsize, selected_method = methods_explored,
                  color1 = color1, case = letters[4:6])
cowplot::ggsave(paste(path,"/corr-d-e", ".pdf", sep =""), fig1, width = 18, height = 20)


fig2 <- plot_fig2(path, showcol = 1)
ggsave(paste(path,"/fig2", ".eps", sep =""), fig2, width = 18, height = 20)


table1 <- createEmpiricalPowerTable(path, showcol = c(1, 4:9), setup = "A1", case = "a")
table2 <- createEmpiricalPowerTable(path, showcol = c(1, 4:9), setup = "A2", case = "a")


#path <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/Simulation2017/Normal2/"
# below is path to the results from the paper
# path <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/old/Power20161223/"
showcol <- c(1, 4, 5, 6, 8, 9)
ta <- createEmpiricalPowerTable(path, setup = "A2", case = "a", showcol = showcol)
tb <- createEmpiricalPowerTable(path, setup = "A2", case = "b", showcol = showcol)
tc <- createEmpiricalPowerTable(path, setup = "A2", case = "c", showcol = showcol)
td <- createEmpiricalPowerTable(path, setup = "A2", case = "d", showcol = showcol)
te <- createEmpiricalPowerTable(path, setup = "A2", case = "e", showcol = showcol)

s <- 1
table3 <- cbind(ta[, s], tb[, s], tc[, s], td[, s], te[, s])
table3

