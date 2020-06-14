library(meaca)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)

path <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/Simulation2017/Normal2"


fig1 <- plot_fig1(path=path)
ggsave(paste(path,"/fig1", ".eps", sep =""), fig1, width = 18, height = 20)

fig2 <- plot_fig2(path, showcol = 1)
ggsave(paste(path,"/fig2", ".eps", sep =""), fig2, width = 18, height = 20)


table1 <- createEmpiricalPowerTable(path, showcol = c(1, 4:9), setup = "A1", case = "a")
table2 <- createEmpiricalPowerTable(path, showcol = c(1, 4:9), setup = "A2", case = "a")


path <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/Power20161223"
path <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/oldFiles/Share/Simulation/Simulation20160318/Normal1"
showcol <- c(1, 4, 5, 6, 8, 9)
ta <- createEmpiricalPowerTable(path, setup = "A2", case = "a", showcol = showcol)
tb <- createEmpiricalPowerTable(path, setup = "A2", case = "b", showcol = showcol)
tc <- createEmpiricalPowerTable(path, setup = "A2", case = "c", showcol = showcol)
td <- createEmpiricalPowerTable(path, setup = "A2", case = "d", showcol = showcol)
te <- createEmpiricalPowerTable(path, setup = "A2", case = "e", showcol = showcol)

s <- 1
table3 <- cbind(ta[, s], tb[, s], tc[, s], td[, s], te[, s])
table3

