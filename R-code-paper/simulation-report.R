rm(list = ls())
# path <- "/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Results/Simulation2017/Normal2"
path <- "simulation-results/gene-20000/n25vs25/DEmean2sd1"
source("R-code-paper/plot-functions.R")

color1= "#0099FF"
color2 = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
           "#F0E442", "#0072B2", "#D55E00", "#CC33CC")
line_values = c("twodash", "dotted", "dotdash", "dashed", 
               "F1", "4C88C488", "12345678", "longdash")
textsize <- c(15, 15,  10, 10)
showcol <- c(1, 3:9)

fig1 <- plot_fig1(path = path, textsize = textsize, showcol = showcol, 
                  color1 = color1, color2 = color2,  
                  case = letters[1:3],
                  line_vaues = line_vaues)
cowplot::ggsave(paste(path,"/corr-a-c", ".pdf", sep =""), fig1, width = 18, height = 20)

fig1 <- plot_fig1(path = path, textsize = textsize, showcol = showcol, 
                  color1 = color1, color2 = color2,  
                  case = c("d", "e", "f"), 
                  line_vaues = line_vaues)
cowplot::ggsave(paste(path,"/corr-d-e", ".pdf", sep =""), fig1, width = 18, height = 20)


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

