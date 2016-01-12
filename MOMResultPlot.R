
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/summarize.results.R")


setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationMOMEnrichment/")
create.hist2(read.table("DE_a_6.txt"), figure.num = "DE_a_6")

create.hist2(read.table("DE_a_50_sample.txt"), figure.num = "DE_a_50_sample")
create.hist2(read.table("DE_a_50_random.txt"), figure.num = "DE_a_50_random")
create.hist2(read.table("DE_a_50_true.txt"), figure.num = "DE_a_50_true" )

create.hist2(read.table("DE_a_50.txt"), figure.num = "DE_a_50")
create.hist2(read.table("DE_b_50.txt"), figure.num = "DE_b_50")
create.hist2(read.table("DE_c_50.txt"), figure.num = "DE_c_50")
create.hist2(read.table("DE_d_50.txt"), figure.num = "DE_d_50")
create.hist2(read.table("DE_e_50.txt"), figure.num = "DE_e_50")
create.hist2(read.table("DE_f_50.txt"), figure.num = "DE_f_50")


create.hist2(read.table("NO_DE_a_50.txt"), figure.num = "NO_DE_a_50")
create.hist2(read.table("NO_DE_c_50.txt"), figure.num = "NO_DE_c_50")
create.hist2(read.table("NO_DE_e_50.txt"), figure.num = "NO_DE_e_50")
create.hist2(read.table("NO_DE_f_50.txt"), figure.num = "NO_DE_f_50")
 

# power: go 20% DE, no_go, 0% DE
create.hist2(read.table("DE_a_50_power.txt"), figure.num = "DE_a_50_power")
create.hist2(read.table("DE_b_50_power.txt"), figure.num = "DE_b_50_power")
create.hist2(read.table("DE_c_50_power.txt"), figure.num = "DE_c_50_power")
create.hist2(read.table("DE_d_50_power.txt"), figure.num = "DE_d_50_power")
create.hist2(read.table("DE_e_50_power.txt"), figure.num = "DE_e_50_power")
create.hist2(read.table("DE_f_50_power.txt"), figure.num = "DE_f_50_power")





setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationRealCor20151207/")


##  standardization 


##  standardization
create.hist2(read.table("DE_a0_50_New.txt"), figure.num = "DE_a0_50")
create.hist2(read.table("DE_a_50_new.txt"), figure.num = "DE_a_50")
create.hist2(read.table("DE_b_50_new.txt"), figure.num = "DE_b_50")
create.hist2(read.table("DE_c_50_new.txt"), figure.num = "DE_c_50")
create.hist2(read.table("DE_d_50_new.txt"), figure.num = "DE_d_50")
create.hist2(read.table("DE_e_50_new.txt"), figure.num = "DE_e_50")
create.hist2(read.table("DE_f_50_new.txt"), figure.num = "DE_f_50")
create.hist2(read.table("DE_g_50_New.txt"), figure.num = "DE_g_50")
create.hist2(read.table("DE_h_50_new.txt"), figure.num = "DE_h_50")


create.hist2(read.table("NO_DE_a_50_New.txt"), figure.num = "NO_DE_a_50")
create.hist2(read.table("NO_DE_c_50_New.txt"), figure.num = "NO_DE_c_50")
create.hist2(read.table("NO_DE_e_50_New.txt"), figure.num = "NO_DE_e_50")
create.hist2(read.table("NO_DE_f_50_New.txt"), figure.num = "NO_DE_f_50")
create.hist2(read.table("NO_DE_g_50_New.txt"), figure.num = "NO_DE_g_50")




create.hist2(read.table("Power_a_50_New.txt"), figure.num = "Power_DE_a_50")
create.hist2(read.table("Power_c_50_New.txt"), figure.num = "Power_DE_c_50")
create.hist2(read.table("Power_e_50_New.txt"), figure.num = "Power_DE_e_50")
create.hist2(read.table("Power_f_50_New.txt"), figure.num = "Power_DE_f_50")
create.hist2(read.table("Power_g_50_New.txt"), figure.num = "Power_DE_g_50")


n <- 1000; nu <- 4;tau <- 0.25
x1 <- 1/rchisq(n, nu)                                   # inverse chi square distribution
s2 <- tau^2*nu*x1                                       # scaled inverse chi-square distribution
stdevs <- sqrt(s2)        



# add GESA 
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationGSEA20151223/")
FigurePath <-"/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"
  
showcol <- c(1, 3:7)
textsize = rep(20,20, 8, 20)

data <- read.table("DE_a0_50.txt")[, showcol]

create.hist2(read.table("DE_a0_50.txt")[, showcol], textsize = textsize, figure.num = "DE_a0_50")
create.hist2(read.table("DE_a_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_a_50")
create.hist2(read.table("DE_b_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_b_50")
create.hist2(read.table("DE_c_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_c_50")
create.hist2(read.table("DE_d_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_d_50")
create.hist2(read.table("DE_e_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_e_50")
create.hist2(read.table("DE_f_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_f_50")
create.hist2(read.table("DE_g_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_g_50")
create.hist2(read.table("DE_h_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_h_50")

pvec <- c(0.01, 0.05, 0.10, 0.20)
A0_NO_DE <- TypeIerror(pvec, read.table("NO_DE_a0_50.txt")[, showcol])
A_NO_DE <- TypeIerror(pvec, read.table("NO_DE_a_50.txt")[, showcol])
C_NO_DE <- TypeIerror(pvec, read.table("NO_DE_c_50.txt")[, showcol])
E_NO_DE <- TypeIerror(pvec, read.table("NO_DE_e_50.txt")[, showcol])
F_NO_DE <- TypeIerror(pvec, read.table("NO_DE_f_50.txt")[, showcol])
G_NO_DE <- TypeIerror(pvec, read.table("NO_DE_g_50.txt")[, showcol])

A0_NO_DE
A_NO_DE
C_NO_DE
E_NO_DE
F_NO_DE
G_NO_DE



create.hist2(read.table("NO_DE_a0_50.txt")[, showcol], textsize = textsize, figure.num ="a0", figname = "NODEA0.eps")
create.hist2(read.table("NO_DE_a_50.txt")[, showcol], textsize = textsize, figure.num = "a", figname = "NODEA.eps")
create.hist2(read.table("NO_DE_b_50.txt")[, showcol], textsize = textsize, figure.num = "NO_DE_b_50")
create.hist2(read.table("NO_DE_c_50.txt")[, showcol], textsize = textsize, figure.num = "c", figname = "NODEC.eps")
create.hist2(read.table("NO_DE_d_50.txt")[, showcol], textsize = textsize, figure.num = "NO_DE_d_50")
create.hist2(read.table("NO_DE_e_50.txt")[, showcol], textsize = textsize, figure.num = "e", figname = "NODEE.eps")
create.hist2(read.table("NO_DE_f_50.txt")[, showcol], textsize = textsize, figure.num = "f", figname = "NODEF.eps")
create.hist2(read.table("NO_DE_g_50.txt")[, showcol], textsize = textsize, figure.num = "g", figname = "NODEG.eps")
create.hist2(read.table("NO_DE_h_50.txt")[, showcol], textsize = textsize, figure.num = "NO_DE_h_50")



pvec <- c(0.01, 0.05, 0.10, 0.20)
A0_DE <- TypeIerror(pvec, read.table("DE_a0_50.txt")[, showcol])
A_DE <- TypeIerror(pvec, read.table("DE_a_50.txt")[, showcol])
C_DE <- TypeIerror(pvec, read.table("DE_c_50.txt")[, showcol])
E_DE <- TypeIerror(pvec, read.table("DE_e_50.txt")[, showcol])
F_DE <- TypeIerror(pvec, read.table("DE_f_50.txt")[, showcol])
G_DE <- TypeIerror(pvec, read.table("DE_g_50.txt")[, showcol])

l1 <- rbind(A0_DE, C_DE, F_DE); l2 <- rbind(A_DE, E_DE, G_DE)
method <- rownames(l1); ab <- data.frame(method, l1, method, l2); ab$method.1 <- ""

library(xtable)
print(xtable(ab, digits = 3), include.rownames = F)



P_a <- TypeIerror(pvec, read.table("Power_a_50.txt")[, showcol])
P_g <- TypeIerror(pvec, read.table("Power_g_50.txt")[, showcol])
P_a
P_g
