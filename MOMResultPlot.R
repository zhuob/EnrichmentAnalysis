
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
textsize = rep(20,20, 6, 20)

data <- read.table("DE_a0_50.txt")[, showcol]

create.hist2(read.table("DE_a0_50.txt")[, showcol], textsize = textsize, figure.num = "DE_a0_50", figname = "DE_a0_50")
create.hist2(read.table("DE_a_50.txt")[, showcol],  textsize = textsize, figure.num = "DE_a_50", figname = "DE_a_50")
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



create.hist2(read.table("NO_DE_a0_50.txt")[, showcol], textsize = textsize, figure.num ="a0", figname = "DEA0.eps")
create.hist2(read.table("NO_DE_a_50.txt")[, showcol], textsize = textsize, figure.num = "a", figname = "DEA.eps")
create.hist2(read.table("NO_DE_b_50.txt")[, showcol], textsize = textsize, figure.num = "NO_DE_b_50")
create.hist2(read.table("NO_DE_c_50.txt")[, showcol], textsize = textsize, figure.num = "c", figname = "DEC.eps")
create.hist2(read.table("NO_DE_d_50.txt")[, showcol], textsize = textsize, figure.num = "NO_DE_d_50")
create.hist2(read.table("NO_DE_e_50.txt")[, showcol], textsize = textsize, figure.num = "e", figname = "DEE.eps")
create.hist2(read.table("NO_DE_f_50.txt")[, showcol], textsize = textsize, figure.num = "f", figname = "DEF.eps")
create.hist2(read.table("NO_DE_g_50.txt")[, showcol], textsize = textsize, figure.num = "g", figname = "DEG.eps")
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



#################################################################################
##  add qusage

## simulation setup  for this data 
nsim <- 1000

size <- 50           # number of samples to be simulated
rho <- c(0.1, 0.05, -0.05)   # correlation for case a, e, f
num_gene <- c(500, 100)
prop <- c(0.2, 0.2)
n_gene <- num_gene[1]
delta <- rnorm(n_gene, 0.5 , 1)  # the DE is very small 



setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationQUsage20160125/")
FigurePath <-"/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/Manuscript/Figures/"

showcol <- c(1, 4:8)
textsize = rep(20,20, 8, 20)

data <- read.table("NO_DE_a0_50.txt")[, showcol]
create.hist2(read.table("NO_DE_a0_50.txt")[, showcol], textsize = textsize, figure.num ="a0", figname = "NODEA0.eps")
create.hist2(read.table("NO_DE_a_50.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "NODEA.eps")
create.hist2(read.table("NO_DE_c_50.txt")[, showcol], textsize = textsize, figure.num ="c", figname = "NODEC.eps")
create.hist2(read.table("NO_DE_e_50.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "NODEE.eps")
create.hist2(read.table("NO_DE_f_50.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "NODEF.eps")
create.hist2(read.table("NO_DE_g_50.txt")[, showcol], textsize = textsize, figure.num ="g", figname = "NODEG.eps")

create.hist2(read.table("DE_a0_50.txt")[, showcol], textsize = textsize, figure.num ="a0", figname = "DEA0.eps")
create.hist2(read.table("DE_a_50.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DEA.eps")
create.hist2(read.table("DE_c_50.txt")[, showcol], textsize = textsize, figure.num ="c", figname = "DEC.eps")
create.hist2(read.table("DE_e_50.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DEE.eps")
create.hist2(read.table("DE_f_50.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "DEF.eps")
create.hist2(read.table("DE_g_50.txt")[, showcol], textsize = textsize, figure.num ="g", figname = "DEG.eps")


create.hist2(read.table("DE_a_50_large.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DEA_LargeRho.eps")
create.hist2(read.table("DE_e_50_large.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DEe_LargeRho.eps")
create.hist2(read.table("DE_f_50_large.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "DEf_LargeRho.eps")


create.hist2(read.table("DE_a_50_12.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DEA_sample12.eps")
create.hist2(read.table("DE_e_50_12.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DEe_sample12.eps")


create.hist2(read.table("DE_a_50_DE10PCT.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DEa_10PCT.eps")
create.hist2(read.table("DE_e_50_DE10PCT.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DEe_10PCT.eps")
create.hist2(read.table("DE_f_50_DE10PCT.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "DEf_10PCT.eps")


##########################################   Jan 28th ##########################
source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/summarize.results.R")
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160128/")
FigurePath <-"/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160128/"

showcol <- c(1, 4:8)
textsize = rep(20,20, 8, 20)
create.hist2(read.table("DE_a_50_DE10PCT.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DEa_10PCT.eps")
create.hist2(read.table("DE_e_50_DE10PCT.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DEe_10PCT.eps")
create.hist2(read.table("DE_f_50_DE10PCT.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "DEf_10PCT.eps")

create.hist2(read.table("DE_a_50_large.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DEA_LargeRho.eps")
create.hist2(read.table("DE_e_50_large.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DEe_LargeRho.eps")
create.hist2(read.table("DE_f_50_large.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "DEf_LargeRho.eps")

create.hist2(read.table("Power_a_50_a.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "Power_a_a.eps")
create.hist2(read.table("Power_e_50_a.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "Power_e_a.eps")
create.hist2(read.table("Power_f_50_a.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "Power_f_a.eps")

create.hist2(read.table("Power_a_50_b.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "Power_a_b.eps")
create.hist2(read.table("Power_e_50_b.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "Power_e_b.eps")
create.hist2(read.table("Power_f_50_b.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "Power_f_b.eps")

create.hist2(read.table("DE_a0_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="a0", figname = "DE_a0_50_10PCT.eps")
create.hist2(read.table("DE_a_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DE_a_50_10PCT.eps")
create.hist2(read.table("DE_c_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="c", figname = "DE_c_50_10PCT.eps")
create.hist2(read.table("DE_e_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DE_e_50_10PCT.eps")
create.hist2(read.table("DE_f_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "DE_f_50_10PCT.eps")
create.hist2(read.table("DE_g_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="g", figname = "DE_g_50_10PCT.eps")




#################################################################################
##  simulate type I error and re-calibrated power ###############################

source("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/summarize.results.R")
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160203/")
FigurePath <-"/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160203/"
showcol <- c(1, 3:8)
textsize = rep(20,20, 8, 20)


## type I error for 20,000 simulations 10% DE

create.hist2(read.table("DE_a0_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="a0", figname = "DE_a0_50_10PCT.eps")
create.hist2(read.table("DE_a_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "DE_a_50_10PCT.eps")
create.hist2(read.table("DE_c_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="c", figname = "DE_c_50_10PCT.eps")
create.hist2(read.table("DE_e_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "DE_e_50_10PCT.eps")
create.hist2(read.table("DE_f_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "DE_f_50_10PCT.eps")
create.hist2(read.table("DE_g_50_10PCT.txt")[, showcol], textsize = textsize, figure.num ="g", figname = "DE_g_50_10PCT.eps")



## recalibrated power for DE 20% VS 10%, with null 10% DE, 20000 simulations for the type I error and 1000 for power
RecalibratePower(read.table("DE_a0_50_10PCT.txt"), read.table("Power_a0_50_20VS10PCT.txt"), colnum = showcol)
RecalibratePower(read.table("DE_a_50_10PCT.txt"), read.table("Power_a_50_20VS10PCT.txt"), colnum = showcol)
RecalibratePower(read.table("DE_c_50_10PCT.txt"), read.table("Power_c_50_20VS10PCT.txt"), colnum = showcol)
RecalibratePower(read.table("DE_e_50_10PCT.txt"), read.table("Power_e_50_20VS10PCT.txt"), colnum = showcol)
RecalibratePower(read.table("DE_f_50_10PCT.txt"), read.table("Power_f_50_20VS10PCT.txt"), colnum = showcol)
RecalibratePower(read.table("DE_g_50_10PCT.txt"), read.table("Power_g_50_20VS10PCT.txt"), colnum = showcol)



## type I error for 20,000 simulations NO DE

create.hist2(read.table("NODE_a0_50.txt")[, showcol], textsize = textsize, figure.num ="a0", figname = "NODE_a0_50.eps")
create.hist2(read.table("NODE_a_50.txt")[, showcol], textsize = textsize, figure.num ="a", figname = "NODE_a_50.eps")
create.hist2(read.table("NODE_c_50.txt")[, showcol], textsize = textsize, figure.num ="c", figname = "NODE_c_50.eps")
create.hist2(read.table("NODE_e_50.txt")[, showcol], textsize = textsize, figure.num ="e", figname = "NODE_e_50.eps")
create.hist2(read.table("NODE_f_50.txt")[, showcol], textsize = textsize, figure.num ="f", figname = "NODE_f_50.eps")
create.hist2(read.table("NODE_g_50.txt")[, showcol], textsize = textsize, figure.num ="g", figname = "NODE_g_50.eps")



## recalibrated power for 10% VS 0%,  20000 simulations for the type I error and 1000 for power
RecalibratePower(read.table("NODE_a0_50.txt"), read.table("Power_a0_50_10VS0PCT.txt"), colnum = showcol, alpha_level = 0.01)
RecalibratePower(read.table("NODE_a_50.txt"), read.table("Power_a_50_10VS0PCT.txt"), colnum = showcol)
RecalibratePower(read.table("NODE_c_50.txt"), read.table("Power_c_50_10VS0PCT.txt"), colnum = showcol)
RecalibratePower(read.table("NODE_e_50.txt"), read.table("Power_e_50_10VS0PCT.txt"), colnum = showcol)
RecalibratePower(read.table("NODE_f_50.txt"), read.table("Power_f_50_10VS0PCT.txt"), colnum = showcol)
RecalibratePower(read.table("NODE_g_50.txt"), read.table("Power_g_50_10VS0PCT.txt"), colnum = showcol)

