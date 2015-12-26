
library(ggplot2)
library(reshape2)

create.hist2 <- function(data, figure.num, textsize = rep(20, 4))
{

  dat_new <- melt(data, value.name = "Pval", variable.name = "Method")
  

  ggplot(dat_new, aes(x=Pval)) + geom_histogram(binwidth=.05, colour="black", fill="white") + 
    facet_grid(Method~ .)  +
    xlim(0, 1) + 
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
          labs(title=figure.num)    
  
}



data <- read.table("DE_a_50.txt") 

library(reshape2)

   #theme(strip.text.x = element_text(size = 15, colour = "red", angle = 90))
  theme(legend.text = element_text(size = 20))


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


create.hist2(read.table("DE_a0_50.txt"), figure.num = "DE_a0_50")
create.hist2(read.table("DE_a_50.txt"), figure.num = "DE_a_50")
create.hist2(read.table("DE_b_50.txt")[, 1:6], figure.num = "DE_b_50")
create.hist2(read.table("DE_c_50.txt")[, 1:6], figure.num = "DE_c_50")
create.hist2(read.table("DE_d_50.txt")[, 1:6], figure.num = "DE_d_50")
create.hist2(read.table("DE_e_50.txt")[, 1:6], figure.num = "DE_e_50")
create.hist2(read.table("DE_f_50.txt")[, 1:6], figure.num = "DE_f_50")
create.hist2(read.table("DE_g_50.txt")[, 1:6], figure.num = "DE_g_50")
create.hist2(read.table("DE_h_50.txt")[, 1:6], figure.num = "DE_h_50")



create.hist2(read.table("NO_DE_a0_50.txt"), figure.num = "NO_DE_a0_50")
create.hist2(read.table("NO_DE_a_50.txt"), figure.num = "NO_DE_a_50")
create.hist2(read.table("NO_DE_b_50.txt"), figure.num = "NO_DE_b_50")
create.hist2(read.table("NO_DE_c_50.txt"), figure.num = "NO_DE_c_50")
create.hist2(read.table("NO_DE_d_50.txt"), figure.num = "NO_DE_d_50")
create.hist2(read.table("NO_DE_e_50.txt"), figure.num = "NO_DE_e_50")
create.hist2(read.table("NO_DE_f_50.txt"), figure.num = "NO_DE_f_50")
create.hist2(read.table("NO_DE_g_50.txt"), figure.num = "NO_DE_g_50")
create.hist2(read.table("NO_DE_h_50.txt"), figure.num = "NO_DE_h_50")


