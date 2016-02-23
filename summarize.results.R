## the functions used to present the results 

library(ggplot2)
library(reshape2)
library(xtable)
library(dplyr)

create.hist2 <- function(data, figure.num = "", textsize = rep(20, 4), figname = "all.eps")
{
  dat_new <- melt(data, value.name = "Pval", variable.name = "Method")
  
  setEPS() 
  postscript(paste(FigurePath, figname, sep=""), width = 10, height = 7)
  
  A <- ggplot(dat_new, aes(x=Pval)) + geom_histogram(binwidth=.05, colour="black", fill="white") + 
    facet_grid(Method~ ., scales = "free")  +
    xlim(0, 1) + 
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    labs(title=figure.num)    
  print(A)
  dev.off()
}


## given a vector of norminal level, produce the actural type I error rate

TypeIerror <- function(norminal = 0.05, data){
  
  actual <- data.frame(matrix(NA, ncol = ncol(data) , nrow = length(norminal)))
  
  for ( i in 1: length(norminal))
  {
    actual[i, ] <- colSums(data <= norminal[i])/nrow(data)
  }
  colnames(actual) <- names(data)
  return(t(actual)) 
}



## the recalibrated power

RecalibratePower <- function(data1, data2, alpha_level = 0.05, colnum = 1:8){
## data1: the type I error table
## data2: the corresponding power table
## alpha_level: what is the desired significance level 
## colnum:  choose the results from which methods you want to report
  data1 <- data1[, showcol]
  data2 <- data2[, showcol]
  adjusted_alpha <- sapply(data1, quantile, alpha_level)
  
  power <- c()
  for ( i in 1: ncol(data1))
  {
    power[i] <- mean(data2[, i] < adjusted_alpha[i])
  }
  names(power) <- colnames(data1)
  names(adjusted_alpha) <- colnames(data1)
  
  return(list(adjusted_alpha = adjusted_alpha, power = power))
  }




### create the power table 


producePower <- function(size, background, coln,  alpha_level= 0.05){
## size: the three levels, 0.05, 0.1, 0.2
## background: has two levels, "BACK0", "BACK10"
    
  setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160211/")
  
  case <- c("a0", "a", "c", "e", "f", "g")
  sim1 <- c("5VS0", "10VS0", "15VS0", "20VS0")
  sim2 <- c("15VS10", "20VS10", "25VS10", "30VS10")

  if (size == 0.05){folder <- "SizePoint05"}
  else if (size == 0.1){folder <- "SizePoint1"}
  else if (size == 0.2){folder <- "SizePoint2"}
  
  if (background == "BACK0") {sim = sim1}
   else {sim = sim2}

  power <- matrix(NA, nrow = length(case), ncol = length(sim))
  for ( i in 1: length(case)){
    for (j in 1:length(sim)){
      file <- paste(folder, "/Power_", case[i], "_", sim[j], "PCT.txt", sep ="")
      data <- read.table(file)[, coln]
      power[i, j] <- mean(data < alpha_level)
    }
  }
  power <- data.frame(power)
  rownames(power) <- case
  colnames(power) <- sim
  return(power)
}


combinePvalue <- function(size, background, coln, case){
  ## size: the three levels, 0.05, 0.1, 0.2
  ## background: has two levels, "BACK0", "BACK10"
  
  setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160211/")
  
#  case <- c("a0", "a", "c", "e", "f", "g")
  sim1 <- c("5VS0", "10VS0", "15VS0", "20VS0")
  sim2 <- c("15VS10", "20VS10", "25VS10", "30VS10")
  
  if (size == 0.05){folder <- "SizePoint05"}
  else if (size == 0.1){folder <- "SizePoint1"}
  else if (size == 0.2){folder <- "SizePoint2"}
  
  if (background == "BACK0") {sim = sim1}
  else {sim = sim2}
  
  p_mat <- data.frame(matrix(NA, nrow = 1000, ncol = length(sim)))
  
    for (j in 1:length(sim)){
      file <- paste(folder, "/Power_", case, "_", sim[j], "PCT.txt", sep ="")
      data <- read.table(file)[, coln]
      p_mat[, j] <- data
    }
  colnames(p_mat) <- sim
  
  p2 <- melt(p_mat, variable.name = "alternative", value.name = "p.value")
  p3 <-p2  %>%
    group_by(alternative) %>%
    summarize(n=n(), power=mean(p.value < 0.05),sd=sd(p.value < 0.05)) %>%
    mutate(se=sd/sqrt(n),LCI=power+qnorm(0.025)*se,UCI=power+qnorm(0.975)*se)
  p3$UCI[p3$UCI >1] = 1; p3$LCI[p3$LCI<0] =0;
  p3$case <- case
  
  return(p3)
  
}



plotPower <- function(size,  background = "BACK0", coln=1, textsize = rep(20, 4), xtext = c("DE percentage", "Power")){
  
  
  r1 <- combinePvalue(size,  background, coln=coln, case = "a0" )
  
  case <- c("a", "c", "e", "f")
  for (i in 1:length(case))
  {
    r2 <- combinePvalue(size, background, coln=coln, case = case[i] )
    r1 <- rbind(r1, r2)
  }
  
  r1$DEpct <- rep(c(0.05, 0.1, 0.15, 0.2), 5)
  if (background == "BACK10") {r1$DEpct <- r1$DEpct + 0.1}
  ggplot(r1, aes(x = as.factor(DEpct), y = power, group=case,  color = case)) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width= 0.1) + 
    geom_line() + 
    geom_point() +
    labs(x=xtext[1], y=xtext[2]) +
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    # scale_y_continuous(labels = percent)  +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1),
           linetype=guide_legend(keywidth = 3, keyheight = 1),
           colour=guide_legend(keywidth = 3, keyheight = 1))
  
}

plotPower(0.05);
plotPower(0.05, background = "BACK10")


powertable <- producePower(0.05, "BACK0", coln=1)[1:5, ]
plotPower(powertable, background = "BACK0")

powertable <- producePower(0.05, "BACK10", coln=1)[1:5, ]
plotPower(powertable, background = "BACK10")





powertable <- producePower(0.05, "BACK0", coln=1)[1:5, ]

## comparing power of different methods

powerDiffMethod <- function(rown, size, background, alpha_level=0.05){
# rown: specifies which case we want to compare
# size: corresponds to the three folders
# background: whether the null is 10% or 0%  
    
  result <- data.frame(matrix(NA, nrow = 8, ncol = 4))
  for ( i in 1:8){
    result[i,] <- producePower(size, background, coln = i, alpha_level = alpha_level)[rown, ]
}
  p1 <- read.table("Power_a0_50_5VS0PCT.txt")
  powertable <- producePower(0.05, "BACK0", coln=1)[1:5, ]
  colnames(result) <- colnames(powertable)
  rownames(result) <- colnames(p1)
  return(result) 
    
}


### create two tables 
p1 <- powerDiffMethod(1, size = 0.05, background = "BACK0", alpha_level = 0.05)
p2 <- powerDiffMethod(1, size = 0.1, background = "BACK10", alpha_level = 0.05)

xtable(p1, digits = 3)
xtable(p2, digits = 3)


## use the calibrated power at level 0.05

showcol <- c(1:8)

power05 <-RecalibratePower(read.table("SizePoint05/TypeIerror_a0_0PCT_SizePoint05.txt"), read.table("SizePoint05/Power_a0_5VS0PCT.txt"), colnum = showcol)
power10 <-RecalibratePower(read.table("SizePoint05/TypeIerror_a0_0PCT_SizePoint05.txt"), read.table("SizePoint05/Power_a0_10VS0PCT.txt"), colnum = showcol)
power15 <-RecalibratePower(read.table("SizePoint05/TypeIerror_a0_0PCT_SizePoint05.txt"), read.table("SizePoint05/Power_a0_15VS0PCT.txt"), colnum = showcol)
power20 <-RecalibratePower(read.table("SizePoint05/TypeIerror_a0_0PCT_SizePoint05.txt"), read.table("SizePoint05/Power_a0_20VS0PCT.txt"), colnum = showcol)
calib_power <- rbind(power05$power, power10$power, power15$power, power20$power)
xtable(t(calib_power), digits =3)

p15 <-RecalibratePower(read.table("SizePoint1/TypeIerror_a0_10PCT_SizePoint1.txt"), read.table("SizePoint1/Power_a0_15VS10PCT.txt"), colnum = showcol)
p20 <-RecalibratePower(read.table("SizePoint1/TypeIerror_a0_10PCT_SizePoint1.txt"), read.table("SizePoint1/Power_a0_20VS10PCT.txt"), colnum = showcol)
p25 <-RecalibratePower(read.table("SizePoint1/TypeIerror_a0_10PCT_SizePoint1.txt"), read.table("SizePoint1/Power_a0_25VS10PCT.txt"), colnum = showcol)
p30 <-RecalibratePower(read.table("SizePoint1/TypeIerror_a0_10PCT_SizePoint1.txt"), read.table("SizePoint1/Power_a0_30VS10PCT.txt"), colnum = showcol)
calib_power2 <- rbind(p15$power, p20$power, p25$power, p30$power)
xtable(t(calib_power2), digits =3)



## type I error simulation plots

FigurePath <-"/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160211/SizePoint05/"
setwd(FigurePath)

showcol <- c(1, 3:8)
textsize = rep(20,20, 6, 20)

case <- "c"
typeIdata <- paste("TypeIerror_", case, "_10PCT_SizePoint1.txt", sep = "")
figure_num <- paste("TypeIerror_", case, "_10PCT", sep = "")
figure_name <- paste("TypeIerror_", case, "_10PCT.eps", sep = "")

create.hist2(read.table(typeIdata)[, showcol], textsize = textsize, figure.num =figure_num, figname = figure_name)

library(gap)

case <- c("a0", "a", "c", "e", "f")
par(mfrow = c(3, 2))
method <- 1
for ( i in case){
  typeIdata <- paste("TypeIerror_", i, "_0PCT_SizePoint05.txt", sep = "")
  x <- read.table(typeIdata)
  qqunif(x[, method], main = paste("Q-Q plot", i, "0PCT") )
}






