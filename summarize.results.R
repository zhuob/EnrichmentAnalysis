## the functions used to present the results 

library(ggplot2)
library(reshape2)

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

RecalibratePower <- function(data1, data2, alpha_level = 0.05, colnum = showcol){
  
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


plotPower <- function(powertable,   textsize = rep(20, 4), xtext = c("DE percentage", "Power")){
  powertable$CorStructure <- rownames(powertable)
  powertable2 <- melt(powertable, id.vars = "CorStructure", value.name = "Power")
  powertable2$DEpct <- rep(c(0.05, 0.1, 0.15, 0.2), each= nrow(powertable))
  
  ggplot(powertable2,  aes(x= DEpct, y = Power, color = CorStructure,  linetype = CorStructure)) + 
    geom_line(size = 1.5) +
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


powertable <- producePower(0.05, "BACK0", coln=1)[1:5, ]
plotPower(powertable)
