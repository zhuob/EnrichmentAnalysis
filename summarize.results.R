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






