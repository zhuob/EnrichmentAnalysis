
library(ggplot2)
library(reshape2)
library(xtable)
library(dplyr)
library(gridExtra)
library(grid)


create.hist2 <- function(data, figure.num = "", textsize = rep(20, 4), figname = "all.eps", print.figure= T)
{
  dat_new <- melt(data, value.name = "Pval", variable.name = "Method")
  
  if(print.figure)   ## show the figure in R window
  {
    
    A <- ggplot(dat_new, aes(x=Pval)) + geom_histogram(binwidth=.05, colour="black", fill="white") + 
      facet_grid(Method~ ., scales = "free")  +
      xlim(0, 1) + 
      theme(legend.position="top", 
            legend.text = element_text(size = textsize[1]),
            plot.title = element_text(size = textsize[2]),  
            strip.text.y = element_text(size = 12, colour = "black", face = "bold"),  # controls the facet_grid
            axis.text=element_text(size=textsize[3]), 
            axis.title=element_text(size=textsize[4],face="bold")) +
      labs(title=figure.num, x= "P-values")    
    print(A)
    
  }
  
  else{            # save the figure in a specified
  setEPS() 
  postscript(paste(FigurePath, figname, sep=""), width = 10, height = 6)
  
  A <- ggplot(dat_new, aes(x=Pval)) + geom_histogram(binwidth=.05, colour="black", fill="white") + 
    facet_grid(Method~ ., scales = "free")  +
    xlim(0, 1) + 
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]),  
          strip.text.y = element_text(size = 12, colour = "black", face = "bold"),  # controls the facet_grid
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    labs(title=figure.num, x= "P-values")       
  print(A)
  dev.off()
  }
}



## maybe a qqplot would be better

qqTypeIerror <- function(data, textsize = rep(20, 4)){
  dat_new <- melt(data, value.name = "Pval", variable.name = "Method")

  A <- ggplot(dat_new) +
    stat_qq(aes(sample = Pval, colour = Method, linetype = Method), dist=qunif, size = 1,
             geom = "line") + 
    geom_abline(slope = 1, intercept = 0) +
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]),  
          strip.text.y = element_text(size = 8, colour = "black", face = "bold"),  # controls the facet_grid
          axis.text=element_text(size=textsize[3], face = "bold"), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    guides(fill = guide_legend(keywidth = 0.8, keyheight = 1),
           linetype=guide_legend(keywidth = 3, keyheight = 1),
           colour=guide_legend(keywidth = 3, keyheight = 1))
  return (A)

}


### Create  qq plot for type I error.   
ArrangeTypeIerror <- function(dat1, dat2, textsize = rep(20, 4), case,  legend= F){
# dat1:  p value matrix for type I error, without DE
# dat2: p value matrix for type I error, with DE
# legend is used for the last case, for all figures.    
  
  # this function is used to extract the legend
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  A1 <- qqTypeIerror(dat1,  textsize = textsize)
  A2 <- qqTypeIerror(dat2, textsize = textsize)
  
  if(legend){    # plot the legend in a separate figure
    mylegend<-g_legend(A1)
    p_fig <- grid.arrange(mylegend , nrow=2,heights=c(3, 1))
      }
  
    else{  ## arrange the plots 
      p_fig <- grid.arrange(arrangeGrob(A1 + theme(legend.position="none"),
                           A2 + theme(legend.position="none"),
                          nrow=1), 
                          nrow=2,   # one for the plot, and one for the label
                          heights=c(10, 1),
                          main=textGrob(paste("(", case, ")", sep =""), vjust= 1,hjust=0,
                                        gp=gpar(fontsize=20, fontface = "bold")))
    }
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
  
  power <- se <- c()
  for ( i in 1: ncol(data1))
  {
    power[i] <- mean(data2[, i] < adjusted_alpha[i])
    se[i] <- sd(data2[, i] < adjusted_alpha[i])/sqrt(nrow(data2))
  }
  names(power) <- colnames(data1)
  names(adjusted_alpha) <- colnames(data1)
  names(se) <- colnames(data1)
  
  return(list(adjusted_alpha = adjusted_alpha, power = power, se = se))
}




### create the power table 


producePower <- function(size, background, coln,  alpha_level= 0.05){
  ## size: the three levels, 0.05, 0.1, 0.2
  ## background: has two levels, "BACK0", "BACK10"
  
  setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/Share/Simulation/SimulationPower20160211/")
  
  case <- c("a0", "a", "c", "e", "f")
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






summarizePower<- function(size, background, coln, case){
  ## size: the three levels, 0.05, 0.1, 0.2
  ## background: has two levels, "BACK0", "BACK10"
  
  #  case <- c("a0", "a", "c", "e", "f", "g")
  sim1 <- c("5VS0", "10VS0", "15VS0", "20VS0")
  sim2 <- c("15VS10", "20VS10", "25VS10", "30VS10")
  
  if (size == 0.05){folder <- "SizePoint05"}
  else if (size == 0.1){folder <- "SizePoint1"}
  else if (size == 0.2){folder <- "SizePoint2"}
  
  if (background == "BACK0") {sim = sim1}
  else {sim = sim2}
  
  p_mat <- data.frame(matrix(NA, nrow = 10000, ncol = length(sim)))
  
  for (j in 1:length(sim)){
    file <- paste(folder, "/Power_", case, "_", sim[j], "PCT.txt", sep ="")
    print(file)
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
  
  updated_case <- c("a", "b", "c", "d", "e")  ##  If I want to change the name of correlation structures
  
  r1 <- summarizePower(size,  background, coln=coln, case = "a0" )
 
  
  case <- c("a", "c", "e", "f")
  for (i in 1:length(case))
  {
    r2 <- summarizePower(size, background, coln=coln, case = case[i] )
    r1 <- rbind(r1, r2)
  }
  
  r1$DEpct <- rep(c(0.05, 0.1, 0.15, 0.2), 5)
  r1$CorStruct[r1$case=="a0"] <- updated_case[1]
  r1$CorStruct[r1$case=="a"] <- updated_case[2]
  r1$CorStruct[r1$case=="c"] <- updated_case[3]
  r1$CorStruct[r1$case=="e"] <- updated_case[4]
  r1$CorStruct[r1$case=="f"] <- updated_case[5]
  
  
  
  if (background == "BACK10") {r1$DEpct <- r1$DEpct + 0.1}
  ggplot(r1, aes(x = as.factor(DEpct), y = power, group=case,  color = CorStruct, linetype = CorStruct)) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width= 0.3, size =1) + 
    geom_line(size=0.8) + 
    geom_point(size = 0.8) +
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




## comparing power (not calibrated) of different methods

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


## the table in the power simulation 
createPowerTable <- function(size = 0.1, showcol =  c(1, 4, 5, 7, 8), background = "BACK0", case = "a0"){
  if (size == 0.05){folder <- "SizePoint05"}
  else if (size == 0.1){folder <- "SizePoint1"}
  else if (size == 0.2){folder <- "SizePoint2"}
  
  
  pw <- "Power_"
  typeIer <- "TypeIerror_"
  t2 <-  c("5VS0PCT", "10VS0PCT", "15VS0PCT", "20VS0PCT"); row_name <- t2
  t3 <- "_0PCT.txt"
  if (background == "BACK10") {
    t2 <- c("15VS10PCT", "20VS10PCT", "25VS10PCT", "30VS10PCT")
    row_name <- t2
    t3 <- "_10PCT.txt"
    }
 
  file1 <- paste(folder, "/", typeIer, case, t3, sep = "" )
  file2 <- paste(folder, "/", pw, case, "_", t2, ".txt", sep ="")
    
    dat1 <- read.table(file1)
    dat2 <- read.table(file2[1])
    power1 <- RecalibratePower(dat1, dat2, colnum = showcol)
    calib_power <- power1$power
    se <- power1$se
      
    for(i in 2:length(file2))
    {
      dat2 <- read.table(file2[i])
      power2 <-RecalibratePower(dat1, dat2, colnum = showcol)
      calib_power <- rbind(calib_power, power2$power)
      se <- rbind(se, power2$se)
    }
     
   rownames(se) <- rownames(calib_power) <- row_name

  
 
  export_tab <- se
  for( i in 1:nrow(se)){
    for ( j in 1:ncol(se))
    {export_tab[i, j] <- paste(format(round(calib_power[i, j], 3), nsmall=3), 
                               "(", format(round(se[i, j], 3), nsmall =3), ")", sep = "")
    }
  }
  
  xtable(t(export_tab))
  
}
  
  
