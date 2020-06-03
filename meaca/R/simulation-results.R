
library(ggplot2)
library(reshape2)
library(xtable)
library(dplyr)
library(gridExtra)
library(grid)
library(cowplot)




qqTypeIerror <- function(data, textsize = rep(20, 4), showcol = c(1, 4:9), color1 = "#0099FF",
                         color2 = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                         "#F0E442", "#0072B2", "#D55E00")){
  
  data <- data[, showcol]
  names(data)[names(data)=="meaca"] <- "MEACA"
  dat_new <- melt(data, value.name = "Pval", variable.name = "Method")

  A <- ggplot(dat_new) +
    stat_qq(aes(sample = Pval, colour = Method, linetype = Method), distribution=qunif, size = 1,
             geom = "line") + 
    geom_abline(slope = 1, intercept = 0, color = color1) +
    theme(legend.position="top", 
          legend.direction = "horizontal",
         # legend.key = element_rect(size = 3),
        #  legend.key.size = unit(1.5, 'lines'),
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]),  
          strip.text.y = element_text(size = 8, colour = "black", face = "bold"),  # controls the facet_grid
          axis.text=element_text(size=textsize[3], face = "bold"), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    guides(fill = guide_legend(keywidth = 0.8, keyheight = 2),
           linetype=guide_legend(keywidth = 8, keyheight = 2),
           colour=guide_legend(keywidth = 8, keyheight = 2))
   A + scale_color_manual(values = color2) + 
     scale_linetype_manual(values=c("twodash", "dotted", "dotdash", "dashed", 
                                    "F1", "4C88C488", "12345678"))
  

}


#' Combines two qq plot in one figure and produce the legend in a separate figure.
#'
#' @title Uniform QQ-plot. 
#' @param dat1  simulated p-value matrix when DE % = 0.
#' @param dat2  simulated p-value matrix for DE % = 10%.
#' @param case which correlation structure to plot. 
#' @param legend if \code{TRUE}, plot the legend; otherwise just plot the qq-plot
#' @param textsize, a vector of four, controlling the text size of legened, title, axis, axis.title, respectively 
#' @return a figure
#' @export
#' 
#' 


ArrangeTypeIerror <- function(path, case, textsize = c(20,20, 8,20), showcol = c(1, 4:9), 
                              legend= F, 
                              color1= "#0099FF",
                              color2 = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                         "#F0E442", "#0072B2", "#D55E00")){
# dat1:  p value matrix for type I error, without DE
# dat2: p value matrix for type I error, with DE
# legend is used for the last case, for all figures.    
  # read the data  
  dat1 <- read.table(paste(path, "/TypeIerror_", case[1], "_0PCT.txt", sep = ""), header = T) 
  dat2 <- read.table(paste(path, "/TypeIerror_", case[1], "_10PCT.txt", sep = ""), header = T) 

  A1 <- qqTypeIerror(dat1, textsize = textsize, showcol = showcol,  color1 = color1, color2=color2)
  A2 <- qqTypeIerror(dat2, textsize = textsize, showcol = showcol, color1=color1, color2=color2)
  
  if(legend){    # plot the legend in a separate figure
    p_fig <- get_legend(A1)
    }
  
    else{  ## arrange the plots 
    p_fig <- plot_grid(A1 + theme(legend.position = "none"), 
                       A2 + theme(legend.position = "none"), ncol = 2, nrow = 1)
    }
  return(p_fig)
}
  


#' uniform qq-plot for type I error simulation.
#'
#' @title Uniform QQ-plot. 
#' @param path  where the simulated p-value matrix locates.
#' @param textsize, a vector of four, controlling the text size of legened, title, axis, axis.title, respectively 
#' @param case   which correlation structures will show
#' @param showcol the methods chosen to compare
#' @param color1 Color for the reference line
#' @param color2 Colors for different methods. It must have the same length as the that of \code{showcol}
#' @return a figure
#' @export
#' 
#' 

plot_fig1 <- function(path, textsize = c(15, 15,  10, 10), case= letters[1:5], showcol = c(1, 4:9),
                              color1= "#0099FF",
                              color2 = c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                                         "#F0E442", "#0072B2", "#D55E00")){
  

  fig_a <- ArrangeTypeIerror(path, case = case[1], textsize = textsize, showcol = showcol, color1 = color1, color2=color2, legend = F)  
  fig_b <- ArrangeTypeIerror(path, case = case[2], textsize = textsize, showcol = showcol, color1 = color1, color2=color2, legend = F)  
  fig_c <- ArrangeTypeIerror(path, case = case[3], textsize = textsize, showcol = showcol, color1 = color1, color2=color2, legend = F)  
  fig_d <- ArrangeTypeIerror(path, case = case[4], textsize = textsize, showcol = showcol, color1 = color1, color2=color2, legend = F)  
  fig_e <- ArrangeTypeIerror(path, case = case[5], textsize = textsize, showcol = showcol, color1 = color1, color2=color2, legend = F)  
  legend_a <- ArrangeTypeIerror(path, case = case[1], textsize = textsize, showcol = showcol, color1 = color1, color2=color2, legend = T)  
  
  fig <- plot_grid(fig_a, fig_b, fig_c, fig_d, fig_e, legend_a, ncol = 1, nrow =5)
  
  fig2 <- ggdraw() + draw_plot(fig, x = 0, y = 0, width= 1, height = 0.9) + 
          draw_plot(legend_a, x = 0, y =0.9, width = 0.7, height = 0.1 )
  
}


#' Summarize the power
#' 
#' @title summerize power. 
#' @param path  where the simulated p-value matrix locates.
#' @param setup if it's group I (\code{A1}) or group II (\code{A2}) simulation. 
#' @param showcol which methods to summarize.
#' @param case the correlation structure.
#' @return a table
#' @export
#' 
summarizePower<- function(path, setup, showcol, case){
  
  sim1 <- c("5VS0", "10VS0", "15VS0", "20VS0")
  sim2 <- c("15VS10", "20VS10", "25VS10", "30VS10")
  
  if (setup == "A1") {sim = sim1}
  else {sim = sim2}
  
  file <- paste(path, "/Power_", case, "_", sim[1], "PCT.txt", sep ="")
  data <- read.table(file, header=T)[, showcol]
  p_mat <- data.frame(matrix(NA, nrow = length(data), ncol = length(sim)))
  
  for (j in 1:length(sim)){
    file <- paste(path, "/Power_", case, "_", sim[j], "PCT.txt", sep ="")
   # print(file)
    data <- read.table(file, header=T)[, showcol]
    p_mat[, j] <- data
  }
  colnames(p_mat) <- paste("S", 1:4, sep = "")

  
  
  p2 <- melt(p_mat, variable.name = "alternative", value.name = "p.value")
  p3 <-p2  %>%
    group_by(alternative) %>%
    summarize(n=n(), power=mean(p.value < 0.05),sd=sd(p.value < 0.05)) %>%
    mutate(se=sd/sqrt(n),LCI=power+qnorm(0.025)*se,UCI=power+qnorm(0.975)*se)
  p3$UCI[p3$UCI >1] = 1; p3$LCI[p3$LCI<0] =0;
  p3$case <- case
  
  return(p3)
  
}




plotPower <- function(path,  setup = "A1", showcol=1, textsize = rep(20, 4), xtext = c("DE percentage", "Power")){
  
  updated_case <- c("a", "b", "c", "d", "e")  ##  If I want to change the name of correlation structures
  
  r1 <- summarizePower(path,  setup, showcol=showcol, case = "a" )
 
  
  case <- c("b", "c", "d", "e")
  for (i in 1:length(case))
  {
    r2 <- summarizePower(path, setup, showcol=showcol, case = case[i] )
    r1 <- rbind(r1, r2)
  }
  
  r1$DEpct <- rep(c(0.05, 0.1, 0.15, 0.2), 5)
  ids <- which(names(r1)=="case") 
  r1 <- as.data.frame(r1)
  names(r1)[ids] <- "Correlation_Structure"
  
  
  if (setup == "A2") {r1$DEpct <- r1$DEpct + 0.1}
  ggplot(r1, aes(x = as.factor(DEpct), y = power, group=Correlation_Structure, color = Correlation_Structure, 
                 linetype = Correlation_Structure)) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width= 0.3, size =1) + 
    geom_line(size=0.8) + 
    geom_point(size = 0.8) +
    labs(x=xtext[1], y=xtext[2]) +
    theme(legend.position="top", 
          legend.text = element_text(size = textsize[1]),
          plot.title = element_text(size = textsize[2]), 
          axis.text=element_text(size=textsize[3]), 
          axis.title=element_text(size=textsize[4],face="bold")) +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1),
           linetype=guide_legend(keywidth = 3, keyheight = 1),
           colour=guide_legend(keywidth = 3, keyheight = 1))
  
  
}


#' Plot figure 2 of the paper.
#'
#' @title Power plot with 95% error bar. 
#' @param path  where the simulated p-value matrix locates.
#' @param showcol which method to summarize (1 for "MEACA").
#' @param textsize a vector of four, controlling the text size of legened, title, axis, axis.title, respectively
#' @param xtext the title for x-axis and y-axis
#' @return a figure
#' @export
#' 
#' 

plot_fig2 <- function(path, showcol = 1, textsize = rep(20, 4), xtext = c("DE percentage", "Power")){
  
  f0 <- plotPower(path, setup = "A1", showcol = showcol, textsize=textsize, xtext = xtext)
  f1 <- f0 + theme(legend.position = "none")
  f2 <- plotPower(path, setup = "A2", showcol = showcol, textsize=textsize, xtext = xtext) + 
    theme(legend.position = "none")
  legend_1 <- get_legend(f0)
  
  ggdraw() + draw_plot(legend_1, x=0, y= 0.94, width = 0.9, height = 0.06) + 
      draw_plot(f1, x = 0, y =  0.47, width = 1, height = 0.47) + 
    draw_plot(f2, x = 0, y =  0, width = 1, height = 0.47)
}







# #' the recalibrated power
# #'
# #' @title calculate the revalibrated power
# #' @param data  simulated p-value matrix.
# #' @param textsize, a vector of four, controlling the text size of legened, title, axis, axis.title, respectively
# #' @return a list
# #' \item{adjusted_alpha}{the adjusted alpha level}
# #' \item{power}{the corresponding powers}
# #' \item{se}{its standard errors}
# #' @export
# #'
# #'

RecalibratePower <- function(data1, data2, alpha_level = 0.05, showcol = 1:8){
  ## data1: the type I error table
  ## data2: the corresponding power table
  ## alpha_level: what is the desired significance level 
  ## showcol:  choose the results from which methods you want to report
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



#' the recalibrated power
#'
#' @title calculate the actual type I error for a nominal alpha 
#' @param data  simulated p-value matrix.
#' @param alpha, a vector or a numeric value that specifies the desired type I error
#' @return the adjusted alpha value.
#' @export
#' 
typeIerror_quantile <- function(data, alpha, showcol = 1){
  
  tf <- data < alpha
  results <- colMeans(tf)[showcol]
  results
}




createPowerTable <- function(folder, showcol, setup = "A1", case = "a"){
  
  pw <- "Power_"
  typeIer <- "TypeIerror_"
  t2 <-  c("5VS0PCT", "10VS0PCT", "15VS0PCT", "20VS0PCT"); row_name <- t2
  t3 <- "_0PCTforcalibration.txt"
  if (setup == "A2") {
    t2 <- c("15VS10PCT", "20VS10PCT", "25VS10PCT", "30VS10PCT")
    row_name <- t2
    t3 <- "_10PCTforcalibration.txt"
    }
 
  file1 <- paste(folder, "/", typeIer, case, t3, sep = "" )
  file2 <- paste(folder, "/", pw, case, "_", t2, ".txt", sep ="")
    
    dat1 <- read.table(file1, header=T)
    dat2 <- read.table(file2[1], header=T)
    power1 <- RecalibratePower(dat1, dat2, showcol = showcol)
    calib_power <- power1$power
    adjusted_alpha <- power1$adjusted_alpha
      
    for(i in 2:length(file2))
    {
      dat2 <- read.table(file2[i], header=T)
      power2 <-RecalibratePower(dat1, dat2, showcol = showcol)
      calib_power <- rbind(calib_power, power2$power)
    }
     rownames(calib_power) <- paste("S", 1:4, sep  = "")
     

  results <- cbind(adjusted_alpha, t(calib_power))
  return(results)  
}




EmpiricalPower <- function(data, alpha_level = 0.05, showcol = 1:8){
  ## data: the power table
  ## alpha_level: what is the desired significance level 
  ## showcol:  choose the results from which methods you want to report
  data <- data[, showcol]
  
  power <- se <- c()
  for ( i in 1: ncol(data))
  {
    power[i] <- mean(data[, i] < alpha_level)
    se[i] <- sd(data[, i] < alpha_level)/sqrt(nrow(data))
  }
  names(power) <- colnames(data)
  names(se) <- colnames(data)
  
  return(list(power = power, se = se))
}
  
#' the table in the emperical power simulation 
#' @title create the emperical power table . 
#' @param path  where the simulated p-value matrix locates.
#' @param showcol for which methods power will be summarized  
#' @param setup if it's group I (\code{A1}) or group II (\code{A2}) simulation. 
#' @return a table
#' @export
#' 
createEmpiricalPowerTable <- function(path, showcol, setup = "A1", case = "a"){
  
  pw <- "Power_"
  typeIer <- "TypeIerror_"
  t2 <-  c("5VS0PCT", "10VS0PCT", "15VS0PCT", "20VS0PCT"); 
  row_name <- t2

    if (setup == "A2") {
    t2 <- c("15VS10PCT", "20VS10PCT", "25VS10PCT", "30VS10PCT")
    row_name <- t2
    }
 
  file <- paste(path, "/", pw, case, "_", t2, ".txt", sep ="")
    
    dat <- read.table(file[1], header=T)
    power1 <- EmpiricalPower(dat, alpha_level = 0.05, showcol = showcol)
    power <- power1$power
      
    for(i in 2:length(file))
    {
      dat <- read.table(file[i], header=T)
      power2 <-EmpiricalPower(dat, showcol = showcol)
      power <- rbind(power, power2$power)
    }
     rownames(power) <- paste("S", 1:4, sep  = "")
     

  results <- t(power)
  return(results)  

  
}  
