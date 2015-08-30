arab <- readRDS("~/Google Drive/Study/Thesis/Correlation/Data/arabidopsis.21.rds")
trt <- readRDS("~/Google Drive/Study/Thesis/Correlation/Data/treatment.21.rds")

source("SimulateLabData.R")
library(NBPSeq)
lab_index <- c(1, 3, 6, 9,12, 14, 16, 19, 21)         # these lab data look similar
arab_1 <- arab[[1]]                                   # the first data set
group_mean <- group.mean(arab_1, trt[[1]])            # calculate seperate mean of each group

max_l <- apply(group_mean, 1, max)                    # the larger of mean for trt/control
min_l <- apply(group_mean, 1, min)                    # the smaller of mean for trt/control

disp_mat <- matrix(NA, nrow=nrow(y), ncol=length(arab))

subset_gene <- names(which(min_l > 90 & max_l < 110)) # select the genes
nb_data <- prepare.nb.data(arab_1)                    # normalization
s <- nb_data$eff.lib.sizes                            # effective library size
x <- model.matrix(~trt[[1]])                          # design matrix
y <- arab_1[which(row.names(arab_1) %in% subset_gene), ] # the subset to be considered
y <- as.matrix(y)
for ( j in 1: nrow(y))
{
  ## disp <- estimate.genewise.disp.mapl(temp_dat, x)
  disp_mat[j, 1] = fit.nb.glm.1u(y[j, ], s, x)$phi      # estimate genewise dispersion
}


for (k in 2:length(arab)) {
  arab_k <- arab[[k]]
  nb_data <- prepare.nb.data(arab_k)                    # normalization
  s <- nb_data$eff.lib.sizes                            # effective library size
  x <- model.matrix(~trt[[k]])                          # design matrix
  y <- arab_k[which(row.names(arab_k) %in% subset_gene), ] # the subset to be considered
  y <- as.matrix(y); disp <- c()
  for ( j in 1: nrow(y))
  {
    if (sum(y[j, ] ) < 30) { disp_mat[j, k] <- NA}        # when all expression values are 0
  #  {disp[j] <- NA}
      else {
      #  disp[j] <- fit.nb.glm.1u(y[j, ], s, x)$phi  
    ## disp <- estimate.genewise.disp.mapl(temp_dat, x)
     disp_mat[j, k] = fit.nb.glm.1u(y[j, ], s, x)$phi      # estimate genewise dispersion
    }
  }
  
}

## if MA plot of different labs have similar level, then we can borrow information across labs
matplot(t(log(disp_mat)), type="l")

dispersion <- as.vector(disp_mat)
labs <- rep(1:length(lab_index), each= length(subset_index))
# there is lab effect, different labs tend to estimate different levels of dispersion
summary(lm(dispersion~ as.factor(labs))) 






