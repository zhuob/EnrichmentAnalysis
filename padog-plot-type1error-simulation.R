
plot_type1error <- function(dat_type = "test", var_col = 1:4){
  dat_name <- paste0("real-data/padog-package/padog-real-data-type1error-simulation-", 
                     dat_type,".csv")
  t1 <- read_csv(dat_name)
  print(dat_name)
  sample_test_genes <- pivot_longer(t1, cols = var_col, names_to = "method")
  
  ggplot(data = sample_test_genes, aes(sample = value)) + 
    stat_qq(distribution = qunif) + 
    geom_abline(intercept = 0, slope = 1) + 
    facet_wrap(~method)
  
}

library(tidyverse)
plot_type1error(dat_type = "resample-test-set")
plot_type1error(dat_type = "resample-back-set")
plot_type1error(dat_type = "all-genes")

plot_type1error(dat_type = "all-genes-v2", var_col = 1:7)
dat_name <- paste0("real-data/padog-package/padog-real-data-type1error-simulation-", 
                   "all-genes-v2",".csv")
t1 <- read_csv(dat_name)

plot_type1error(dat_type = "all-genes-v3-row-bootstrap", var_col = 1:9)
