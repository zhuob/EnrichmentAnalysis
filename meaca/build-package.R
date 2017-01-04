library(devtools)
library(roxygen2)
# create("meaca")
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/meaca/meacaPackage")
document()
check()
build()
setwd("..")
install("meacaPackage")
install.packages("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/meaca/meaca_0.0.0.9000.tar.gz", repos = NULL, type="source")
 
pack <- "meaca"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))
