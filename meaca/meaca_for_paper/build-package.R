library(devtools)
library(roxygen2)
setwd("/Users/Bin/Google Drive/Study/Thesis/Correlation/EnrichmentAnalysis/meaca/meaca_for_paper/meaca")
document()
check()
build()
setwd("..")
install("meaca")

pack <- "MEQLEA"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))
