rm(list = ls())
library(devtools)
library(roxygen2)
# create("meaca")
setwd("meaca")
devtools::load_all()
# devtools::test()
devtools::document()
devtools::check()

Sys.getenv("PATH")
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "C:\\Program Files\\MiKTeX 2.9\\miktex\\bin\\x64", sep=.Platform$path.sep))
devtools::build(manual = TRUE)

setwd("..")
install("meaca")
install.packages("meaca/meaca_0.0.0.9000.tar.gz", repos = NULL, type="source")
 
pack <- "meaca"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))
