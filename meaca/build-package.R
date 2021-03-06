rm(list = ls())
library(devtools)
library(roxygen2)
# create("meaca")
devtools::load_all("./meaca")
# devtools::test()
devtools::document("./meaca")
devtools::check("./meaca")

Sys.getenv("PATH")
devtools::build("./meaca", manual = TRUE)

# devtools::install_local("../meaca_0.1.1.tar.gz", dependencies = NA, upgrade = "never")
install.packages("meaca_0.2.3.tar.gz", repos = NULL, type = "source")
 
pack <- "meaca"
path <- find.package(pack)
if(file.exists(paste0(pack, ".pdf"))){
  file.remove(paste0(pack, ".pdf"))
}
system(paste(file.path(R.home("bin"), "R"), "CMD", "Rd2pdf", shQuote(path)))

