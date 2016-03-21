
#' read the gene sets of the MsigDB format.
#'
#' @title Convert gene sets to lists
#' @param msigdb gene set ensemble downloaded from broad institute.
#' @return a list 
#' \item{total}{ number of gene sets contained.}
#' \item{size}{a numerical vector containing the size of each gene set.}
#' \item{gene_set}{ a list. The first element is the set name. From the third element each containing members of the gene set.}
#' @export
#' @examples


read_gene_set <- function(msigdb){
  ## read the gene sets from MsigDB   
  gene_set <- readLines(msigdb)                                        # the gene sets
  temp <- gene_set
  gs.line <- list()
  
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric") 
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line[[i]] <- noquote(unlist(strsplit(temp[[i]], "\t")))
  }
  return(list(total = length(temp.size.G), size = temp.size.G, gene_set = gs.line))
}



#' read the gene sets of the MsigDB format.
#'
# #' @title Read Cls files for treatment labels.
# #' @param msigdb gene set ensemble downloaded from broad institute.
# #' @return a list 
# #' \item{total} { number of gene sets contained.}
# #' \item{size} {a numerical vector containing the size of each gene set.}
# #' \item{gene_set} { a list. The first element is the set name. From the third element each containing members of the gene set.}
# #' @export
# #' @examples


GSEA.ReadClsFile <- function(file = "NULL") { 
  #
  # Reads a class vector CLS file and defines phenotype and class labels vectors for the samples in a gene expression file (RES or GCT format)
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  cls.cont <- readLines(file)
  num.lines <- length(cls.cont)
  class.list <- unlist(strsplit(cls.cont[[3]], " "))
  s <- length(class.list)
  t <- table(class.list)
  l <- length(t)
  phen <- vector(length=l, mode="character")
  phen.label <- vector(length=l, mode="numeric")
  class.v <- vector(length=s, mode="numeric")
  for (i in 1:l) {
    phen[i] <- noquote(names(t)[i])
    phen.label[i] <- i - 1
  }
  for (i in 1:s) {
    for (j in 1:l) {
      if (class.list[i] == phen[j]) {
        class.v[i] <- phen.label[j]
      }
    }
  }
  return(list(phen = phen, class.v = class.v))
}


#' read the expression data of Gct format.
#'
#' @title Convert gene sets to lists
#' @param file  the expression data files. 
# #' @return a list 
# #' \item{total} { number of gene sets contained.}
# #' \item{size} {a numerical vector containing the size of each gene set.}
# #' \item{gene_set} { a list. The first element is the set name. From the third element each containing members of the gene set.}
# #' @export
# #' @examples



GSEA.Gct2Frame2 <- function(filename = "NULL") { 
  #
  # Reads a gene expression dataset in GCT format and converts it into an R data frame
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  content <- readLines(filename)
  content <- content[-1]
  content <- content[-1]
  col.names <- noquote(unlist(strsplit(content[1], "\t")))
  col.names <- col.names[c(-1, -2)]
  num.cols <- length(col.names)
  content <- content[-1]
  num.lines <- length(content)
  
  
  row.nam <- vector(length=num.lines, mode="character")
  row.des <- vector(length=num.lines, mode="character")
  m <- matrix(0, nrow=num.lines, ncol=num.cols)
  
  for (i in 1:num.lines) {
    line.list <- noquote(unlist(strsplit(content[i], "\t")))
    row.nam[i] <- noquote(line.list[1])
    row.des[i] <- noquote(line.list[2])
    line.list <- line.list[c(-1, -2)]
    for (j in 1:length(line.list)) {
      m[i, j] <- as.numeric(line.list[j])
    }
  }
  ds <- data.frame(m)
  names(ds) <- col.names
  row.names(ds) <- row.nam
  return(ds)
}
