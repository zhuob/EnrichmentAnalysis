
#' read the gene sets of the MsigDB format.
#'
#' @title Convert gene sets to lists
#' @param msigdb gene set ensemble downloaded from broad institute.
#' @return a list 
#' \item{total}{ number of gene sets contained.}
#' \item{size}{a numerical vector containing the size of each gene set.}
#' \item{gene_set}{ a list. The first element is the set name. From  the third element, each corresponds to a member gene of the gene set.}
#' @export
# #' @examples


read_gene_set <- function(msigdb){
  ## read the gene sets from MsigDB   
  gene_set <- readLines(msigdb)                                        # the gene sets
  temp <- gene_set
  gs.line <- list()
  set_name <- c()
  
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric") 
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
    gs.line[[i]] <- noquote(unlist(strsplit(temp[[i]], "\t")))[-(1:2)]
    set_name[i] <- noquote(unlist(strsplit(temp[[i]], "\t")))[1]
  }
  return(list(total = length(temp.size.G), size = temp.size.G, 
              set_name=set_name,  gene_set = gs.line))
}

