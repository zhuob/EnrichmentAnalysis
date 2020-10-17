
#' Convert RNA-Seq count data to format that MEACA requires 
#'
#' @param y the expression count matrix, columns being samples, rows being genes
#' @param group a vector of treatment label (e.g., 0 for control, 1 for treatment).
#'
#' @return a matrix of transformed data
#' @export
#' @seealso [edgeR::camera.DGEList()]
#' @examples
#' mu <- matrix(10, 100, 4)
#' group <- factor(c(0,0,1,1))
#' design <- model.matrix(~group)
#' set.seed(123)
#' library(edgeR)
#' y0 <- matrix(rnbinom(100*4, mu=mu, size=10),100,4)
#' y <- DGEList(counts=y0, group=group)
#' y <- estimateDisp(y, design)
#' 
#' iset1 <- 1:10
#' camera.DGEList(y, iset1, design)
#' 
#' # the Pvalue should be the same 
#' y2 <- transform_rna_seq(y = y0, group = group)
#' camera(y2, iset1, design)
#' 
transform_rna_seq <- function(y, group){
  
  if(ncol(y) != length(group)){
    stop("number of samples must be equal to number of group labels")
  }
  
  design <- model.matrix(~group)
  # create a DGEList object
  y <- edgeR::DGEList(counts = y, group = group)
  # estimate dispersions
  y <- edgeR::estimateDisp(y, design)
  
  y <- edgeR:::.zscoreDGE(y = y, design = design, contrast = ncol(design))
  
  return(y)
}

