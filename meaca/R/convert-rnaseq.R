
#' Perform count matrix transformation using edgeR procedure 
#'
#' @param y the expression count matrix, columns being samples, rows being genes
#' @param group a vector of treatment label (e.g., 0 for control, 1 for treatment).
#'
#' @return a matrix of transformed data 
#' @export
#' @seealso \code{\link[edgeR]{camera.DGEList}}
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
#' y2 <- transform_count_edgeR(y = y0, group = group)
#' camera(y2, iset1, design)
#' 
transform_count_edgeR <- function(y, group){
  
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




#' @title  Perform variance stabilizing transformation for the count matrix data
#' @details It is absolutely critical that the columns of the count matrix and
#'   the rows of the column data (information about samples) are in the same
#'   order. For more details, see \url{https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html}
#' @param y the expression count matrix, columns being samples, rows being genes
#' @param group a vector of treatment label (e.g., 0 for control, 1 for treatment).
#' @param ... other parameters used in \code{\link[DESeq2]{varianceStabilizingTransformation}}
#' @return a matrix of transformed data 
#' @export
#'
#' @examples
#' y <- matrix(rbinom(6000, 20, 0.4), nrow = 1000)
#' group <- c(0, 0, 0, 1, 1, 1)
#' yr <- transform_count_vst(y = y, group = group)

transform_count_vst <- function(y, group, ...){
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html  
  
  if(ncol(y) != length(group)){
    stop("number of samples must be equal to number of group labels")
  }
  
  samples <- data.frame(condition = as.factor(group))
  rownames(samples) <- paste("sample", 1:length(group), sep = "")
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = y, colData = samples, 
                                        design = ~ condition)  

  vsd <- DESeq2::varianceStabilizingTransformation(object = dds, ...)
  
  trans_dat <- SummarizedExperiment::assay(vsd)
  
  return(trans_dat)
}

