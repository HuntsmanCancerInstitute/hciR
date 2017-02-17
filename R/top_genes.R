#' Get top genes
#'
#' Return top genes in DESeq results
#'
#' @param res an annotated DESeq2 results file
#' @param top  Number of top genes to display in matrix
#' @param alpha Adjusted p-value cutoff
#' @param basemean basemean cutoff
#' @param log2FC absolute value log2 fold change cutoff
#' @param sort_fc Sort top genes by fold changes and get top n/2 up and down-regulated, default is to sort by adjusted p-value
#'
#' @return A tibble of top genes
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  top_genes(res)
#'  top_genes(res[[1]], sort_fc=TRUE)
#' }
#' @export

top_genes <- function(res, top=40, alpha = 0.05,  basemean, log2FC, sort_fc=FALSE){

    ## TO DO  code for DataFrame -  fix for tibbles
     x <- subset(res, padj <= alpha )
    if(!missing(basemean))  x <- subset(res, baseMean > basemean )
    if(!missing(log2FC))  x <- subset(res, abs(log2FoldChange) > log2FC )
    if(nrow(x) == 0 ) stop("No rows matching cutoffs")
    if(sort_fc){
      x1 <- x[order(x$log2FoldChange, decreasing=TRUE), ]
      x2 <- x[order(x$log2FoldChange), ]
      x <- rbind( head(x1, top/2),  utils::head(x2, top/2))
   }else{
      # sort by p-adjusted
      x <- x[order(x$padj), ]
      x <- head(x, top)
   }
   x
}
