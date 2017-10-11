#' Get top genes
#'
#' Return top genes in DESeq results
#'
#' @param res an annotated DESeq2 results file
#' @param top  Number of top genes to display in matrix
#' @param padjust Adjusted p-value cutoff
#' @param basemean basemean cutoff
#' @param log2FC absolute value log2 fold change cutoff
#' @param sort_fc Sort by fold changes and get top n/2 up and down-regulated,
#'  default is to sort by adjusted p-value
#'
#' @return A tibble of top genes
#'
#' @author Chris Stubben
#'
#' @examples
#' data(pasilla)
#'  top_genes(pasilla$results)
#' @export

top_genes <- function(res, top=40, padjust = 0.05,  basemean, log2FC, sort_fc=FALSE){
    ## TO DO  code for DataFrame -  fix for tibbles
     x <- subset(res, padj <  padjust )
    if(!missing(basemean))  x <- subset(res, baseMean > basemean )
    if(!missing(log2FC))  x <- subset(res, abs(log2FoldChange) > log2FC )
    if(nrow(x) == 0 ) stop("No rows matching cutoffs")
    if(sort_fc){
      ## sort largest fold changes OR p-value ?
      ##x1 <- x[order(x$log2FoldChange, decreasing=TRUE), ]
      ##x2 <- x[order(x$log2FoldChange), ]
      x1 <- subset(res, log2FoldChange < 0)
      x1 <- x1[order(x1$padj), ]
      x2 <- subset(res, log2FoldChange > 0)
      x2 <- x2[order(x2$padj), ]
      x <- rbind( utils::head(x1, top/2),  utils::head(x2, top/2))
   }else{
      # sort by p-adjusted
      x <- x[order(x$padj), ]
      x <- utils::head(x, top)
   }
   x
}
