#' Get a count matrix of top genes
#'
#' Join top genes in DESeq results to rlog or other counts, mainly for plotting gene heatmaps.
#'
#' @param res an annotated DESeq2 results file
#' @param rld rlog or other counts in DESeqTransform object
#' @param by join count rownames by this column number in results, default 1
#' @param top  Number of top genes to display in matrix
#' @param basemean basemean cutoff
#' @param log2FC absolute value log2 fold change cutoff
#' @param col_names a column name in colData(rld) to use as column names
#' @param row_names a column name in results to use as row names, default gene_name
#' @param difference Subtract the row mean from counts to display the difference from the gene's average
#' @param sort_fc Sort top genes by fold changes and get top n/2 up and down-regulated, default is to sort by adjusted p-value
#'
#' @return A tibble with colData attribute
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- top_counts(res, rld)
#'  pheatamp(x)
#'  pheatmap_genes(x)
#'  # or with annotation bar, re-ordered branches and color scale midpoint at zero
#'  pheatmap_genes(x, "Trt", border=NA )
#' }
#' @export

top_counts <- function(res, rld, top=40, alpha = 0.05,  basemean, log2FC,
   col_names, trt, row_names="gene_name", difference=TRUE, sort_fc=FALSE, ...){
    if(!class(rld) == "DESeqTransform") stop("rld shoud be a DESeqTransform")
    rldx <- assay(rld)
    colx <- data.frame( colData(rld) , drop=FALSE)
    ## rename columns in heatmap
    if(!missing(col_names)){
         n <- as.character( colData(rld)[[col_names]] )
         if(is.null(n)) stop("No column matching ", col_names, " in colData(rld)")
         colnames(rldx) <- n
         rownames(colx) <- n
      }
    ## CUTOFFS
     x <- subset(res, padj <= alpha )
    if(!missing(basemean))  x <- subset(res, baseMean > basemean )
    if(!missing(log2FC))  x <- subset(res, abs(log2FoldChange) > log2FC )
   if(nrow(x) == 0 ) stop("No rows matching cutoffs")
   if(sort_fc){
      x1 <- x[order(x$log2FoldChange, decreasing=TRUE), ]
      x2 <- x[order(x$log2FoldChange), ]
      x <- rbind( head(x1, top/2),  head(x2, top/2))
   }else{
      # sort by p-adjusted
      x <- x[order(x$padj), ]
      x <- head(x, top)
   }
     ## match column 1 in results to count rownames
     n <- match(x[[ by ]], rownames(rldx))
     if(all(is.na(n))) stop("Column ", by, " in results and rownames in counts do not match")
     if(any(is.na(n))) stop(sum(is.na(n)) , " result rows not in count matrix" )
     mat <- rldx[n, ]
     # subtract the row mean
     if(difference) mat <- mat - rowMeans(mat)
     ## gene name by default as id
      # add_column( as_tibble(x), id=rownames(x) , .before=1)
     mat <- bind_cols( tibble(id=x[[ row_names ]]), as_tibble(mat) )
    attr(mat, "colData") <- colx
    mat
}
