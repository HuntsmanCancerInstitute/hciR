#' Plot heatmap
#'
#' Plot pheatmap with annoation bar, reordered branches, and color scale midpoint at zero.
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param intgroup one or more column names for \code{annotation_col} bar
#' @param  color a divergent color palette name or vector of colors
#' @param midpoint0 center color scale midpoint at zero
#' @param dendsort reorder branches using \code{dendsort} package
#' @param \dots additional options passed to \code{pheatmap}
#'
#' @return A pheatmap
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- top_counts( res[[1]], rld)
#' pheatmap_genes(x)
#' pheatmap_genes(x, "Trt", border=NA )
#' pheatmap_genes( top_counts(res[[1]], rld , col_label="Name", sort_fc=TRUE)
#' }
#' @export

 pheatmap_genes <-   function( x, intgroup, color="RdYlBu",  midpoint0 = TRUE, dendsort=TRUE,  ...){

    if(length(color)==1)  color <- rev( colorRampPalette( brewer.pal(11, color))(255) )
    df <- NA
    if(!missing( intgroup)){
        df <- attr(x, "colData")[, intgroup, drop=FALSE]
       df[,1] <- paste0(df[,1], "   ")  # hack to fix right margin
    }
      x <- as.matrix(x)
      brks <- NA
      if(midpoint0){
          n <- max(abs(range(x)))
          brks <- seq(-n, n, length=255)
      }
      ## use  dendsort to reorder branches
       if(dendsort){
             callback = function(hc, ...){dendsort(hc)}
            pheatmap(x, color, clustering_callback = callback, annotation_col=df, breaks=brks,  ...)
         }else{
            pheatmap(x, color, annotation_col=df, breaks=brks,  ...)
         }
}
