#' Plot heatmap of top genes
#'
#' Plot heatmap of top genes with annoation bars, reordered branches, and color scale midpoint at zero.
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param intgroup one or more column names for pheatmap \code{annotation_col} bar
#' @param output  pheatmap or d3heatmap
#' @param palette a color palette name or vector of colors
#' @param dendsort reorder branches using \code{dendsort} package
#' @param midpoint0 center color scale midpoint at zero
#' @param border pheatmap border color, default none
#' @param \dots additional options passed to \code{pheatmap}
#'
#' @return A pheatmap
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- top_counts( res, rld)
#' plot_genes(x)
#' plot_genes(x, "Trt", palette="RdBu" )
#' plot_genes( top_counts(res, rld, sort_fc=TRUE) )
#' }
#' @export

 plot_genes <-   function( x, intgroup, output="pheatmap", palette="RdYlBu", dendsort=TRUE, midpoint0 = TRUE, border=NA,  ...){
   clrs <- palette
   if(length(palette)==1){
       # reverse divergent color palette
       if(palette %in% c("BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral")){
          clrs <- rev( colorRampPalette( brewer.pal(11, palette))(255) )
       }else{
          clrs <- colorRampPalette( brewer.pal(9, palette))(255)
       }
   }
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
   if(output == "pheatmap"){
      ## use  dendsort to reorder branches
      if(dendsort){
         callback <- function(hc, ...){ dendsort(hc)  }
         pheatmap(x, clrs, clustering_callback = callback, annotation_col=df, breaks=brks, border=border,  ...)
      }else{
         pheatmap(x, clrs, annotation_col=df, breaks=brks, border=border,  ...)
      }
   }else{
      if(dendsort){
          ##  need to flip rows to match pheatmap??
          d3heatmap(x, reorderfun= function(d,w) dendsort(d) , colors = clrs, ...)
      }else{
         d3heatmap(x, colors = clrs, ...)
      }
   }
}
