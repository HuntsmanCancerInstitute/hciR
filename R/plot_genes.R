#' Plot heatmap of top genes
#'
#' Plot heatmap of top genes with annoation bars, reordered branches, and color scale midpoint at zero.
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param intgroup one or more column names for pheatmap \code{annotation_col} bar
#' @param output  pheatmap or d3heatmap
#' @param palette a color palette name or vector of colors, diverging color palettes are reversed
#' @param dendsort reorder branches using \code{dendsort} package
#' @param scale  scale values, default diff will substract the row mean from each value.
#'   Other options are none, row and column as described in \code{heatmap}
#' @param midpoint0 if scale = diff, then center color scale at zero
#' @param max_scale if scale = diff, then max value for color scale, default is max(abs(range(x)))
#' @param border pheatmap border color, default NA
#' @param \dots additional options passed to \code{pheatmap} or \code{d3heatmap}
#'
#' @return A pheatmap or d3heatmap
#'
#' @author Chris Stubben
#'
#' @examples
#' data(pasilla)
#'  x <- top_counts(pasilla$results, pasilla$rlog)
#' plot_genes(x)
#' plot_genes(x, "condition", palette="RdBu")
#' plot_genes(top_counts(pasilla$results, pasilla$rlog, sort_fc=TRUE) )
#' plot_genes(top_counts(pasilla$results, pasilla$rlog, top=200), output = "d3")
#' @export

 plot_genes <-  function( x, intgroup, output="pheatmap", palette="RdBu", dendsort=TRUE,
      scale="diff", midpoint0 = TRUE, max_scale = NA, border=NA,    ...){
   clrs <- palette
   if(length(palette)==1){
       # reverse divergent color palette
       if(palette %in% c("BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral")){
          clrs <- rev( grDevices::colorRampPalette( RColorBrewer::brewer.pal(11, palette))(255) )
       }else{
          clrs <- grDevices::colorRampPalette( RColorBrewer::brewer.pal(9, palette))(255)
       }
   }
   df <- NA
   if(!missing( intgroup)){
       df <- attr(x, "colData")[, intgroup, drop=FALSE]
       for(i in 1:ncol(df))  df[,i] <- paste0(df[,i], "    ")  # hack to fix right margin
   }
   ## convert tibble to matrix
   x <- as_matrix(x)
   brks <- NA
   if(scale == "diff"){
      # subtract the row mean ...
      x <- x - rowMeans(x)
      scale <- "none"   #for pheatmap
      if(midpoint0){
          n <- max_scale
          if(is.na(n)) n <- max(abs(range(x)))
          brks <- seq(-n, n, length=255)
      }
   }
   if(output == "pheatmap"){
      ## use  dendsort to reorder branches
      if(dendsort){
         callback <- function(hc, ...){ dendsort::dendsort(hc)  }
         pheatmap::pheatmap(x, clrs, clustering_callback = callback, annotation_col=df,
               scale=scale, breaks=brks, border_color=border,  ...)
      }else{
         pheatmap::pheatmap(x, clrs, annotation_col=df, scale=scale, breaks=brks, border_color=border,  ...)
      }
   }else{
      ## Breaks for d3heatmap ?  or scale = "row"
      if(dendsort){
          ##  need to flip rows to match pheatmap??
          d3heatmap::d3heatmap(x, reorderfun= function(d,w) dendsort::dendsort(d), colors = clrs, scale=scale, ...)
      }else{
         d3heatmap::d3heatmap(x, colors = clrs, scale=scale, ...)
      }
   }
}
