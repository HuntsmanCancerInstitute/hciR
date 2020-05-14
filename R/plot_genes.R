#' Plot heatmap of top genes
#'
#' Plot heatmap of top genes with annoation bars, reordered branches, and color
#' scale midpoint at zero.
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param intgroup one or more column names for pheatmap \code{annotation_col}
#' bar
#' @param output  pheatmap or d3heatmap
#' @param palette RColorBrewer palette name, vector of colors, or "RdGn" for
#' Red-Green color scale.
#' @param dendsort reorder branches using \code{dendsort} package
#' @param scale  scale values, default none.  diff will substract the row mean from
#' each value. Other options are none, row and column as described in
#' \code{heatmap}
#' @param midpoint0 if scale = diff, then center color scale at zero
#' @param max_scale if scale = diff, then max value for color scale, default is
#' max(abs(range(x)))
#' @param border pheatmap border color, default NA
#' @param cluster_cols cluster columns
#' @param \dots additional options passed to \code{pheatmap} or \code{d3heatmap}
#'
#' @return A pheatmap or d3heatmap
#'
#' @author Chris Stubben
#'
#' @examples
#' x <- top_counts(pasilla$results, pasilla$rlog)
#' plot_genes(x, c("condition", "type"), scale="row", annotation_names_col=FALSE)
#' plot_genes(x, output = "d3")
#' @export

plot_genes <-  function( x, intgroup, output="pheatmap", palette="RdBu",
  dendsort=TRUE, scale="none", midpoint0 = TRUE, max_scale = NA, border=NA,
  cluster_cols=TRUE, ...){
   clrs <- palette
   if(length(clrs)==1) clrs <- palette255(clrs)
   df <- NA
   if(!missing( intgroup)){
       df <- attr(x, "colData")[, intgroup, drop=FALSE]
       if(!cluster_cols){
          if(!is.factor(df[[1]])) df[,1] <- factor(df[,1], levels = unique(df[[1]]))
       }
       for(i in 1:ncol(df)){
             # hack to fix right margin by padding with spaces
            if(is.character(df[,i])) df[, i] <- as.factor(df[, i])
            if(is.factor(df[,i])){
                levels(df[, i]) <- paste0(levels(df[, i]), "    ")
            }
       }
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
         callback <- function(hc, ...){ dendsort::dendsort(hc) }
         pheatmap::pheatmap(x, clrs, clustering_callback=callback, scale=scale,
            annotation_col=df, breaks=brks, border_color=border, cluster_cols = cluster_cols, ...)
      }else{
         pheatmap::pheatmap(x, clrs, annotation_col=df, scale=scale,
             breaks=brks, border_color=border, cluster_cols = cluster_cols,  ...)
      }
   }else{
      ## Breaks for d3heatmap ?  or scale = "row"
      if(dendsort){
         ##  need to flip rows to match pheatmap??
         d3heatmap::d3heatmap(x, reorderfun=function(d,w) dendsort::dendsort(d),
             colors = clrs, scale=scale, ...)
      }else{
         d3heatmap::d3heatmap(x, colors = clrs, scale=scale, ...)
      }
   }
}
