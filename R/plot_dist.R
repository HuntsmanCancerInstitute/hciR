#' Plot heatmap of sample distances
#'
#' Plot sample distances using rlog or other counts.
#'
#' @param rld a matrix or DESeqTransform with assay slot
#' @param intgroup one or more column names in colData(rld) for pheatmap \code{annotation_col}
#' @param output  pheatmap or d3heatmap
#' @param palette a color palette name or vector of colors
#' @param diagNA set the diagonal to NA, default TRUE
#' @param font_size x and yaxis font size for d3heatmap
#' @param \dots Additional options like colors passed to \code{d3heatmap} or \code{pheatmap}
#'
#' @return A pheatmap or d3heatmap
#'
#' @author Chris Stubben
#'
#' @examples
#' data(pasilla_dds)
#' rld <- rlog(pasilla_dds)
#' plot_dist(rld, c("condition", "type") )
#' plot_dist(rld, palette="Blues", diagNA=FALSE)
#' plot_dist(rld, output = "d3")
#' @export

 plot_dist <-   function( rld, intgroup, output="pheatmap", palette="RdYlBu", diagNA = TRUE, font_size=12, ...){
     if(class(rld)[1] == "DESeqTransform"){
        d1 <-   stats::dist(t( SummarizedExperiment::assay(rld) ))
     }else{
        d1 <- stats::dist(t(rld))
     }
     sample_dist <- as.matrix(d1)
     ## coloring the diagonal often skews the color scale
     if(diagNA) diag(sample_dist) <- NA
     ## clrs
     clrs <- palette

     if(length(palette) == 1){
        # reverse divergent color palette
        if(palette %in% c("BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral")){
           clrs <- rev( grDevices::colorRampPalette( RColorBrewer::brewer.pal(11, palette))(255) )
        }else{
           clrs <- grDevices::colorRampPalette( RColorBrewer::brewer.pal(9, palette))(255)
        }
     }
     if(output == "pheatmap"){
        ##  dendsort to reorder branches
        callback <- function(hc, ...){dendsort::dendsort(hc)}
         df <- NA
        if(!missing(intgroup)){
              df <- as.data.frame( SummarizedExperiment::colData(rld)[, intgroup, drop=FALSE])
              df[,1] <- paste0(df[,1], "  ")  # hack to fix right margin
        }
        pheatmap::pheatmap(sample_dist, color=clrs, clustering_callback = callback,
                 clustering_distance_rows=d1,
                 clustering_distance_cols=d1, annotation_col=df, ...)
     }else{
          ##  dendsort to reorder branches
        dend <- dendsort::dendsort( stats::as.dendrogram( stats::hclust(d1) ) , isReverse=TRUE)
        d3heatmap::d3heatmap(sample_dist, Rowv=dend, Colv=dend, colors = clrs,
            xaxis_font_size = font_size, yaxis_font_size = font_size, ...)
     }
}