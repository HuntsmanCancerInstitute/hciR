#' Plot heatmap of sample distances
#'
#' Plot sample distances using rlog or other counts.
#'
#' @param rld a matrix or DESeqTransform with assay slot
#' @param intgroup a character vector of names in colData(rld) for pheatmap \code{annotation_col}
#' @param output  pheatmap or d3heatmap
#' @param palette RColorBrewer palette name
#' @param revclrs reverse color palette
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
#' plot_sample_dist(rld, c("condition", "type") )
#' plot_sample_dist(rld, palette="Blues", revclrs=TRUE, diagNA=FALSE)
#' plot_sample_dist(rld, output = "d3")
#' @export

 plot_sample_dist <-   function( rld, intgroup, output="pheatmap", palette="RdYlBu", revclrs= FALSE, diagNA = TRUE, font_size=12, ...){
     if(class(rld)[1] == "DESeqTransform"){
        d1 <-   dist(t( assay(rld) ))
     }else{
        d1 <- dist(t(rld))
     }
     sample_dist <- as.matrix(d1)
     ## coloring the diagonal often skews the color scale
     if(diagNA) diag(sample_dist) <- NA
     ## clrs
     n <- 9
     if(palette %in% c("BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral")) n <- 11
     clrs <- colorRampPalette( brewer.pal(n, palette))(255)
     if(revclrs) clrs <- rev(clrs)

     if(output == "pheatmap"){
        ##  dendsort to reorder branches
        callback = function(hc, ...){dendsort(hc)}
        if(!missing(intgroup)){
              df <- as.data.frame(colData(rld)[, intgroup, drop=FALSE])
              df[,1] <- paste0(df[,1], "  ")  # hack to fix right margin
        }else{
             df <- NA
        }
        pheatmap(sample_dist, clustering_callback = callback, col=clrs,
                 clustering_distance_rows=d1,
                 clustering_distance_cols=d1, annotation_col=df, ...)
     }else{
          ##  dendsort to reorder branches
        dend <- dendsort( as.dendrogram( hclust(d1) ) , isReverse=TRUE)
        d3heatmap(sample_dist, Rowv=dend, Colv=dend, colors = clrs,
            xaxis_font_size = font_size, yaxis_font_size = font_size, ...)
     }
}
