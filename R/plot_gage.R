#' Plot gage_all output
#'
#' Plot heatmap with enrichment scores by contrast and gene set from \code{\link{gage_all}}
#'
#' @param x list from gage_all
#' @param \dots other options passed to \code{pheatmap}
#' @author Chris Stubben
#' @examples
#' \dontrun{
#'   data(msig)
#'   x <- gage_all(res, msig$KEGG)
#'   plot_gage(x)
#'  }

plot_gage <- function(x, ...){
   y <- dplyr::bind_rows(x, .id = "contrast")
   ## fix KEGG names in all CAPS
   y$name <- stringr::str_to_title(gsub("_", " ", y$name))
   CAPs <- c(Dna= "DNA", Trna="tRNA", Tca="TCA", And="and")
   y$name <- stringr::str_replace_all(y$name, CAPs)
   z <- dplyr::select(y, contrast, name, stat.mean) %>%
         tidyr::spread(contrast, stat.mean)
   clrs <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(255)
   pheatmap::pheatmap(as_matrix(z), color = clrs, cluster_rows=FALSE, cluster_cols=FALSE, ...)
}
