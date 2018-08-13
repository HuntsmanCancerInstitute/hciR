#' Plot gage_all output
#'
#' Plot heatmap with enrichment scores by contrast and gene set from \code{\link{gage_all}}
#'
#' @param x list from gage_all
#' @param trim trim long set names, default more than 70 characters
#' @param n_sets display contrasts sharing n or more sets for n > 1.  If n = 1,
#' then only plot unique sets.  If missing, then plots all sets, default.
#' @param cluster_cols cluster columns, default FALSE
#' @param \dots other options passed to \code{pheatmap}
#' @author Chris Stubben
#' @examples
#' \dontrun{
#'   library(hciRdata)
#'   x <- gage_all(res, msig_pathways$KEGG)
#'   plot_gage(x)
#'  }

plot_gage <- function(x, trim=70, n_sets, cluster_cols=FALSE, ...){
   y <- dplyr::bind_rows(x, .id = "contrast")
   ## order columns by order in list (or alphabetical)
   y$contrast <- factor(y$contrast, levels= names(x))
   y$name <- ifelse(nchar(y$name) > trim, paste0(substr(y$name, 1, trim-2), "..."), y$name)
   z <- dplyr::select(y, contrast, name, stat.mean) %>%
         tidyr::spread(contrast, stat.mean)
   if(!missing(n_sets)){
     n <- apply(z[, -1], 1, function(x) sum(!is.na(x)))
     if(n_sets ==1){
        z <- filter(z, n == 1)
     }else{
        z <- filter(z, n >= n_sets)
     }
   }
   clrs <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(255)
   ## too many NAs to cluster
   z <- as_matrix(z)
   n1 <- max(abs(z), na.rm=TRUE)
   brks <- seq(-n1, n1, length = 255)
   pheatmap::pheatmap(z, color = clrs, breaks = brks, cluster_cols=cluster_cols, cluster_rows=FALSE, ...)
}
