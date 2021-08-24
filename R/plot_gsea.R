#' Plot GSEA output
#'
#' Plot heatmap with enrichment scores by contrast
#'
#' @param x list with two or more \code{link{read_gsea}} tables
#' @param trim trim long set names, default more than 70 characters
#' @param n_sets display contrasts sharing n or more sets for n > 1.  If n = 1,
#' then only plot unique sets.  If missing, then plots all sets, default.
#' @param nes plot NES (or ES if FALSE)
#' @param \dots other options passed to \code{pheatmap}
#' @author Chris Stubben
#' @examples
#' \dontrun{
#'   data(msig)
#'   x <- list( X1= read_gsea("gsea1.txt"), X2 = read_gsea("gsea2.txt") )
#'   plot_gsea(x)
#'  }
#' @export

plot_gsea <- function(x, trim=70, n_sets, nes=TRUE, ...){
   if(is.data.frame(x)) stop("A list of read_gsea tables is required")
   y <- dplyr::bind_rows(x, .id = "contrast")
   ## order columns by order in list (or alphabetical)
   y$contrast <- factor(y$contrast, levels= names(x))
   y$name <- ifelse(nchar(y$name) > trim,
                 paste0(substr(y$name, 1, trim-2), "..."), y$name)
   ## ES or NES ?
   if(nes){
   z <- dplyr::select(y, contrast, name, nes) %>%
         tidyr::spread(contrast, nes)
   }else{
   z <- dplyr::select(y, contrast, name, es) %>%
         tidyr::spread(contrast, es)
   }
   if(!missing(n_sets)){
     n <- apply(z[, -1], 1, function(x) sum(!is.na(x)))
     if(n_sets ==1){
        z <- dplyr::filter(z, n == 1)
     }else{
        z <- dplyr::filter(z, n >= n_sets)
     }
   }
   clrs <- grDevices::colorRampPalette(rev(
               RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(255)
   ## too many NAs to cluster
   z <- as_matrix(z)
   message(nrow(z) , " total sets")
   n1 <- max(abs(z), na.rm=TRUE)
   brks <- seq(-n1, n1, length = 255)
   pheatmap::pheatmap(z, color = clrs, breaks = brks, cluster_cols=FALSE,
      cluster_rows=FALSE, ...)
}
