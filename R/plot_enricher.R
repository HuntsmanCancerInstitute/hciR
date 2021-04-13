#' Plot enricher output
#'
#' Plot heatmap with enrichment scores by contrast
#'
#' @param x list from \code{link{fgsea_all}}
#' @param trim trim long names, default more than 70 characters
#' @param sets display contrasts sharing n or more sets for n > 1.  If n = 1,
#' then only plot unique sets.  If missing, then plots all sets, default.
#' @param cluster_row Cluster dendrogram rows, default FALSE for an alphabetical list
#' @param cluster_row Cluster dendrogram columns, default FALSE
#' @param \dots other options passed to \code{pheatmap}
#' @author Chris Stubben
#' @examples
#' \dontrun{
#'   library(hciRdata)
#'   fc <- write_gsea_rnk(res, write=FALSE)
#'   x <- fgsea_all(fc, msig_pathways$KEGG, FDR= 0.25)
#'   plot_fgsea(x)
#'  }
#' @export

plot_enricher <- function(x, trim=70, sets, cluster_row=FALSE,  cluster_col=FALSE, ...){
   if(is.data.frame(x)) stop("A list of from enricher_all is required")
   y <- dplyr::bind_rows(x, .id = "contrast")
   ## order columns by order in list (or alphabetical)
   y$contrast <- factor(y$contrast, levels= names(x))
   y$pathway <- ifelse(nchar(y$pathway) > trim,
                   paste0(substr(y$pathway, 1, trim-2), "..."), y$pathway)
   z <- dplyr::select(y, contrast, pathway, p.adjust) %>%
        dplyr::mutate(p.adjust = -log10(p.adjust)) %>%
         tidyr::spread(contrast, p.adjust)

   if(!missing(sets)){
     n <- apply(z[, -1], 1, function(x) sum(!is.na(x)))
     if(sets ==1){
        z <- filter(z, n == 1)
     }else{
        z <- filter(z, n >= sets)
     }
   }
   clrs <- grDevices::colorRampPalette(
              c("white", RColorBrewer::brewer.pal(n = 9, name = "Reds")[-1]))(255)
   ## too many NAs to cluster
   z <- as_matrix(z)
   z[is.na(z)] <- 0
   message(nrow(z) , " total sets")
   pheatmap::pheatmap(z, color = clrs, cluster_cols=cluster_col,
       cluster_rows=cluster_row, ...)
}
