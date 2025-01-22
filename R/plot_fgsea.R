#' Plot fGSEA output
#'
#' Plot heatmap with enrichment scores by contrast
#'
#' @param x list from \code{link{fgsea_all}}
#' @param trim trim long names, default more than 70 characters
#' @param sets display contrasts sharing n or more sets for n > 1.  If n = 1,
#' then only plot unique sets.  If missing, then plots all sets, default.
#' @param nes plot NES (or ES if FALSE)
#' @param cluster_row Cluster dendrogram rows, default is an alphabetical list
#' @param cluster_col Cluster dendrogram columns, default FALSE
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

plot_fgsea <- function(x, trim=70, sets, nes=TRUE, cluster_row=FALSE,  cluster_col=FALSE, ...){
   if(is.data.frame(x)) stop("A list from fgsea_all is required")
   y <- dplyr::bind_rows(x, .id = "contrast")
   ## order columns by order in list (or alphabetical)
   y$contrast <- factor(y$contrast, levels= names(x))
   y$pathway <- ifelse(nchar(y$pathway) > trim,
                   paste0(substr(y$pathway, 1, trim-2), "..."), y$pathway)
   ## ES or NES ?
   if(nes){
   z <- dplyr::select(y, contrast, pathway, NES) %>%
         tidyr::spread(contrast, NES)
   }else{
   z <- dplyr::select(y, contrast, pathway, ES) %>% 
         tidyr::spread(contrast, ES)
   }
   if(!missing(sets)){
     n <- apply(z[, -1], 1, function(x) sum(!is.na(x)))
     if(sets ==1){
        z <- dplyr::filter(z, n == 1)
     }else{
        z <- dplyr::filter(z, n >= sets)
     }
   }
   clrs <- grDevices::colorRampPalette(
               rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(255)
   ## too many NAs to cluster
   z <- as_matrix(z)
   z[is.na(z)] <- 0
   message(nrow(z) , " total sets")
   n1 <- max(abs(z), na.rm=TRUE)
   brks <- seq(-n1, n1, length = 255)
   pheatmap::pheatmap(z, color = clrs, breaks = brks, cluster_cols=cluster_col,
       cluster_rows=cluster_row, ...)
}
