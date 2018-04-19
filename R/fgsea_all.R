#' Run fgsea on all DESeq2 result tables
#'
#' A wrapper for fgsea to convert ouput into tibbles
#'
#' @param fc A list of ranked fold changes from \code{write_gsea_rnk}
#' @param gsets gene set
#' @param FDR FDR cutoff, default 0.1
#' @param nperm Number of permutations, default 10000
#' @param \dots Additional options passed to \code{fgsea}
#'
#' @return A list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' data(msig)
#' fc <- write_gsea_rnk(res, write=FALSE)
#' x <- fgsea_all(fc, gsets=msig$KEGG)
#' }
#'
#' @export

fgsea_all <- function(fc, gsets, FDR = 0.1, nperm=10000, ...){
   n <- length(fc)
   fgsea_res <- vector("list", n)
   vs <- names(fc)
   names(fgsea_res) <-  vs
   ## padded for message
   vs1 <- sprintf(paste0("%-", max(nchar(vs))+2, "s"), paste0(vs, ":") )
   message( "Using FDR < ", FDR)
   for(i in 1:n){
      f1 <- fgsea::fgsea(pathways = msig$KEGG,
                        stats = fc[[i]],
                        nperm=nperm, ...)
      f1 <- tibble::as_tibble(f1)
      f1 <- dplyr::filter(f1, padj < FDR)
      if(nrow(f1) == 0){
          message( i, ". ", vs1[[i]], " No enriched sets found")
      }else{
         f1$enriched <- "negative"
         f1$enriched[f1$ES>0] <- "positive"
         f1 <- dplyr::arrange(f1, enriched, padj, pval)
      message( i, ". ", vs1[[i]], " ", nrow(f1), " enriched sets (",
         sum(f1$enriched=="positive"), " positive, ", sum(f1$enriched=="negative"), " negative)")
      }
      fgsea_res[[i]] <-  f1
   }
   if(length( fgsea_res) == 1 ) fgsea_res <- fgsea_res[[1]]
   fgsea_res
}
