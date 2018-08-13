#' Run GAGE on all DESeq2 result tables
#'
#' A wrapper for gage to convert ouput into tibbles
#'
#' @param res DESeq results
#' @param gsets gene set
#' @param FDR FDR cutoff, default 0.05
#' @param \dots Additional options passed to \code{gage}
#'
#' @return A list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(hciRdata)
#' x <- gage_all(res, gsets=msig_pathways$KEGG)
#' }
#'
#' @importFrom dplyr filter
#' @export

gage_all <- function(res, gsets, FDR = 0.1, ...){
   ##  if results are a tibble (since simplify=TRUE by default)
   if( class(res)[1] != "list"){
         n <- attr(res, "contrast")
         res <- list(res)
         names(res) <- n
   }
   n <- length(res)
   gage_res <- vector("list", n)
   vs <- names(res)
   names(gage_res) <-  vs
   ## padded for message
   vs1 <- sprintf(paste0("%-", max(nchar(vs))+2, "s"), paste0(vs, ":") )
   message( "Using q-value (FDR) < ", FDR)

   for(i in 1:n){
       logfc <- res[[i]]$log2FoldChange
       names(logfc) <- res[[i]]$gene_name
       if( "human_homolog" %in% colnames(res[[i]]))  names(logfc) <- res[[i]]$human_homolog
       x <- gage::gage(logfc, gsets, ...)
       y <- lapply(x, function(x1){
           tibble::as_tibble(data.frame(
           name = rownames(x1),
           x1, stringsAsFactors=FALSE))  })
      y$stats <- NULL
      x <- dplyr::bind_rows(y, .id="enriched")
      x <- dplyr::filter(x, q.val <  FDR)
      if(nrow(x) == 0){
          message( i, ". ", vs1[[i]], " No enriched sets found")
      }else{
         x$enriched[x$enriched == "less"] <- "negative"
         x$enriched[x$enriched == "greater"] <- "positive"

      message( i, ". ", vs1[[i]], " ", nrow(x), " enriched sets (",
         sum(x$enriched=="positive"), " positive, ", sum(x$enriched=="negative"), " negative)")
      }
      gage_res[[i]] <-  x
   }
   if(length( gage_res) == 1 ) gage_res <- gage_res[[1]]
   gage_res
}
