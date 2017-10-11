#' Write DESeq results to GSEA rank file
#'
#' Writes gene name (or human homolog) and log2 fold change sorted in
#'  descending order to tab-delimited file.  Duplicate gene name with the lowest
#'  absolute fold change value are removed.
#'
#' @param res a list of DESeq results
#'
#' @return Tab-delimited file with gene name and log2 fold change
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   write_gsea_rnk(res)
#' }
#' @export

write_gsea_rnk <- function(res){
   # needs list as input
   if(is.data.frame( res )){
      n <- attr(res, "contrast")
      res <- list(res)
      names(res) <- n
   }
   ## add check for mouse or human
   n <- length(res)
   for(i in 1:n){
        y <- res[[i]]
        vs <- gsub( "\\.* ", "_", names(res[i]))
        vs <- gsub("_+_", "_", vs, fixed=TRUE)
         ## add txt for GNomEx
      outfile <- paste0( gsub("/", "", vs), ".rnk")
      if("human_homolog" %in% colnames(y)){
          x <- dplyr::filter( y, human_homolog != "") %>%
                dplyr::select(human_homolog, log2FoldChange)
         names(x)[1] <- "gene_name"
      }else{
          x <- dplyr::filter( y, gene_name != "") %>%
                dplyr::select(gene_name, log2FoldChange)
      }
     ## remove duplicates
     x <- dplyr::arrange(x, gene_name, desc( abs(log2FoldChange)))
      n <- duplicated(x$gene_name)
     if( i == 1) message( "Removing ", sum(n), " duplicate genes")
      x <- x[!n,]
      x <-  dplyr::arrange(x, desc( log2FoldChange))
      message("Saving ", outfile)
      readr::write_tsv(x, outfile, col_names=FALSE)
   }
}
