#' Read GSEA enrichment files
#'
#' Reads both positive and negative enrichment results from Broad's GSEA into a single table
#'
#' @param gsea_dir Path to GSEA output files
#' @param FDR FDR cutoff, default 0.1
#' @param gsea_file GSEA output file
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' x <- read_gsea("~/gsea_home/output/dec15/my_analysis.GseaPreranked.1513361740375")
#' }
#'
#' @export

read_gsea <-function( gsea_dir, FDR = 0.1){
    # timestamp should be part of directory name..
   timestamp <- gsub(".*?([0-9]+)$", "\\1", gsea_dir)
   na_pos <- paste0("gsea_report_for_na_pos_", timestamp, ".xls")
   gsea1 <- paste0( gsea_dir, "/", na_pos)
   if(!file.exists(gsea1) ) stop("No file found matching ", na_pos)
   x1 <- read_gsea_file(gsea1)
   x2 <- read_gsea_file(gsub("na_pos", "na_neg", gsea1))
   x <- dplyr::bind_rows( list( positive = x1, negative = x2 ), .id="enriched") %>%
     dplyr::filter(fdr < FDR)
   if(nrow(x) == 0){
       message( "No enriched sets found")
   }else{
   message(  nrow(x), " enriched sets (",
      sum(x$enriched=="positive"), " positive, ", sum(x$enriched=="negative"), " negative)")
   }
   x
}

#' @describeIn read_gsea Read a single GSEA output file
#' @export

read_gsea_file <- function(gsea_file){
   ## suppress Missing column names filled in: 'X12' [12]
   x <- suppressWarnings( readr::read_tsv(gsea_file, col_types="c--idddddi--"))
   names(x) <- c("name", "size", "es", "nes", "p", "fdr", "fwer", "rank")
   # split set from name: KEGG_,  REACTOME_ , etc
   x <- tidyr::separate(x, name, into=c("set", "name"), sep="_", extra="merge")
   x
}
