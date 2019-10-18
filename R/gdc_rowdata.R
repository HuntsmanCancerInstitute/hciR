#' GDC row data
#'
#' Joins the \code{rowData} slot from \code{GDCprepare} to Ensembl gene annotations
#'
#' @param gdc a SummarizedExperiment object from \code{GDCprepare} in the \code{TCGAbiolinks} package.
#' @param ensembl Ensembl human gene annotations, version 87 is best.
#'
#' @note Currently works with harmonized data with Ensembl gene ids.
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' brca <- GDCprepare(query)
#' library(hciRdata)
#' genes <- gdc_rowdata(brca, human90)
#' genes
#' }
#' @export


gdc_rowdata <- function(gdc, ensembl){
   genes <- tibble::as_tibble(as.data.frame(SummarizedExperiment::rowData(gdc)))
   if(names(genes)[1] != "ensembl_gene_id"){
      message("Missing ensembl_gene_id in Summarized Experiment")
   }else{
       # drop original_ensembl_gene_id
      genes <- dplyr::left_join(genes[, 1:2], ensembl, by=c(ensembl_gene_id="id"))
   }
   genes
}
