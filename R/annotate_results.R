#' Annotate DESeq2 results
#'
#' Add BioMart annotations to a DESeq2 result table
#'
#' @param result a DESeqResults object with Ensembl IDs as row names
#' @param biomart annotations from \code{read_biomart} with id column matching
#' row names in results
#' @param add a vector of biomart columns to add to result table, default
#' gene_name, biotype and description
#' @param id position or name of column in biomart that matches rows names in
#' the DESeq results
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  annotate_results(res, human98)
#' }
#' @export

annotate_results <- function(result, biomart, add, id = 1 ){
   n1 <- match( rownames(result), biomart[[id]])
   if(missing(add)){
        add <- c("gene_name", "biotype", "description")
        if("human_homolog" %in% names(biomart)) add <- c(add, "human_homolog")
    }
   if(all(is.na(n1))) stop("Rownames in results do not match column ", id, " in biomart table")
   if(any(is.na(n1))) message(sum(is.na(n1)), " rows in results are missing from biomart table")
   #use data.frame to avoid tibble error ... Each variable must be a 1d atomic vector or list.
   res1 <- data.frame(id=rownames(result), biomart[n1,  add], result, stringsAsFactors=FALSE)
   tibble::as_tibble(res1)
}
