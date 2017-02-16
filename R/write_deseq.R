#' Write DESeq results to an Excel file
#'
#' Write DESeq result files, raw counts, rlog values, normalized counts.
#'
#' @param dds a DESeqDataset objects object with Ensembl IDs as row names
#' @param result_all a list from \code{results_all}
#' @param rld a DESeqTransform obect with rlog values
#' @param biomart annotations from \code{read_biomart}
#' @param file file name
#' @param \dots additional options passed to \code{annotate_results}
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  write_deseq(dds, res, rld, hg)
#' }
#' @export

write_deseq <- function(dds, result_all, rld, biomart, file = "DESeq.xlsx", ...){

   ##  checks
   if(!class(result_all) == "list") stop("result_all should be a list of DESeqResults")

   ## 1. summary
   sum1 <-  bind_rows(lapply(result_all, summary_deseq), .id= "contrast")
   # write.xlsx does not like tibbles, so use as.data.frame
   sum1 <- as.data.frame( sum1)
   # write.xlsx requires a named list for writing mulitple worksheets
   sum1 <- list("summary" = sum1)

    res1 <- lapply(result_all, as.data.frame )

   ## write.xlsx replaces space with .
   names(res1)  <- gsub( "\\.? ", "_", names(res1))

   # sample data in colData ... drop replaceable
   samp1 <- as.data.frame(colData(dds))
   samp1$replaceable <- NULL

   DESeq_tables <-  c(
     sum1,
     res1,
   list(
     "raw_counts" = counts(dds),
     "normalized" = counts(dds, normalized=TRUE),
     "rlog"       = assay(rld),
     "samples"    = samp1,
     "Ensembl"    = as.data.frame(biomart))
  )
   message("Saving ", length(DESeq_tables), " worksheets to ", file)
 #DESeq_tables
   write.xlsx(DESeq_tables, file = file, rowNames=c(FALSE, rep(FALSE, length(res1)),TRUE,TRUE,TRUE, FALSE,FALSE))
}
