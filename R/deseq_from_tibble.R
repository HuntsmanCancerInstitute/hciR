#' Run DESeq from tibbles
#'
#' Creates DESeqDataSet object using count and sample tibbles as input and optionally runs DESeq
#'
#' @param counts a count table with feature ids in row 1
#' @param samples a sample table with names matching count table column names in row 1
#' @param run_deseq Run \code{\link{DESeq}}, default TRUE
#' @param \dots additional options like the design formula are passed to \code{DESeqDataSetFromMatrix}
#' @return A DESeqDataSet object
#'
#' @note This function first runs \code{\link{sort_counts}} to check and
#' reorder count columns by the first column in samples, and then \code{DESeqDataSetFromMatrix}
#' and optionally \code{DESeq}.  To match DESeq2 versions < 1.16, set run_deseq = FALSE and then
#' run \code{DESeq(dds, betaPrior =TRUE)}'
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    deseq_from_tibble(counts, samples, design = ~ trt)
#' }
#' @export

deseq_from_tibble <- function( counts, samples, run_deseq = TRUE, ...){
   counts <- sort_counts(counts, samples)
   ## use hciR::as_matrix to convert to matrix
   counts <- as_matrix(counts)
   ## round for RSEM?
   if (any(round(counts[1:10,]) != counts[1:10,])) {
     message("Rounding counts (from RSEM?)")
     counts <- round( counts, 0)
   }
   mode(counts) <- "integer"
   dds <- DESeq2::DESeqDataSetFromMatrix(counts, samples, ...)
   if(run_deseq) dds <- DESeq2::DESeq(dds)
   dds
}
