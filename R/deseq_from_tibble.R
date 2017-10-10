#' Run DESeq from tibbles
#'
#' Creates DESeqDataSet object using count and sample tibbles as input and run DESeq

#' @param count_tbl a count table with feature ids in row 1
#' @param sample_tbl a sample table with names matching count table column names in row 1
#' @param \dots additional options like design formula passed to \code{DESeqDataSetFromMatrix}
#'
#' @return A DESeqDataSet object
#'
#' @note This function first runs \code{\link{sort_counts}} to check and
#' reorder count_tbl columns by the first column in sample_tbl, and then \code{DESeqDataSetFromMatrix}
#' and \code{DESeq}
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    deseq_from_tibble(count_tbl, sample_tbl, design = ~ trt)
#' }
#' @export

deseq_from_tibble <- function( count_tbl, sample_tbl,  ...){
   count_tbl <- sort_counts(count_tbl, sample_tbl)
   ## use as_matrix to convert to matrix
   counts <- as_matrix(count_tbl)
   ## round for RSEM?
   counts <- round( counts, 0)
   mode(counts) <- "integer"
   ## convert to factor ??
   samples <- factor_trts(sample_tbl )
   dds <- DESeq2::DESeqDataSetFromMatrix(counts, samples, ...)
   dds <- DESeq2::DESeq(dds)
   dds
}

#' @describeIn deseq_from_tibble Avoid warnings about factors
#' @param samples a sample table
factor_trts <- function(samples){
   for(i in 1:ncol(samples)){
      if(class(samples[[i]])=="character" ){
         if(sum(duplicated(samples[[i]])) > 0 ){
             samples[[i]] <- factor(samples[[i]])
         }
      }
   }
   samples
}
