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
   ## use as.matrix.tbl to convert to matrix
   counts <- as.matrix(count_tbl)
   ## round for RSEM
   counts <- round( counts, 0)
   mode(counts) <- "integer"
   ## convert to factor ??
   samples <- factor_trts(sample_tbl )
   dds <- DESeq2::DESeqDataSetFromMatrix(counts, samples, ...)
   dds <- DESeq2::DESeq(dds)
   dds
}


## avoid printing in markdown files... Warning ...some variables in design formula are characters, converting to factors
factor_trts <- function(samples){
   for(i in 1:ncol(samples)){
      if(class(samples[[i]])=="character" ){
         if(sum(duplicated(samples[[i]])) > 1 ){
             samples[[i]] <- factor(samples[[i]])
         }
      }
   }
   samples
}
