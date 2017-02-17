#' Regularized log transformation
#'
#' This is a wrapper for \code{rlog} and avoids loading the DESeq2 package.
#'
#' @param a DESeqDataSet
#' @param \dots additional options passed to \code{rlog}
#'
#' @return A DESeqTransform object
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rlog_dds( dds)
#' }
#' @export

rlog_dds <- function( dds,  ...){
   rld <- DESeq2::rlog(dds, ...)
   rld
}
