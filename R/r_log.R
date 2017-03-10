#' Regularized log transformation
#'
#' This is a wrapper for \code{rlog} and avoids loading the DESeq2 package for basic workflows.
#'
#' @param dds a DESeqDataSet
#' @param \dots additional options passed to \code{rlog}
#'
#' @return A DESeqTransform object
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rlog_dds(dds)
#' }
#' @export

r_log <- function( dds,  ...){
   rld <- DESeq2::rlog(dds, ...)
   rld
}
