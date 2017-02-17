#' Summarize DESeq2 results
#'
#' Summarize DESeq2 results into a \code{data.frame}
#'
#' @param object a DESeqResults object
#'
#' @return A data.frame
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' summary_deseq(res1)
#' }
#' @export


summary_deseq <-  function(object){
   if(class(object)[1] == "DESeqResults") {
        alpha      <- S4Vectors::metadata(object)$alpha
        ## if using tibble from results_all or  annotate_results
      }else{
          alpha <- attr(object, "alpha")
      }
   notallzero <- sum(object$baseMean > 0)
    outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
   # check if  "ihwResult" %in% names(metadata(object)) ??
    filt <- sum(!is.na(object$pvalue) & is.na(object$padj))

   up <- sum(object$padj < alpha & object$log2FoldChange > 0, na.rm = TRUE)
 down <- sum(object$padj < alpha & object$log2FoldChange < 0, na.rm = TRUE)
  x1 <- data.frame( summary = c("up-regulated", "down-regulated", "outliers", "low counts"),
                    count = c( up, down, outlier, filt),
                 percent = round(c( up, down, outlier, filt) /notallzero *  100, 2)
                  )
  x1
}
