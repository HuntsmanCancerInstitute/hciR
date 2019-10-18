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
#' summary_deseq(pasilla$results)
#' @export


summary_deseq <-  function(object){
   if(class(object)[1] == "DESeqResults") {
        alpha  <- S4Vectors::metadata(object)$alpha
        ft     <- S4Vectors::metadata(object)$filterThreshold
        ## displayed in print.DESeqResults?
        contrast <- ""
        ## if using tibble from results_all or  annotate_results
      }else{
          alpha <- attr(object, "alpha")
          contrast <- attr(object, "contrast")
          ft <- attr(object, "filterThreshold")
      }
   if(!is.null(ft))   ft <- round(ft, 1)
   if(is.null(alpha)) alpha <- 0.05
   notallzero <- sum(object$baseMean > 0)
    outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
   # check if  "ihwResult" %in% names(metadata(object)) ??
    filt <- sum(!is.na(object$pvalue) & is.na(object$padj))


   up <- sum(object$padj < alpha & object$log2FoldChange > 0, na.rm = TRUE)
 down <- sum(object$padj < alpha & object$log2FoldChange < 0, na.rm = TRUE)
  x1 <- tibble::tibble( summary = c("up-regulated", "down-regulated", "outliers", paste("low counts <", ft) ),
                    count = c( up, down, outlier, filt),
                 percent = round(c( up, down, outlier, filt) /notallzero *  100, 2)
                  )
   message(contrast)
  x1
}
