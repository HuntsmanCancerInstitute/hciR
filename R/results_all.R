#' Extract all results from a DESeq analysis
#'
#' Extract all possible contracts and result tables from a DESeq analysis. 
#' Currently supports simple designs with a single variable.
#'
#' @param object a DESeqDataSet
#' @param alpha the significance cutoff for the adjusted p-value cutoff (FDR)
#' @param \dots additional options passed to \code{results}
#'
#' @return A list of DESeqResults objects for each contrast
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   res <- results_all(dds)
#' }
#' @export

results_all <- function( object, alpha = 0.05, ...){
   message("Using adjusted p-value < ", alpha)
   n <- as.character(design(object))
   ## [1] "~"   "condition"
   if(length(n) > 2) stop("The design has multiple variables and only simple designs are currently supported")
   trt <- n[2]

   n <- levels( colData(object)[[trt]])
   contrast <- combn(n, 2)
   res <- vector("list", length(n))
   vs <- apply(contrast, 2, paste, collapse = " vs. ")
   names(res) <- vs
   ## padded for message
   vs1 <- sprintf(paste0("%-", max(nchar(vs))+2, "s"), paste0(vs, ":") )

   for(i in seq_along(n)){
       res1 <- results(object, contrast = c( trt, contrast[1,i], contrast[2,i] ), alpha = alpha, ...)
        x <- summary_DESeq(res1)
        message(i, ". ", vs1[i], x[1,2], " up and ", x[2,2], " down regulated" )
      res[[i]] <- res1
   }
   res
}
