#' Extract and annotate a single DESeq results table
#'
#' Extract and annotate a single contrast.  See \code{\link{results_all}} to annotate all contrasts.
#'
#' @param dds a DESeqDataSet
#' @param biomart annotations from \code{read_biomart} with column 1 matching row names in results
#' @param exp name of experimental trt group in the numerator
#' @param ref name of reference trt group in the denominator
#' @param trt treatment group in the design formula
#' @param alpha the significance cutoff for the adjusted p-value cutoff (FDR)
#' @param lfcShrink  shrink fold changes using \code{lfcShrink} for DESeq2 version >= 1.16
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(hciRdata)
#' res <- results2(pasilla$dds, fly98, "treated", "untreated", trt="condition")
#' }
#' @export

results2 <- function(dds, biomart, exp, ref, trt="trt", alpha =0.05, lfcShrink=TRUE){
      vs <- paste(c(exp, ref), collapse= " vs. ")
      res1 <- DESeq2::results(dds, contrast = c(trt, exp, ref), alpha= alpha)
      if(lfcShrink) res1 <- DESeq2::lfcShrink(dds, contrast = c(trt, exp, ref), res = res1)
      ft <- S4Vectors::metadata(res1)$filterThreshold
      Columns <- c("gene_name", "biotype", "chromosome",  "description")
      if("human_homolog" %in% names(biomart)) Columns <- c(Columns, "human_homolog")
      res1 <- annotate_results(res1, biomart, Columns)
       x <- suppressMessages(summary_deseq(res1))
      message(vs, ": ", x[1, 2], " up and ", x[2, 2],  " down regulated")
      attr(res1, "contrast") <- vs
      attr(res1, "alpha") <- alpha
      attr(res1, "filterThreshold") <- ft
      res1
    }
