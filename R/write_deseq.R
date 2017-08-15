#' Write DESeq results
#'
#' Write DESeq result files, raw counts, rlog values, normalized counts to an Excel file.
#'
#' @param result_all a list from \code{results_all}
#' @param dds a DESeqDataset object with count tables
#' @param rld a DESeqTransform obect with rlog values
#' @param biomart annotations from \code{read_biomart}
#' @param txt_files write results to separate txt files, mainly for IPA input
#' @param file file name
#' @param \dots additional options passed to \code{annotate_results}
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  write_deseq(res, dds, rld, hsa)
#' }
#' @export

write_deseq <- function(result_all, dds, rld, biomart, txt_files = FALSE, file = "DESeq.xlsx", ...){

   ##  if results are a tibble (since simplify=TRUE by default)
   if(!class(result_all)[1] == "list"){
         n <- attr(result_all, "contrast")
         result_all <- list(result_all)
         names(result_all) <- n
   }
   if(txt_files){
      res <- result_all
      for (i in 1:length(res)){
         vs <- gsub( "\\.* ", "_", names(res[i]))
         vs <- gsub("_+_", "_", vs, fixed=TRUE)
         vs <- paste0(vs, ".txt")
         message( "Saving ",  vs)
         readr::write_tsv(res[[i]], vs)
      }
   }else{
     if( !class(dds)[1] == "DESeqDataSet") stop("dds should be a DESeqDataSet object")
     if( !class(rld)[1] == "DESeqTransform") stop("rld should be a DESeqTransform object")

   ## 1. summary
   sum1 <-  dplyr::bind_rows(lapply(result_all, summary_deseq), .id= "contrast")
   # write.xlsx does not like tibbles, so use as.data.frame
   sum1 <- as.data.frame( sum1)
   # write.xlsx requires a named list for writing mulitple worksheets
   sum1 <- list("summary" = sum1)
   res1 <- lapply(result_all, as.data.frame )

   ## write.xlsx replaces space with .
   names(res1)  <- gsub( "\\.? ", "_", names(res1))
   ## forward slash will cause Excel errors
   names(res1)  <- gsub( "/", "", names(res1))

   # sample data in colData ... drop replaceable
   samp1 <- as.data.frame(SummarizedExperiment::colData(dds))
   # treatments is a list is TCGA biolinks
    samp1 <- samp1[, colnames(samp1) != "treatments"]
    samp1$replaceable <- NULL

   DESeq_tables <-  c(
     sum1,
     res1,
   list(
     "raw_counts" = DESeq2::counts(dds),
     "normalized" = DESeq2::counts(dds, normalized=TRUE),
     "rlog"       = SummarizedExperiment::assay(rld),
     "samples"    = samp1,
     "Ensembl"    = as.data.frame(biomart))
  )
   message("Saving ", length(DESeq_tables), " worksheets to ", file)
  # DESeq_tables
   openxlsx::write.xlsx(DESeq_tables, file = file, rowNames= sapply(DESeq_tables, is.matrix) )
  }
}
