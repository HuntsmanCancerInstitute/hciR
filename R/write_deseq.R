#' Write DESeq results
#'
#' Write DESeq result files, raw counts, rlog values, normalized counts to an Excel file.
#'
#' @param result_all a list from \code{results_all}
#' @param dds a DESeqDataset object with count tables
#' @param rld a DESeqTransform obect with rlog values
#' @param biomart annotations from \code{read_biomart}
#' @param sets table with set intersections, optional
#' @param fpkm matrix of fpkm values, optional
#' @param text_files write results to separate txt files, mainly for IPA input
#' @param cutoff if text_files is TRUE, only write genes below this p-value cutoff
#' @param file file name
#' @param \dots additional options passed to \code{annotate_results}
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  write_deseq(res, dds, rld, human98)
#' }
#' @export

write_deseq <- function(result_all, dds, rld, biomart, sets, fpkm,
   text_files = FALSE, cutoff, file = "DESeq.xlsx", ...){

   ##  if results are a tibble (since simplify=TRUE by default)
   if(!class(result_all)[1] == "list"){
         n <- attr(result_all, "contrast")
         result_all <- list(result_all)
         names(result_all) <- n
   }
   if(text_files){
      res <- result_all
      names(res)  <- gsub( "/", "", names(res))

      for (i in 1:length(res)){
         vs <- gsub( "\\.* ", "_", names(res[i]))
         vs <- gsub("_+_", "_", vs, fixed=TRUE)
         vs <- paste0(vs, ".txt")
         message( "Saving ",  vs)
         if(!missing(cutoff))  res[[i]] <- filter(res[[i]], padj < cutoff)
         readr::write_tsv(res[[i]], vs)
      }
   }else{
     if( !class(dds)[1] == "DESeqDataSet") stop("dds should be a DESeqDataSet object")
     if( !class(rld)[1] == "DESeqTransform") stop("rld should be a DESeqTransform object")

   ## 1. summary  - will print contrast name (which may change if longer than 32 characters)
   sum1 <- suppressMessages( dplyr::bind_rows(lapply(result_all, summary_deseq), .id= "contrast"))
   # write.xlsx does not like tibbles, so use as.data.frame
   sum1 <- as.data.frame( sum1)
   # write.xlsx requires a named list for writing mulitple worksheets
   sum1 <- list("summary" = sum1)
   res1 <- lapply(result_all, as.data.frame )

   ## write.xlsx replaces space with .
   names(res1)  <- gsub( "\\.? ", "_", names(res1))
   ## forward slash will cause Excel errors
   names(res1)  <- gsub( "/", "", names(res1))
   # also not allowed  \ * [ ] : ?

   # Worksheet names cannot exceed 31 characters.
   if(any(nchar(names(res1)) > 31)){
	  message("Note: some contrast names are longer than the 31 characters and will be truncated")
      z <- substr(names(res1), 1, 31)
	  if(any(duplicated(z))) stop("Warning: names(result_all) is not unique after truncating. Please shorten the contrast names for Excel")
	  names(res1) <- z
   }
   for(i in names(res1)) message(i)

   # sample data in colData ... drop replaceable
   samp1 <- as.data.frame(SummarizedExperiment::colData(dds))
   # TCGAbiolinks has list columns!
    samp1 <- samp1[, colnames(samp1) != "treatments"]
    n <- which(sapply(samp1, class) == "AsIs")
    for (i in n) samp1[[i]] <- unlist(samp1[[i]])

    samp1$replaceable <- NULL
  Counts <-    list(
       "raw_counts" = DESeq2::counts(dds),
       "log2_norm"  = SummarizedExperiment::assay(DESeq2::normTransform(dds)),
       "rlog"       = SummarizedExperiment::assay(rld))

     if(!missing(fpkm)){
       if(class(fpkm)[1] == "list"){
          Counts <- c( Counts, fpkm)
       }else{
          Counts <- c( Counts, list("fpkm"= fpkm))
      }
   }
   Meta <- list( "samples" = samp1)
   if(!missing(biomart)){
        Meta <- c(Meta, list( "ensembl"= as.data.frame(biomart)))
        # add gene names to count tables?
        Counts <- lapply(Counts, function(y) dplyr::right_join(biomart[,1:3],
                tibble::rownames_to_column( as.data.frame(y), "id"), by="id"))
   }
   DESeq_tables <-  c( sum1, res1, Counts, Meta)
   if(!missing(sets)) DESeq_tables <-  c( sum1, res1, list(Sets=sets), Counts, Meta)

   message("Saving ", length(DESeq_tables), " worksheets to ", file)
  # DESeq_tables
   openxlsx::write.xlsx(DESeq_tables, file = file, rowNames= sapply(DESeq_tables, is.matrix) )
  }
}
