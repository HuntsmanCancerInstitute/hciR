#' Read a DESeq Excel file from write_deseq
#'
#' @param file write_deseq output file
#' @param object type of object to create, either results, rlog or log2_norm
#'
#' @return a list of tibbles or DESeqTransform object
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  res <- read_deseq("DESeq.xlsx")
#'  rld <- read_deseq("DESeq.xlsx", "rlog")
#' }
#' @export

read_deseq <- function(file,  object="results"){
   if(object =="results"){
     ## load all worksheets with _vs_ ???
    wk <- readxl::excel_sheets(file)
    n <- grep("_vs_", wk)
    if(length(n)==0) stop("No results with _vs_ in worksheet names")
      res <- vector("list", length(n))
      names(res) <- gsub("_vs_", " vs. ", wk[n])
      for(i in seq_along(n))  {
          message("Loading ", names(res)[i])
          res[[i]] <- read_excel(file, sheet= n[i])
      }
   } else{
      s1 <- readxl::read_excel(file, sheet="samples")
      message("Loading ", object, " worksheet")
      r1 <- read_excel(file, sheet= object)
      ## delete gene_name and biotype
      n2 <- which(colnames(r1) %in% c("gene_name", "biotype"))
      if(length(n2) !=2) message("Note: Gene name and biotype columns are missing?")
      r1 <- as_matrix(r1[, -n2])
      s1 <- SummarizedExperiment::SummarizedExperiment(r1, colData = s1)
      res <- DESeq2::DESeqTransform(s1)
   }
   res
}
