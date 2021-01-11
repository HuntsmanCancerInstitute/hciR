#' Write DESeq results to GSEA rank file
#'
#' Writes gene name (or human homolog) and log2 fold change sorted in
#' descending order to a list of vectors or tab-delimited file.
#' Duplicate gene name with the lowest absolute fold change are removed.
#'
#' @param res a list of DESeq results
#' @param write write a file (default) or return a list of named vectors
#' @param protein_coding only write protein_coding genes
#' @param na_pvalue remove genes with NA p-values (extreme count outliers and
#' low mean normalized counts)
#'
#' @return Tab-delimited file with gene name and log2 fold change
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   write_gsea_rnk(res)
#' }
#' @export

write_gsea_rnk <- function(res, write=TRUE, protein_coding = TRUE, na_pvalue = TRUE, file = NULL){
   # needs list as input
   if(is.data.frame(res) ){
      n <- attr(res, "contrast")
      res <- list(res)
      names(res) <- n
   }
   ## add check for mouse or human
   n <- length(res)
   rnk <- vector("list", n)
   names(rnk) <- names(res)
   for(i in 1:n){
      y <- res[[i]]
      vs <- gsub( "\\.* ", "_", names(res[i]))
      vs <- gsub("_+_", "_", vs, fixed=TRUE)
      ## add txt for GNomEx
      ## name output file
      if(is.null(file)){
          outfile <- paste0( gsub("/", "", vs), ".rnk")
      } else {
          outfile = as.character(file)
      }

      if(protein_coding && "biotype" %in% names(y)){
           y <- dplyr::filter(y, biotype == "protein_coding")
      }
      if("padj" %in% names(y) ){
         nna <- sum(is.na(y$padj))
         if(na_pvalue & nna>0){
           message("Removing ", nna,  " genes with NA p-values")
           y <- dplyr::filter(y, !is.na(padj))
        }
     }
      if("human_homolog" %in% colnames(y)){
          ## include NAs?  extreme outliers and low mean normalized count
          x <- dplyr::filter( y, human_homolog != "") %>%
                dplyr::select(human_homolog, log2FoldChange)
          names(x)[1] <- "gene_name"
          ## split comma-separated lists!
          x <- tidyr::separate_rows(x, gene_name, sep=",")
      }
      else{
          x <- dplyr::filter( y, gene_name != "") %>%
                dplyr::select(gene_name, log2FoldChange)
      }
      ## remove duplicates
      x <- dplyr::arrange(x, gene_name, dplyr::desc( abs(log2FoldChange)) )
      n <- duplicated(x$gene_name)
      if( i == 1) message( "Removing ", sum(n), " duplicate genes")
      x <- x[!n,]
      x <- dplyr::arrange(x, dplyr::desc( log2FoldChange))
      if(write){
          message("Saving ", outfile)
          readr::write_tsv(x, outfile, col_names=FALSE)
      }else{
         # named vector
         y <- x[[2]]
         names(y) <- x[[1]]
         rnk[[i]] <- y
      }
   }
   if(!write){
     #  if(length(rnk)==1)  rnk <- rnk[[1]]
      rnk
   }
}
