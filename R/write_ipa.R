#' Write DESeq results to text files for bulk uploading into IPA
#'
#' @param result_all a list from \code{results_all}
#' @param run add run ID to file name
#' @param filename a name for the output file, string, e.g., "ipa_formatted_results.txt"
#' 
#' @return Text files
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  write_ipa(res)
#' }
#' @export

write_ipa <- function(result_all, run, filename = NULL){

   ##  if results are a tibble (since simplify=TRUE by default)
   if(!class(result_all)[1] == "list"){
         n <- attr(result_all, "contrast")
         result_all <- list(result_all)
         names(result_all) <- n
   }
   res <- result_all
   names(res)  <- gsub( "/", "", names(res))
if(!missing(run)) names(res)  <- paste(run, names(res))
   
   for (i in 1:length(res)){
      # name output file
      if(filename == NULL){
        vs <- gsub( "\\.* ", "_", names(res[i]))
        vs <- gsub("_+_", "_", vs, fixed=TRUE)
        vs <- paste0(vs, ".txt")
      } else {
          vs = filename
      }
      message( "Saving ",  vs)
      z <- dplyr::select(res[[i]], id, baseMean, log2FoldChange, padj)
      readr::write_tsv(z, vs)
   }
}
