#' Convert GAGE output to tibbles
#'
#' Converts gage ouput (list of matrices) into list of tibbles
#' and parses matrix rownames into id and name columns.
#'
#' @param gage_results Results from \code{gage}
#'
#' @return A list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' x <- gage(fc, gsets=kegg.sets.mm)
#' as_tibble_gage(x)
#' }
#'
#' @export

as_tibble_gage <- function(gage_results){
   lapply(gage_results, function(x1){
      tbl_df(data.frame(
         id  = gsub("([^ ]+) .*", "\\1", rownames(x1)),
         name= gsub("[^ ]+ (.*)", "\\1", rownames(x1)),
         x1, stringsAsFactors=FALSE))
    })
}
