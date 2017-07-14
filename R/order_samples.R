#' Order sample names
#'
#' Order character vectors by removing non-numeric characters.  
#'
#' @param samples a vector of sample names
#'
#' @return A vector
#'
#' @author Chris Stubben
#'
#' @examples
#'  ids <- c("135X1", "135X21", "135X2", "135X10")
#'  order(ids)
#'  order_samples(ids)
#' @export

order_samples <- function( samples ){
    order(as.numeric(gsub("[^0-9]", "",  samples )))
}
