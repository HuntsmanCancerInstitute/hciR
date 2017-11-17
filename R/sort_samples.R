#' Sort sample names
#'
#' Sort character vectors by removing non-numeric characters.  This works for very simple names,
#'  try \code{mixedsort} in \code{ggtools} for more complicated examples
#'
#' @param samples a vector of sample names
#'
#' @return A vector
#'
#' @author Chris Stubben
#'
#' @examples
#'  ids <- c("135X1", "135X21", "135X2", "135X10")
#'  sort(ids)
#'  sort_samples(ids)
#' @export

sort_samples <- function(samples){
   samples[order_samples(samples)]
}

#' @describeIn sort_samples Order sample names
#' @export
order_samples <- function( samples ){
    order(as.numeric(gsub("[^0-9]", "",  samples )))
}
