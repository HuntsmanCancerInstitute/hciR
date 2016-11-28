#' Sort sample names
#'
#' Sort character vectors by removing non-numeric characters.  This works for very simple cases,
#'  try \code{mixedsort} in \code{ggtools} for more complicated examples
#'
#' @param samples a vector of sample names
#'
#' @return A vector
#'
#' @author Chris Stubben
#'
#' @examples
#'  samples <- c("135X1", "135X2", "135X11")
#'  sort(samples)
#'  sample_sort(samples)
#' @export

sample_sort <- function(samples){
   samples[sample_order(samples)]
}

sample_order <- function( samples ){
    order(as.numeric(gsub("[^0-9]", "",  samples )))
}
