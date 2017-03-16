#' Convert tibbles to a matrix
#'
#' Mainly used to convert count tables to matrix for Bioconductor functions that
#' require a matrix as input.
#'
#' @param x a tibble with identifiers in column 1 (assigned as rownames) and
#'  numbers in remaining columns
#' @return A matrix
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rlog(as_matrix(count_tbl))
#' }
#' @export


as_matrix <- function(x, ...){
   if(!tibble::is_tibble(x) ) stop("x must be a tibble")   
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}
