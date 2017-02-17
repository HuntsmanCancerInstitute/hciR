#' Convert tbl_df to a matrix
#'
#' Mainly used to convert count tables to matrix for \code{rlog} and other
#' functions that require a matrix as input.
#'
#' @param x a tbl_df with identifiers in column 1 (assigned as rownames) and
#'  numbers in remaining columns
#' @param \dots additional options for \code{as.matrix}
#' @return A matrix
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'   rlog(as.matrix(fc_tbl))
#' }
#' @export


as.matrix.tbl_df <- function(x, ...){
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}
