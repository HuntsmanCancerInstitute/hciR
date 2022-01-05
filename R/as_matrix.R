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
#' x <- top_counts(pasilla$results, pasilla$rlog, top = 25)
#' x
#' as_matrix(x)
#' @export


as_matrix <- function(x) {
  if (!tibble::is_tibble(x)) stop("x must be a tibble")
  y <- as.matrix.data.frame(x[, -1])
  rownames(y) <- x[[1]]
  y
}

#' @describeIn as_matrix Convert matrix to tibble
#' @param var name of column, default id
#' @export

as_tbl <- function(x, var = "id") {
  if (!is.matrix(x)) stop("x must be a matrix")
  tibble::as_tibble(tibble::rownames_to_column(data.frame(x, check.names = FALSE), var))
}
