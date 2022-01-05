#' Sort columns of a count matrix by rows in sample data
#'
#' \code{DESeqDataSetFromMatrix} requires count columns sorted by the rows in sample data and only checks
#' if sample row names match count column names (and tibbles do not have row names).  This matches column
#' names in a count table to the first column in sample data and reorders columns to match.
#'
#' @param count_tbl count matrix, data.frame or tibble
#' @param sample_tbl sample tabl with a column containing the columns names in count matrix
#' @param id position of column in sample data that matches column names in the
#'   count table, default is the first column
#'
#' @return A count matrix with columns re-ordered to match the sample table
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' count_tbl <- sort_counts(counts, samples)
#' }
#' @export

sort_counts <- function(count_tbl, sample_tbl, id = 1) {
  cname <- deparse(substitute(count_tbl))
  orig_count_tbl <- count_tbl
  if (class(count_tbl)[1] != "matrix") {
    count_tbl <- as_matrix(count_tbl)
    if (!is.numeric(count_tbl)) stop("Count table is not numeric.
  Counts should be the first option and you used ", shQuote(cname))
  }
  # match first column in sample data by default
  if (ncol(count_tbl) != nrow(sample_tbl)) {
    stop("count_tbl should have one more column than sample_tbl rows")
  }
  n <- match(sample_tbl[[id]], colnames(count_tbl))
  if (any(is.na(n))) stop("Column names in count_tbl do not match sample names in column ", id)
  if (all(diff(n) == 1)) {
    # message("counts are already sorted by samples")
    c1 <- orig_count_tbl
  } else {
    message("Reordering columns in counts to match samples")
    if (dplyr::is.tbl(orig_count_tbl)) {
      c1 <- orig_count_tbl[, c(1, n + 1)]
    } else {
      c1 <- orig_count_tbl[, n]
    }
  }
  c1
}
