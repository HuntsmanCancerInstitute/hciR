#' Filter top count matrix
#'
#' Remove columns from count matrix returned by \code{link{top_counts}} by subsetting the colData attribute
#'
#' @param x a tibble from \code{link{top_counts}}
#' @param \dots additional options passed to \code{subset}
#'
#' @return A count matrix
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  filter_top_counts( top_counts(res, rld),  type == "paired-end")
#' }
#' @export

filter_top_counts <- function(x , ...){
      s1 <-  subset( attr(x, "colData"), ... )
      x1 <- x[, c(1, which( colnames(x) %in% rownames(s1) )) ]
     attr(x1, "colData") <- s1
     x1
}
