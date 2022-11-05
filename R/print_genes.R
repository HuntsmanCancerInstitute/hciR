#' Print the scaled values from plot_genes and optionally add branches from cutree
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param cut number of branches to cut
#'
#' @return tibble with id, branch and rlog values used by \code{\link{plot_genes}}
#'
#' @author Chris Stubben
#'
#' @examples
#' x <- top_counts(pasilla$results, pasilla$rlog)
#' plot_genes(x, c("condition", "type"), scale="row", annotation_names_col=FALSE)
#' print_genes(x, cut=5)
#' @export

print_genes <- function(x, cut){
  # scale by rows
  y <- scale(t(as_matrix(x)))
  x1 <- as_tbl(t(y))
  attr(x1, "colData") <- attr(x, "colData")
  p <- plot_genes(x1, scale="none", silent=TRUE)
  ## reorder rows and columns in matrix to match heatmap
  y2 <- y[p$tree_col$order, p$tree_row$order]
  y2 <- round(y2,3)
  ## convert to tibble (easier to print and save)
  y2 <- as_tbl(t(y2))
  if(!missing(cut)){
     ## cut column dendrogram
     n <- dendextend::cutree(p$tree_row, cut, order=FALSE)
     y2 <- tibble::add_column(y2, branch=n, .after="id")
  }
  y2
}
