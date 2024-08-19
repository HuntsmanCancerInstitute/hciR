#' Plot counts for a single gene
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param intgroup column name in \code{attr(x, "colData")} to use for grouping
#' @param n total number of genes to plot, default 25
#' @param geom geom type, default jitter or point, boxplot and violin
#' @param ylab y-axis label
#' @param reorder level names to reorder factor levels on x-axis
#' @param ncol number of columsn for facet_wrap
#' @param scales scales for facet_wrap, default fixed
#' @param \dots additional options passed to geom_jitter or other geom
#'
#' @note See \code{link{plot_interactions}} to plot many genes.
#'
#' @return A ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#'  x <- top_counts(pasilla$results, pasilla$rlog, top=12)
#'  plot_counts(x, "condition")
#' @export

plot_counts <- function(x, intgroup, n=25, geom="jitter", ylab = "Log2 counts", reorder, ncol=NULL, scales="fixed", ...){
  if(!tibble::is_tibble(x)) stop("x should be a tibble from top_counts")
  s1  <- attr(x, "colData")
  if(is.null(s1)) stop("x should be a tibble from top_counts")
  if(!intgroup %in% names(s1))  stop('intgroup are missing from attr(x, "colData")')

  x1 <- t( as_matrix(x))
  if(ncol(x1) > n){
      # message("Displaying top ", n, " genes from top_counts")
      x1 <- x1[, 1:n]
  }
  y <- data.frame(s1[, intgroup, drop=FALSE], x1)
  z <- tidyr::gather(y, -dplyr::all_of(intgroup), key="gene2", value = "count")
   ## reorder?
   if(!missing(reorder)) z[[1]] <- factor(z[[1]], levels = reorder)

 # Note: Using an external vector in selections is ambiguous.
   #  Use `all_of(intgroup)` instead of `intgroup` to silence this message.
 if(geom=="boxplot"){
   ggplot2::ggplot(z, ggplot2::aes_string(x = names(z)[1], y = "count", color=names(z)[1])) +
    ggplot2::geom_boxplot(show.legend = FALSE, ...) +
       ggplot2::facet_wrap(~gene2, ncol=ncol, scales=scales) + ggplot2::ylab( ylab)
 }else if(geom=="violin"){
   ggplot2::ggplot(z, ggplot2::aes_string(x = names(z)[1], y = "count", color=names(z)[1])) +
    ggplot2::geom_violin(show.legend = FALSE, ...) +
     ggplot2::geom_boxplot(show.legend = FALSE, width=0.1) +
       ggplot2::facet_wrap(~gene2, ncol=ncol, scales=scales) + ggplot2::ylab( ylab)
 }else if(geom == "jitter"){
   ggplot2::ggplot(z, ggplot2::aes_string(x = names(z)[1], y = "count", group =1)) +
    ggplot2::geom_jitter( ...) +
      ggplot2::stat_summary( fun="mean", geom="line") +
       ggplot2::facet_wrap(~gene2, ncol=ncol, scales=scales) + ggplot2::ylab( ylab)
}else{
  ggplot2::ggplot(z, ggplot2::aes_string(x = names(z)[1], y = "count", group =1)) +
   ggplot2::geom_point( ...) +
     ggplot2::stat_summary( fun="mean", geom="line") +
      ggplot2::facet_wrap(~gene2, ncol=ncol, scales=scales) + ggplot2::ylab( ylab)
}



}
