#' Plot counts for many genes
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param intgroups two column names in \code{attr(x, "colData")} to use for
#' grouping
#' @param ylab y-axis label
#' @param scaled scale counts
#' @param n total number of genes to plot
#' @param max_scale maximum values for scaled values in case of outliers
#' @param reorder level names to reorder factor levels on x-axis
#' @param \dots additional options like ncol or nrow passed to facet_wrap
#'
#' @return A ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#'  x <- top_counts( pasilla$results, pasilla$rlog, top=25)
#'  # no interaction from single or paired end as expected
#'  plot_interactions(x,  c("condition", "type"))
#' @export

plot_interactions <- function(x, intgroups, ylab = "scaled rlog", scaled = TRUE,
 n=40, max_scale, reorder, ...){
   if(length(intgroups) != 2) stop( "intgroups should be a vector of two names")
   if(!tibble::is_tibble(x)) stop("x should be a tibble from top_counts")
   s1  <- attr(x, "colData")
   if(is.null(s1)) stop("x should be a tibble from top_counts")
   if(!all( intgroups %in% names(s1))){
       stop('intgroups are missing from attr(x, "colData")')
   }
   x1 <- t( as_matrix(x))
   if(scaled) x1 <- scale(x1)
   if(!missing(max_scale)){
       x1[x1 > max_scale] <-  max_scale
       x1[x1 < -max_scale] <-  -max_scale
   }
   if(ncol(x1) > n){
       # message("Displaying top ", n, " genes from top_counts")
       x1 <- x1[, 1:n]
   }

   y <- data.frame(s1[, intgroups], x1)
   z <- tidyr::gather(y, -intgroups, key="gene2", value = "rlog")
    ## reorder?
    if(!missing(reorder)) z[[1]] <- factor(z[[1]], levels = reorder)

	# Note: Using an external vector in selections is ambiguous.
    #  Use `all_of(intgroups)` instead of `intgroups` to silence this message.

  ## add option to drop geom_poins (only fitted line)
   ggplot2::ggplot(z, ggplot2::aes_string(x = names(z)[1], y = "rlog",
      group =names(z)[2], shape =names(z)[2], color=names(z)[2])) +
     ggplot2::geom_point() +
       ggplot2::stat_summary( fun="mean", geom="line") +
        ggplot2::facet_wrap(~gene2, ...) + ggplot2::ylab( ylab)

}
