#' Plot counts for many genes
#'
#' @param x a tibble from \code{\link{top_counts}}
#' @param intgroups two column names in \code{attr(x, "colData")} to use for grouping
#' @param ylab y-axis label
#' @param scaled scale counts
#' @param n total number of genes to plot
#' @param reorder level names to reorder factor levels on x-axis
#'
#' @return A ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- top_counts( res, rld, top=25)
#'  plot_interactions(x,  c("time", "trt"), ylab = "scaled rlog")
#' }
#' @export

plot_interactions <- function(x, intgroups, ylab = "count", scaled = TRUE, n=40, reorder){
   if(length(intgroups) != 2) stop( "intgroups should be a vector of two column names")
   if(!tibble::is_tibble(x)) stop("x should be a tibble from top_counts")
   s1  <- attr(x, "colData")
   if(is.null(s1)) stop("x should be a tibble from top_counts")
   if(!all( intgroups %in% names(s1))) stop('intgroups are missing from attr(x, "colData")')
   x1 <- t( as_matrix(x))
   if(scaled) x1 <- scale(x1)
   if(ncol(x1) > n){
      message("Displaying top ", n, " genes from top_counts")
       x1 <- x1[, 1:n]
   }
   x1 <- data.frame( x1)
   y <- split(x1,  s1[, intgroups] )
   y <- lapply(y, colMeans, na.rm=TRUE)
   z <- dplyr::bind_rows( lapply(y, tibble::as_tibble), .id="stage") %>%
     tidyr::separate(stage, into= intgroups, sep="\\.")
    ## reorder?
    if(!missing(reorder)) z[[1]] <- factor(z[[1]], levels = reorder)
   z <-    dplyr::mutate(z, gene= rep( names(y[[1]]),  length(y) )) %>%
      filter(!is.na(value))
   ggplot2::ggplot(z, ggplot2::aes_string(x = intgroups[1], y = "value",
      group =intgroups[2], shape =intgroups[2], color=intgroups[2])) +
       ggplot2::geom_line() +  ggplot2::geom_point() +
        ggplot2::facet_wrap(~gene) + ggplot2::ylab( ylab)

}
