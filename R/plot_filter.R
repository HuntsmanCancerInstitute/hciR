#' Plot count matrix to check filter cutoff
#'
#' Plot the total number of features removed from count matrix using different
#' low count cutoffs for both maximum and total reads
#'
#' @param x count matrix, data.frame or tibble
#' @param n maximum cutoff, default is to display 0 to 20
#'
#' @return A plot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' plot_filter(counts)
#' }
#' @export

plot_filter <- function(x, n = 20) {
  if (dplyr::is.tbl(x)) x <- as_matrix(x)
  n1 <- rowSums(x)
  n2 <- apply(x, 1, max, na.rm = TRUE)
  x1 <- table(factor(n1[n1 <= n], levels = 0:n))
  x2 <- table(factor(n2[n2 <= n], levels = 0:n))
  y <- cbind(x1, x2)
  z <- apply(y, 2, cumsum)
  graphics::matplot(rownames(z), z,
    pch = c(17, 19), col = c("red", "blue"),
    xlab = "Count cutoff", ylab = "Total features removed"
  )
  graphics::legend("bottomright", c("max", "total"),
    pch = c(19, 17),
    col = c("blue", "red"), bty = "n", inset = 0.1, title = "Filter"
  )
}
