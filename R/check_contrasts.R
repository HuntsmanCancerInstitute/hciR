#' Check contrasts
#'
#' Check all possible contrasts used by \code{\link{results_all}}.  Following
#' DESeq2 conventions, control group should be listed first and then pairswise
#' comparisons are selected in reverse order, so A, B, C levels will return
#' C vs. B, C vs. A, and B vs. A
#'
#' @param trt a vector
#' @param vs either compare all vs. all (default) or a all vs. a specific treatment, or see note.
#'
#' @note If you combine factors of interest into a single group following section 3.3 in the DESeq2 vignette,
#'  you can set vs = "combined" to limit the comparisons to those within the first group.
#'
#' @return A vector of contrasts
#'
#' @author Chris Stubben
#'
#' @examples
#' # check_contrasts(samples$trt)
#' trt <- c("B", "A", "C")
#' # sorted alphabetically and then reverse order combn
#' check_contrasts(trt)
#' # or specify factor levels with control first
#' check_contrasts(factor(trt, levels = c("C", "B", "A")))
#' # combine 2 trt groups  (66 possible contrasts)
#' trt <- paste(1:4, rep(c("A", "B", "C"), each = 4))
#' check_contrasts(trt, vs = "combined")
#' @export

check_contrasts <- function(trt, vs = "all") {
  if (is.factor(trt)) {
    n <- levels(trt)
  } else {
    # DESeq2 will create factor for character strings in alphabetical order
    n <- sort(unique(trt))
  }
  n <- rev(n)
  contrast <- utils::combn(n, 2)
  if (vs == "combined") {
    ## if two columns are combined into a single trt group, compare within first group
    n1 <- apply(contrast, 2, function(x) length(unique(gsub("[ _-].+", "", x))) == 1)
    contrast <- contrast[, n1]
  } else if (vs %in% n) {
    contrast <- rbind(vs, n[n != vs])
  }
  vs <- apply(contrast, 2, paste, collapse = " vs. ")
  message(length(vs), " contrasts: ")
  vs
}
