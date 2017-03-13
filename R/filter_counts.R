#' Filter count matrix
#'
#' Remove features with zero counts and less than the low count cutoff.  The
#'  default method is a maximum-based filter unless sum=TRUE
#'
#' @param count_tbl a count matrix, data.frame or tibble
#' @param n low count cutoff, default 1
#' @param sum Use total read counts to apply the cutoff.  The default is the maximum reads.
#'
#' @return A count matrix
#'
#' @author Chris Stubben
#'
#' @examples
#' c1 <- matrix(c(0,0,0,2,12,0,0,1,1,0,0,1,1,1,0), ncol=3)
#' c1
#' filter_counts(c1)
#' filter_counts(c1, sum=TRUE)
#' @export

filter_counts <- function(count_tbl,  n=1,  sum=FALSE ){
   orig_count_tbl <- count_tbl
   if( class(count_tbl)[1] != "matrix") count_tbl <- as_matrix(count_tbl)

   n1 <- rowSums(count_tbl) == 0
   message( "Removed ", sum(n1), " features with 0 reads")
   if(sum){
     n2 <- rowSums(count_tbl) <= n  & !n1
     message( "Removed ", sum(n2), " features with <=", n, " total reads")
  }else{
     n2 <- apply(count_tbl, 1, max) <= n  & !n1
     message( "Removed ", sum(n2), " features with <=", n, " maximum reads")
  }
  c1 <- orig_count_tbl[!(n1|n2), , drop=FALSE]
  c1
}
