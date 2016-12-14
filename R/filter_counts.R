#' Filter count matrix
#'
#' Remove features with zero counts and less than the low count cutoff.  The
#'  default method is a maximum-based filter unless sum=TRUE
#'
#' @param count_matrix a count matrix
#' @param n low count cutoff, default 1
#' @param sum Use total read counts to apply the cutoff.  The default is the maximum reads.
#'
#' @return A count matrix
#'
#' @author Chris Stubben
#'
#' @examples
#' c1 <- matrix(c(0,0,0,2,0,0,1,1,0,1,1,1), ncol=3)
#' c1
#' filter_counts(c1)
#' filter_counts(c1, sum=TRUE)
#' @export

filter_counts <- function(count_matrix,  n=1,  sum=FALSE ){
   n1 <- rowSums(count_matrix) == 0
   message( "Removed ", sum(n1), " features with 0 reads")
   c1 <- count_matrix[!n1, , drop=FALSE]
   if(sum){
     n2 <- rowSums(c1) <= n
     message( "Removed ", sum(n2), " features with <=", n, " total reads")
  }else{
     n2 <- apply(c1, 1, max) <= n
     message( "Removed ", sum(n2), " features with <=", n, " maximum reads")
  }
  c1 <- c1[!n2, , drop=FALSE]
  c1
}
