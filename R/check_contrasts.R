#' Check contrasts
#'
#' Check all possible contrasts used by \code{\link{results_all}}
#'
#' @param trt a vector
#' @param vs either compare all vs. all (default) or a all vs. a specific treatment, or see note.
#' @param vs2 position of specific treatment in contrast vector, set FALSE for specific treatment vs all
#'
#' @note If you combine factors of interest into a single group following section 3.3 in the DESeq2 vignette,
#'  you can set vs = "combined" to limit the comparisons to those within the first group.
#'
#' @return A vector of contrasts
#'
#' @author Chris Stubben
#'
#' @examples
#'  # check_contrasts( samples$trt)
#'  trt <- factor(c( "control", "heart", "lung"))
#'  check_contrasts(trt)
#'  trt <- factor(trt, levels = c(  "heart", "lung", "control"))
#'  check_contrasts(trt, control_first =FALSE)
#'  data.frame( vs=check_contrasts(trt), control_first =FALSE)
#'  check_contrasts(trt, vs = "control")
#'  check_contrasts(trt, vs = "control", vs2=FALSE)
#'  check_contrasts(factor(trt, levels = trt)  )
#'  # combine 2 trt groups  (66 possible contrasts, 12 within 1st group).
#'  trt <- paste( 1:4, rep(c("A", "B", "C"), each=4))
#'  check_contrasts( trt, vs = "combined")
#' @export

check_contrasts <- function( trt, vs="all", vs2 = TRUE, control_first = FALSE){
   if(is.factor(trt)){
      n <- levels( trt )
   }else{
	   # DESeq2 will create factor for characters
      n <- sort(unique(trt))
   }
   if(control_first) n <- rev(n)
   contrast <- utils::combn(n, 2)
   if( vs == "combined"){
      ## if two columns are combined into a single trt group, compare within first group
      n1 <- apply(contrast, 2, function(x) length(unique( gsub("[ _-].+", "", x)))==1)
      contrast <- contrast[, n1]
   }else if( vs %in% n){
      if(vs2){
         contrast <- rbind( n[n!=vs], vs)
      }else{
         contrast <- rbind( vs, n[n!=vs])
      }
   }
   vs <- apply(contrast, 2, paste, collapse = " vs. ")
   message(length(vs), " contrasts: ")
   vs
}
