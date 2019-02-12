#' Sort sample names
#'
#' Sort samples names using the first and last run of numbers in a string
#'
#' @param samples a vector of sample names
#'
#' @return A vector
#'
#' @author Chris Stubben
#'
#' @examples
#'  ids <- c("135X1", "135X21", "135X2", "135X10", "81X5", "81X15")
#'  sort(ids)
#'  sort_samples(ids)
#' @export

sort_samples <- function(samples){
   samples[order_samples(samples)]
}

#' @describeIn sort_samples Order sample names
#' @export
order_samples <- function( samples ){
     ## check if number [Character] number
     need_to_sort <- all(grepl("^[0-9]+[A-Za-z]+[0-9]+$", samples))
     if(need_to_sort){
       n1 <-  gsub("[A-Za-z]+.*", "",  samples )
       n2 <-  gsub(".*[A-Za-z]+", "", samples )
       ## check for _ or other characters
       n1 <- as.numeric(gsub("[^0-9]", "", n1))
       n2 <- as.numeric( gsub("[^0-9]", "", n2))
       n <-  order( n1,n2 )
   }else{
      n <- order(samples)
   }
   n
}
