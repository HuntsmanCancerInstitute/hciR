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
     ## number [Character] number
     n1 <- gsub("[A-Za-z].*", "",  samples )
     n2 <-  gsub(".*[A-Za-z]", "", samples )
    #  sprintf("%06s", n2) will pad zeros on Mac and blank on Linux
    order( as.numeric(paste0(n1,sprintf("%06d", as.numeric(n2)))) )
}
