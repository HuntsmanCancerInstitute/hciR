#' Summarize tibble columns
#'
#' Summarize tibbles with many columns like TCGA metadata.  Summaries for large tibbles
#' with many rows will be slow and not recommended.
#'
#' @param x a tibble
#' @param y a vector
#'
#' @return A tibble with column names, types, number of unique values and NAs,
#' minimum, maximum and top three values.
#'
#' @author Chris Stubben
#'
#' @examples
#'  summary_tibble(mtcars)
#'  #drop_empty_columns( tibble( a=1:5, b=NA, c="x", d="") )
#' @export

summary_tibble <- function(x){
   tibble::tibble(
      column = colnames(x),
       class = vapply(x, tibble::type_sum, character(1)),   # code from glimpse
      unique = sapply(x, function(y) length( stats::na.omit( unique(y)))),
         NAs = sapply(x, function(y) sum(is.na(y) ) ),
        ### suppressWarnings to avoid :  min(c(NA, NA), na.rm=TRUE)
        min  = suppressWarnings( apply(x, 2, min, na.rm=TRUE )),
        max  = suppressWarnings( apply(x, 2, max, na.rm=TRUE )),
        ## will be slow with many rows...
       top3 = sapply(x, top3)
   )
}

#' @describeIn summary_tibble Top three values
#' @export
top3 <- function(y){
   z <- sort( table(y), decreasing = TRUE)
   if(length(z) > 3) z <- z[1:3]
   z <- paste0(names(z), " (", z, ")")
   paste(z, collapse=", ")
}

#' @describeIn summary_tibble Drop empty columns
#' @export
drop_empty_columns <- function(x){
   n1 <- apply(x, 2, function(y) all(is.na(y) | y=="") )
   if(sum(n1) > 0){
      message("Dropping ", sum(n1), " columns: ", paste( colnames(x)[n1], collapse=", ") )
      x <- x[, !n1]
   }
   x
}
