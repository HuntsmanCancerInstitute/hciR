#' Sort the columns of a count matrix by rows in sample data
#'
#' Creating a DESeqDataSet from a count matrix requires columns sorted by the rows in the sample dataset.
#'  This reorders the columns in the count matrix to match sample data.
#'
#' @param count_matrix a count matrix
#' @param sample_data a sample table
#' @param id position of column in sample data that matches column names in the
#'   count matrix, default is the first column
#'
#' @return A re-ordered count matrix
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' count_matrix <-  sort_counts(count_matrix, sample_data)
#' }
#' @export

sort_counts <- function(count_matrix, sample_data, id=1){
   # match first column in sample data by default
   if( ncol(count_matrix) !=  nrow(sample_data) ){
         stop("sample_data should have the same number of rows as columns in the count_matrix")
   }
   n <- match(sample_data[[id]], colnames(count_matrix ) )
   if(any(is.na(n))) stop( "Column names in count_matrix do not match names in column ", id)
   if(all(diff(n)==1) ){
       message("Count_matrix is already sorted by sample_data")
   }else{
      message("Reordering columns in count_matrix to match sample_data")
      count_matrix <- count_matrix[, n]
    }
    count_matrix
}
