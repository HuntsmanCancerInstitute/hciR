#' Parse sample from file names
#'
#' Sample names are parsed from file names without the extension.  If the file name is not unique,
#'    the parent directory is used.
#'
#' @param files file name with path
#'
#' @return A vector
#'
#' @author Chris Stubben
#'
#' @examples
#'  # File name or parent directory should be unique
#'  sample_extract(c("align1/1355X1.counts", "align2/1355X2.counts"))
#'  sample_extract(c("align1/1355X1/Log.out", "align2/1355X2/Log.out"))
#' @export


sample_extract <- function(files){
  ## capture sample in file name  or path
  x <- lapply( strsplit(files, "/"), rev)
  samples <- sapply(x, "[", 1 )
  # remove file extension
  samples <- gsub("\\..*", "", samples)
  if( any(duplicated(samples )) ){
      # use parent directory, 13555X2/Log.final.out
      samples <- sapply(x, "[", 2 )
  }
  if( any(duplicated(samples )) ){
     stop("Sample names are not unique:  \n  ", paste(files, collapse="\n  "), call.=FALSE)
  }
  samples
}
