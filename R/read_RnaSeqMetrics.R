#' Read CollectRnaSeqMetrics output
#'
#' Reads CollectRnaSeqMetric output files and creates metrics and coverage tables.
#'
#' The path and pattern options are used by \code{\link{list.files}}.
#'
#' @param path the path name to CollectRnaSeqMetrics output files; the default
#'        corresponds to the working directory.
#' @param pattern a regular expression, default ".txt$". Only file names which
#'        match the regular expression will be loaded.
#' @param split a regular expression to parse sample name from file name, default is name before _ or . "[_.]".
#' @param \dots additional options passed to \code{\link{list.files}}
#'
#' @return A list with coverage and stats data.frames
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    # read all *.txt files in out/ directory
#'    read_RnaSeqMetrics( "out")
#'    # read output files matching Control-1 in current directory
#'    read_RnaSeqMetrics( pattern="Control-1")
#'    # read all files in out/ directory and all subdirectories
#'    read_RnaSeqMetrics(out", "*", recursive=TRUE)
#' }

read_RnaSeqMetrics<- function( path=".", pattern="\\.txt$",  split= "[_.]", ...){
  n <- list.files(path, pattern, ...)
  if(length(n)==0) stop("No files found")

  out1 <- vector("list", length(n))
  out2 <- vector("list", length(n))

  ## if recursive=TRUE in list.files, then drop subdirectory
  n1 <- gsub(".+/", "", n)
  # use  first part of file name as sample
  samples <- gsub(paste0( "(.+?)", split, ".*") , "\\1", n1)
  ## check if not unique
  if( length(unique(samples)) < length(samples) ){
    message("Warning: first part of file name before split is not unique")
     samples <-  n1
  }

  for(i in seq_along(n)){
    x <- readLines( paste(path, n[i], sep="/"))

    n1 <- grep("^PF_BASES", x)
  # warn and skip instead?
    if( length(n1) == 0)  stop("Line starting with PF_BASES not found in ", n[i] )

    zz <- textConnection(x[n1:(n1+1)])
    x1 <- read.delim(zz)
    close(zz)
    ## long format
    y <- data.frame(Sample= samples[i],
             Key= paste0(substr(names(x1), 1,1), tolower(substring(names(x1), 2) )),
             Value= as.vector(unlist(x1)), stringsAsFactors=FALSE)
    y <- subset(y, !is.na(Value))
    out1[[i]] <- y

    n2 <- grep("^normalized_position", x)
    if( length(n2) == 0)  stop("Line starting with normalized_position not found in ", n[i] )

    zz <- textConnection(x[n2:length(x)])
    x2 <- read.delim(zz)
    close(zz)
    names(x2) <- c("Position", "Coverage")
    out2[[i]] <- cbind(Sample=samples[i], x2)
  }
  message("Read ", length(n), " files")
  stats <- do.call("rbind", out1)
  coverage <- do.call("rbind", out2)
  list(metrics=stats, coverage=coverage)
}
