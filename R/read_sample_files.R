#' Read and combine common output files
#'
#' Read output files and parse the file name to add sample IDs in the first column
#'
#' @param path the path to output files, the default corresponds to the working directory.
#' @param pattern regular expression for file name matching
#' @param delim separator in output files, default table
#' @param \dots additional options such as col_names passed to \code{read_delim}.
#'
#' @note Sample names are parsed from file names without extensions.  If the file name is not unique,
#'    the parent directory is used.
#'
#' @return A list with coverage and stats data.frames
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    #FeatureCounts summary (second column name with *.bam is always unique, so skip and assign)
#'    fc <- read_sample_files(".summary$", skip=1, col_names=c("status", "count"))
#'  filter(fc, count!=0) %>%
#'    hchart("bar", x=sample, y=count, group=status) %>%
#'     hc_plotOptions(bar = list(stacking = "normal"))
#' }
#' @export

read_sample_files <- function(path=".", pattern="\\.counts$", delim="\t",  ...){
   outfiles <- list.files(path, pattern, recursive=TRUE, full.names=TRUE)
   if(length(outfiles) == 0) stop("No ", pattern, " files found in ", path, call.=FALSE)
   samples <- sample_extract(outfiles)

   out1 <- vector("list", length(outfiles))
   for(i in seq_along(outfiles)){
       message("Reading ", outfiles[i])
       x <- suppressMessages( readr::read_delim(outfiles[i], delim=delim, ...) )
       out1[[i]] <- tibble::add_column (x, sample= samples[i], .before=1)
   }
   dplyr::bind_rows(out1)
}
