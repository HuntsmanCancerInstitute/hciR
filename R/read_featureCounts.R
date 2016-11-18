#' Read featureCounts output files
#'
#' Reads featureCount count or alignment summary files and optionally reshape into
#'  wide format.
#'
#' @param path the path to featureCounts output files, the default
#'        corresponds to the working directory.
#' @param pattern regular expression for file name matching, deafult .counts and .summary
#' @param reshape reshape into wide format with samples in rows (count matrix)
#' @param stats read stats tables, default counts
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  ft_counts <- read_featureCounts(".")
#'  fc <- read_featureCounts(".", stats=TRUE)
#'  fc
#'  fc <- read_featureCounts(".", stats=TRUE, reshape=FALSE)
#'  fc
#'  filter(fc, count!=0) %>%
#'    hchart("bar", x=sample ,  y=count, group=status) %>%
#'      hc_plotOptions(bar = list(stacking = "normal"))
#' }
#' @export

read_featureCounts <- function( path=".", pattern, reshape=TRUE, stats=FALSE){
   if(stats){
      if(missing(pattern))  pattern <- "\\.summary$"
     # second column name is always unique, so skip header and assign column names
      fc <- read_sample_files(path, pattern, col_names=c("status", "count"), skip=1)
      if(reshape)  fc <- dplyr::filter(fc, count!=0) %>% tidyr::spread(status, count)
   }
   else{
     if(missing(pattern))  pattern <- "\\.counts$"
      fc <- read_sample_files(path, pattern, col_names=c("geneid",	"chr", "start",	"end", "strand", "length", "count"), skip=2)
      if(reshape)  fc  <- dplyr::select(fc, sample, geneid, count) %>% tidyr::spread(sample, count)
   }
   fc
}
