#' Read samtools idxstats output
#'
#' Read mapped reads by sequence in samtools idxstats output
#'
#' @param path the path to idxstats output files, the default
#'        corresponds to the working directory.
#' @param pattern regular expression for file name matching, default .idxstats
#' @param reshape reshape mapped reads into wide format with samples in rows
#' @param fragment keep sequence names like GL000008.2 with 8 or more characters
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' read_idxstats(pattern = "idxstats.txt$")
#' }
#' @export

read_idxstats <- function(path = ".", pattern, reshape = FALSE, fragment = FALSE) {
  if (missing(pattern)) pattern <- "\\.idxstats$"
  # second column name is always unique, so skip header and assign column names
  idx <- read_sample_files(path, pattern, col_names = c("sequence", "length", "mapped", "unmapped"))
  # drop fragments
  if (!fragment) idx <- subset(idx, nchar(sequence) < 8)
  # add reads per kb and percent mapped
  idx <- dplyr::filter(idx, sequence != "*") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      reads_kb = round(mapped / length * 1000, 1),
      percent_map = round(mapped / sum(mapped) * 100, 3)
    )
  if (reshape) idx <- dplyr::select(idx, sample, sequence, mapped) %>% tidyr::spread(sequence, mapped)
  idx
}
