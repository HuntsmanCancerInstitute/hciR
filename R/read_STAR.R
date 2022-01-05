#' Read STAR log files
#'
#' Read STAR Log.final.out files and optionally reshape into wide format.
#'
#' @param path the path to STAR log files, the default corresponds to the working directory.
#' @param pattern regular expression for file name matching, default .final.out
#' @param reshape reshape percent mapping into wide format with samples in rows
#'
#' @note Reading output files requires a unique sample identifier in either the file name or parent directory
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' x <- read_STAR(pattern = ".star.out$")
#' reads <- c(
#'   "Uniquely mapped", "Mapped to multiple loci",
#'   "Mapped to too many loci", "Unmapped reads"
#' )
#' y <- filter(x, stat %in% reads) %>% mutate(stat = factor(stat, levels = reads))
#' hchart(y, "bar", x = sample, y = value, group = stat) %>%
#'   hc_plotOptions(bar = list(stacking = "normal")) %>%
#'   hc_colors(c("#437bb1", "#7cb5ec", "#f7a35c", "#b1084c")) %>%
#'   hc_yAxis(reversedStacks = FALSE)
#' read_STAR(pattern = ".star.out$", reshape = TRUE)
#' }
#' @export

read_STAR <- function(path = ".", pattern, reshape = FALSE) {
  if (missing(pattern)) pattern <- "\\.final.out$"
  # suppress warnings about subheaders....  Warning: 4 parsing failures.
  x <- suppressWarnings(read_sample_files(path, pattern, col_names = c("stat", "value"), skip = 5, trim_ws = TRUE))

  x <- dplyr::filter(x, !is.na(value)) %>% # drop subheader with NA value ... UNIQUE READS:
    dplyr::mutate(
      stat = gsub(" \\|$", "", stat), # remove pipe from end of statistic
      stat = gsub("Number of reads m", "M", stat), # shorten  Number of reads mapped to...
      stat = gsub(" reads number", "", stat), #  shorten Uniquely mapped reads number
      value = gsub("%", "", value), # drop % from value
      value = as.numeric(value)
    ) # change value to numeric

  # ADD unmapped reads  -
  mapped <- c("Uniquely mapped", "Mapped to multiple loci", "Mapped to too many loci")
  # see http://stackoverflow.com/questions/40749742/add-missing-subtotals-to-each-group-using-dplyr
  x <- dplyr::group_by(x, sample) %>%
    dplyr::summarize(
      value = value[stat == "Number of input reads"] - sum(value[stat %in% mapped]),
      stat = "Unmapped reads"
    ) %>%
    dplyr::bind_rows(x) %>%
    dplyr::select(sample, stat, value) %>% # back to original column order
    dplyr::arrange(sample)

  ##  split unmapped reads into too many mismatches, too short and  other?
  #  divide %unmapped too short  by total %unmapped  and muliply by Unmapped reads to get estimate

  if (reshape) {
    n <- c(mapped, "Unmapped reads", "Number of input reads")
    x <- dplyr::filter(x, stat %in% n) %>%
      dplyr::mutate(stat = factor(stat, levels = n)) %>%
      tidyr::spread(stat, value) %>%
      dplyr::mutate_each(dplyr::funs(as.integer), -1)
  }
  x
}
