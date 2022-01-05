#' Read CollectRnaSeqMetrics output
#'
#' Reads CollectRnaSeqMetric output files and creates metrics and coverage tables.
#'
#' @param path the path name to CollectRnaSeqMetrics output files, the default
#'        corresponds to the working directory.
#' @param pattern a regular expression, default ".txt$". Only file names which
#'        match the regular expression will be loaded.
#'
#' @return A list with coverage and stats tables
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' # read all *.txt files in out/ directory
#' read_RnaSeqMetrics("out")
#' # read output files matching Control-1 in current directory
#' read_RnaSeqMetrics(pattern = "Control-1")
#' }
#' @export

read_RnaSeqMetrics <- function(path = ".", pattern = "\\.txt$") {
  outfiles <- list.files(path, pattern, recursive = TRUE, full.names = TRUE)
  if (length(n) == 0) stop("No ", pattern, " files found in ", path, call. = FALSE)
  samples <- sample_names(outfiles)
  samples <- gsub("_metrics$", "", samples)

  out1 <- vector("list", length(outfiles))
  out2 <- vector("list", length(outfiles))

  for (i in seq_along(outfiles)) {
    x <- readr::read_lines(outfiles[i])

    n1 <- grep("^PF_BASES", x)
    # warn and skip instead?
    if (length(n1) == 0) stop("Line starting with PF_BASES not found in ", n[i])
    x1 <- readr::read_tsv(paste(x[n1:(n1 + 1)], collapse = "\n"))
    names(x1) <- tolower(names(x1))
    ## long format, convert to numeric, remove NAs and add sample name
    out1[[i]] <- tidyr::gather(x1) %>%
      dplyr::mutate(value = as.numeric(value)) %>%
      dplyr::filter(!is.na(value)) %>%
      tibble::add_column(sample = samples[i], .before = 1)

    n2 <- grep("^normalized_position", x)
    if (length(n2) == 0) stop("Line starting with normalized_position not found in ", n[i])
    # normalized_position All_Reads.normalized_coverage
    x2 <- readr::read_tsv(paste(x[n2:length(x)], collapse = "\n"), skip = 1, col_names = c("position", "coverage"))
    out2[[i]] <- tibble::add_column(x2, sample = samples[i], .before = 1)
  }
  message("Read ", length(outfiles), " files")
  list(metrics = dplyr::bind_rows(out1), coverage = dplyr::bind_rows(out2))
}
