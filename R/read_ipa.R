#' Read IPA Expression Analysis
#'
#' Read an IPA Expression Analysis export all and optionally
#' write to an  Excel file with each table in a separate tab.
#'
#' @param file IPA output file (created using Export All on the main Summary tab)
#' @param excel write to Excel file
#' @param mylists include optional My Pathways and My Lists
#'
#' @return a list of tibbles or Excel file
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' x <- read_ipa("RPR2_vs_RPM.txt")
#' names(x)
#' read_ipa("RPR2_vs_RPM.txt", excel = TRUE)
#' }
#' @export

read_ipa <- function(file, excel = FALSE, mylists = FALSE) {
  x <- readr::read_lines(file)
  ## Find start of tables
  n <- grep("for My Projects", x)
  # save results
  n1 <- length(n)
  ## check...
  if (n1 == 0) stop("No My Projects headers. File should be an IPA output file from Export All")
  z <- vector("list", n1)
  ## get names for Excel tabs
  TABS <- gsub(" for My .*", "", x[n])
  names(z) <- TABS
  ## add length of x
  n <- c(n, length(x) + 1)
  # Loop through tables...
  for (i in 1:n1) {
    start <- n[i] + 1
    end <- n[i + 1] - 2
    if (!mylists) if (TABS[i] %in% c("My Lists", "My Pathways")) next
    if (end - start == 0) next
    message("Loading ", TABS[i])
    y <- x[start:end]
    ## fix 1 row with  3 columns
    if (i == 1) y <- gsub("\tand consider Both up", " and consider Both up", y)
    ## first table with Analysis Details returns Warning: 5 parsing failures.
    z[[i]] <- suppressWarnings(readr::read_tsv(paste(gsub("\t$", "", y), collapse = "\n")))
  }
  ## fix Pathways, rename columns, add total matches, set size

  z[["Canonical Pathways"]] <- dplyr::rename(z[["Canonical Pathways"]],
    Pathway = `Ingenuity Canonical Pathways`, pValue = `-log(p-value)`
  ) %>%
    dplyr::mutate(
      pValue = 10^(-pValue),
      #  zScore = as.numeric(sprintf("%.3f", zScore)),
      zScore = round(zScore, 3),
      N = nchar(gsub("[^,]", "", Molecules)) + 1,
      setSize = round(N / Ratio, 0)
    ) %>%
    dplyr::arrange(pValue)
  ## NaN read as #NUM! in Excel
  z[["Canonical Pathways"]]$zScore[is.nan(z[["Canonical Pathways"]]$zScore)] <- NA
  ## add adjusted p-value
  z[["Canonical Pathways"]]$padj <- p.adjust(z[["Canonical Pathways"]]$pValue, method = "BH")

  # drop Analysis name and empty ID
  z[["Upstream Regulators"]] <- z[["Upstream Regulators"]][, -(1:2)]
  z[["Networks"]] <- z[["Networks"]][, -2]
  # drop missing tables , usually My lists
  z <- z[!sapply(z, is.null)]
  if (excel) {
    z <- lapply(z, data.frame, check.names = FALSE)
    outfile <- gsub(".txt", ".xlsx", file)
    message("Saved to ", outfile)
    openxlsx::write.xlsx(z, file = outfile)
  } else {
    z
  }
}
