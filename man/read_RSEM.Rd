% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_RSEM.R
\name{read_RSEM}
\alias{read_RSEM}
\title{Read RSEM output files}
\usage{
read_RSEM(
  path = ".",
  gene = TRUE,
  value = "expected_count",
  reshape = TRUE,
  stats = FALSE
)
}
\arguments{
\item{path}{the path to RSEM output files, the default corresponds to the working directory.}

\item{gene}{read gene or transcript (isoform) result files, default gene}

\item{value}{a string with the RSEM column name to populate cell values, default expected_count}

\item{reshape}{reshape into wide format with samples in rows (a count matrix).}

\item{stats}{read stat files, default counts}
}
\value{
A tibble in long or wide format if reshape=TRUE
}
\description{
Reads RSEM counts or stats files and optionally reshape into wide format.
}
\note{
The cnt and model files in the stats directory vary depending on RSEM options
and the parser may fail.
}
\examples{
\dontrun{
 # count matrix
rsem_counts <- read_RSEM( "Alignments")
 # read TPM
tpm <- read_RSEM("Alignments", value="TPM")
 # reshape uses alignments stats only (rows 1-3 in *.cnt files)
read_RSEM( ".", stats=TRUE)
}
}
\references{
\url{https://github.com/deweylab/RSEM/blob/master/cnt_file_description.txt} and
\url{https://github.com/deweylab/RSEM/blob/master/model_file_description.txt} for output format
}
\author{
Chris Stubben
}
