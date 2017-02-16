% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/top_genes.R
\name{top_genes}
\alias{top_genes}
\title{Get top genes}
\usage{
top_genes(res, top = 40, alpha = 0.05, basemean, log2FC, col_names, trt,
  row_names = "gene_name", difference = TRUE, sort_fc = FALSE, ...)
}
\arguments{
\item{res}{an annotated DESeq2 results file}

\item{top}{Number of top genes to display in matrix}

\item{alpha}{Adjusted p-value cutoff}

\item{basemean}{basemean cutoff}

\item{log2FC}{absolute value log2 fold change cutoff}

\item{sort_fc}{Sort top genes by fold changes and get top n/2 up and down-regulated, default is to sort by adjusted p-value}
}
\value{
A tibble of top genes
}
\description{
Return top genes in DESeq results
}
\examples{
\dontrun{
 top_genes(res)
 top_genes(res[[1]], sort_fc=TRUE)
}
}
\author{
Chris Stubben
}
