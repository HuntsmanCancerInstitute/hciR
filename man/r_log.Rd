% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r_log.R
\name{r_log}
\alias{r_log}
\title{Regularized log transformation}
\usage{
r_log(dds, ...)
}
\arguments{
\item{dds}{a DESeqDataSet}

\item{\dots}{additional options passed to \code{rlog}}
}
\value{
A DESeqTransform object
}
\description{
This is a wrapper for \code{rlog} and avoids loading the DESeq2 package for
simple workflows
}
\examples{
\dontrun{
  rlog_dds(dds)
}
}
\author{
Chris Stubben
}
