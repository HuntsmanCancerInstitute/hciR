% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_volcano.R
\name{plot_volcano}
\alias{plot_volcano}
\title{Volcano plot}
\usage{
plot_volcano(
  res,
  pvalue_cutoff,
  foldchange_cutoff,
  max_pvalue = 200,
  radius = 3,
  ggplot = TRUE,
  palette = "RdBu",
  missing = FALSE,
  labelsize = 3,
  ...
)
}
\arguments{
\item{res}{Annotated DESeq results table from results_all.}

\item{pvalue_cutoff}{either a single p-value cutoff on the y-axis or two
p-values to label down- and up-regulated genes. The default is 1.3
(corresponding to padj = 0.05) for highcharts and no labels for ggplot.}

\item{foldchange_cutoff}{either the absolute value of the log2 fold change
cutoff or negative and positive fold changes to label genes, default is 2 for
highcharts and no labels for ggplot.}

\item{max_pvalue}{y-axis limit, maximum value on a -10 log10 y-axis scale,
the default is 200 (padj < 1e-200), so genes below this cutoff are assigned
the maximum p-value.}

\item{radius}{point size, default 3}

\item{ggplot}{plot ggplot version, default TRUE}

\item{palette}{RColorBrewer palette name, vector of colors, or "RdGn" for
ggplot}

\item{missing}{Replace missing gene names with Ensembl IDs}

\item{labelsize}{label size, default 3}

\item{\dots}{other options passed to \code{hc_chart} or \code{geom_text_repel}}
}
\value{
A highchart or ggplot.
}
\description{
Plot fold changes and adjusted p-values in an interactive volcano plot or
ggplot
}
\details{
For ggplot, the results should not be sorted by p-value or fold change to
avoid stacking close overlapping points.  Labels are added using
\code{ggrepel}, so avoid labeling too many points (200 is the limit).
}
\examples{
plot_volcano(pasilla$results, pvalue=c(35,25), foldchange=2.5)
}
\author{
Chris Stubben
}
