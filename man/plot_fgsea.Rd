% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_fgsea.R
\name{plot_fgsea}
\alias{plot_fgsea}
\title{Plot fGSEA output}
\usage{
plot_fgsea(
  x,
  trim = 70,
  sets,
  nes = TRUE,
  cluster_row = FALSE,
  cluster_col = FALSE,
  ...
)
}
\arguments{
\item{x}{list from \code{link{fgsea_all}}}

\item{trim}{trim long names, default more than 70 characters}

\item{sets}{display contrasts sharing n or more sets for n > 1.  If n = 1,
then only plot unique sets.  If missing, then plots all sets, default.}

\item{nes}{plot NES (or ES if FALSE)}

\item{cluster_row}{Cluster dendrogram rows, default is an alphabetical list}

\item{cluster_col}{Cluster dendrogram columns, default FALSE}

\item{\dots}{other options passed to \code{pheatmap}}
}
\description{
Plot heatmap with enrichment scores by contrast
}
\examples{
\dontrun{
  library(hciRdata)
  fc <- write_gsea_rnk(res, write=FALSE)
  x <- fgsea_all(fc, msig_pathways$KEGG, FDR= 0.25)
  plot_fgsea(x)
 }
}
\author{
Chris Stubben
}
