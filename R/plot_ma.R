#' Highchart version of MA-plot
#'
#' Plot mean normalized counts and fold changes in an interactive MA-plot
#'
#' @param res Annotated DESeq results table from results_all
#' @param signif significance level, default 0.05
#' @param  baseMean normalized count cutoff for labeling points, default 10000
#' @param foldchange absolute value of log2 fold change cutoff for labeling
#' points, default 2
#' @param radius highchart point size, default 3
#' @param size ggplot point size, default 1
#' @param alpha ggplot alpha transparency, default 0.3
#' @param sig_alpha use transparency for significant genes, default TRUE
#' @param ylab y-axis label
#' @param ggplot plot ggplot version
#' @param \dots other options like width passed to \code{hc_chart}
#'
#' @return A highchart or ggplot. Only points above the baseMean or log2 fold
#' change cutoffs have mouseover labels.
#'
#' @author Chris Stubben
#'
#' @examples
#' plot_ma(pasilla$results, ggplot = TRUE)
#' @export

plot_ma <- function(res, signif = 0.05, baseMean = 10000, foldchange = 2, radius = 3, size = 1,
                    alpha = 0.3, sig_alpha = TRUE, ylab = "Log2 Fold Change", ggplot = TRUE, ...) {
  if (!tibble::is_tibble(res)) {
    if (is.list(res)) {
      message("Plotting the first table in the list")
      res <- res[[1]]
    } else {
      stop("Results should be a tibble")
    }
  }
  ## limma top_tibble? TO DO: don't use log10 for AveExpr!
  n <- match(c("Gene.Symbol", "logFC"), names(res))
  if (!any(is.na(n))) {
    names(res)[n] <- c("gene_name", "log2FoldChange")
    res$gene_name <- gsub(", .*", "", res$gene_name)
  }
  x <- dplyr::filter(res, !is.na(log2FoldChange))
  ## plot(log10(x$baseMean), x$log2FoldChange, col=rgb(0,0,1,.3), pch=19)
  if (ggplot) {
    # limma top tibble
    if ("AveExpr" %in% names(x)) {
      x$x1 <- x$AveExpr
      xlab1 <- "Average Expression"
    } else {
      x$x1 <- log10(x$baseMean)
      xlab1 <- "Log10 Mean Normalized Counts"
    }
    # TO DO - label genes
    x <- dplyr::mutate(x, de = ifelse(padj < signif & !is.na(padj), "sig", "ns"))
    # option for transparency on significant genes
    if (!sig_alpha) alpha <- ifelse(x$de == "sig", 1, alpha)
    ggplot2::ggplot(data = x, ggplot2::aes(x = x1, y = log2FoldChange, color = de, shape = de)) +
      ggplot2::geom_point(alpha = alpha, size = size, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = c("gray20", "red")) +
      ggplot2::xlab(xlab1) +
      ggplot2::ylab(ylab) +
      ggplot2::theme_light()
  } else {
    ### Grouping column for enableMouseTracking
    if ("AveExpr" %in% names(res)) stop("Only ggplot=TRUE for limma top table")
    x$sig <- ifelse(x$baseMean > baseMean | abs(x$log2FoldChange) > foldchange,
      "Y", "N"
    )
    n <- sum(x$sig == "Y")
    if (n == 0) stop("No points above cutoffs, need to fix hchart")
    message(
      "Adding mouseover labels to ", n, " genes (",
      round(n / nrow(x) * 100, 1), "%)"
    )
    ## use Ensembl ID if gene_name is missing ?
    n <- is.na(x[["gene_name"]]) | x[["gene_name"]] == ""
    if (sum(n) > 0) {
      message("Missing ", sum(n), " gene names, using Ensembl IDs instead")
      x[["gene_name"]][n] <- x[["id"]][n]
    }
    highcharter::hchart(x, "scatter", highcharter::hcaes(log10(baseMean),
      log2FoldChange,
      group = sig, value = gene_name
    ),
    color = "rgba(0,0,255, 0.3)", enableMouseTracking = c(FALSE, TRUE),
    showInLegend = FALSE, marker = list(radius = radius)
    ) %>%
      highcharter::hc_tooltip(
        pointFormat = "{point.value}",
        headerFormat = ""
      ) %>%
      highcharter::hc_xAxis(
        title = list(text = "Log10 Mean Normalized Counts"),
        gridLineWidth = 1, tickLength = 0, startOnTick = "true",
        endOnTick = "true"
      ) %>%
      highcharter::hc_yAxis(title = list(text = ylab)) %>%
      highcharter::hc_chart(zoomType = "xy", ...) %>%
      highcharter::hc_exporting(enabled = TRUE, filename = "MA-plot")
  }
}
