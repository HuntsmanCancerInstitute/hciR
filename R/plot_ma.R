#' Highchart version of MA-plot
#'
#' Plot mean normalized counts and fold changes in an interactive MA-plot
#'
#' @param res Annotated DESeq results table from results_all
#' @param  baseMean normalized count cutoff for labeling points, default 10000
#' @param log2FoldChange absolute value of log2 fold change cutoff for labeling points, default 2
#' @param radius point size, default 2
#' @param ggplot plot ggplot version
#' @param \dots other options like width passed to \code{hc_chart}
#'
#' @return A highchart. Only points above the baseMean and log2 fold change cutoffs have mouseover labels.
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' plot_ma(res)
#' }
#' @export

plot_ma <- function(res, baseMean = 10000 , log2FoldChange = 2, radius=2, ggplot=FALSE, ...){
   if(!is_tibble(res)){
      if(is.list(res)){
        message("Plotting the first table in the list")
        res <- res[[1]]
      }else{
         stop("Results should be a tibble")
      }
   }
   x <- dplyr::filter(res, !is.na(log2FoldChange))
   ## plot(log10(x$baseMean), x$log2FoldChange, col=rgb(0,0,1,.3), pch=19)

   if(ggplot){
   ggplot2::ggplot(data=x, aes(x=log10(baseMean), y= log2FoldChange )) +
        geom_point(color="blue", alpha=0.3, size=1) +
        xlab("Log10 Mean Normalized Counts") + ylab("Log2 Fold Change")
  }else{
   ### Grouping column for enableMouseTracking
   x$sig = ifelse( x$baseMean > baseMean | abs(x$log2FoldChange)> log2FoldChange, "Y", "N")
   n <- sum(x$sig == "Y")
   if(n ==0) stop("No points above cutoffs, need to fix hchart to plot this case")
    message("Adding mouseover labels to ", n, " genes (",  round( n/nrow(x)*100, 1), "%)")
    ## use Ensembl ID if gene_name is missing ?
     n <- is.na(x[["gene_name"]]) | x[["gene_name"]] == ""
     if(sum(n)>0 ){
         message("Missing ", sum(n), " gene names, using Ensembl IDs instead")
         x[["gene_name"]][n] <- x[["id"]][n]
      }
    highcharter::hchart(x, "scatter", highcharter::hcaes( log10(baseMean), log2FoldChange,
            group = sig, value = gene_name), color = 'rgba(0,0,255, 0.3)',
             enableMouseTracking = c(FALSE, TRUE), showInLegend=FALSE, marker = list(radius = radius)) %>%
        highcharter::hc_tooltip( pointFormat = "{point.value}", headerFormat = "") %>%
         highcharter::hc_xAxis(title = list(text = "Log10 Mean Normalized Counts"), gridLineWidth = 1, tickLength = 0, startOnTick = "true", endOnTick = "true") %>%
         highcharter::hc_yAxis(title = list(text = "Log2 Fold Change")) %>%
         highcharter::hc_chart(zoomType = "xy", ...) %>%
         highcharter::hc_exporting(enabled=TRUE, filename = "MA-plot")
   }
}
