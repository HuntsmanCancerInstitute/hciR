#' Highchart version of volcano plot
#'
#' Plot fold changes and adjusted p-values in an interactive volcano plot
#'
#' @param res Annotated DESeq results table from results_all
#' @param  padj p-value cutoff for labeling points, default 0.05
#' @param log2FoldChange absolute value of log2 fold change cutoff for labeling points, default 2
#' @param radius point size, default 4
#' @param ggplot plot ggplot version
#' @param \dots other options like width passed to \code{hc_chart}
#'
#' @return A highchart. Only points above the adjusted p-value and log2 fold change cutoffs have mouseover labels.
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' plot_volcano(res)
#' }
#' @export

plot_volcano <- function(res, padj = 0.05 , log2FoldChange = 2, radius=4, ggplot=FALSE, ...){
   if(!tibble::is_tibble(res)){
      if(is.list(res)){
        message("Plotting the first table in the list")
        res <- res[[1]]
      }else{
         stop("Results should be a tibble")
      }
   }
   x <- dplyr::filter(res, !is.na(padj))
   ## center at zero...
   fc <- max(abs(x$log2FoldChange), na.rm=TRUE)

   if(ggplot){
   ggplot2::ggplot(data=x, ggplot2::aes(x=log2FoldChange, y= -log10(padj) )) +
        ggplot2::geom_point(alpha=0.4, size=1.5, color= "blue") +
        ggplot2::xlab("Log2 Fold Change") + ggplot2::ylab("-Log10 Adjusted P-value")  +
        ggplot2::xlim( -fc, fc)
  }else{
   ### Grouping column for enableMouseTracking
   x$sig = ifelse( x$padj  > padj & abs(x$log2FoldChange) < log2FoldChange, "N", "Y")
   n <- sum(x$sig == "Y")
   if(n ==0) stop("No points above cutoffs, need to fix hchart to plot this case")
    message("Adding mouseover labels to ", n, " genes (",  round( n/nrow(x)*100, 1), "%)")
    ## use Ensembl ID if gene_name is missing ?
     n <- x$sig=="Y" & (is.na(x[["gene_name"]]) | x[["gene_name"]] == "")
     if(sum(n)>0 ){
         message("Missing ", sum(n), " gene names, using Ensembl IDs instead")
         x[["gene_name"]][n] <- x[["id"]][n]
      }

    highcharter::hchart(x, "scatter", highcharter::hcaes(log2FoldChange,  -log10(padj),
                 group = sig, value = gene_name), color = 'rgba(0,0,255, 0.3)',
             enableMouseTracking = c(FALSE, TRUE), showInLegend=FALSE, marker = list(radius = radius)) %>%
        highcharter::hc_tooltip( pointFormat = "{point.value}", headerFormat = "") %>%
         highcharter::hc_xAxis(title = list(text = "Log2 Fold Change"), gridLineWidth = 1, tickLength = 0, startOnTick = "true", endOnTick = "true" , min= -fc, max=fc) %>%
         highcharter::hc_yAxis(title = list(text = "-Log10 Adjusted P-value")) %>%
         highcharter::hc_chart(zoomType = "xy", ...) %>%
         highcharter::hc_exporting(enabled=TRUE, filename = "volcano")
      }
}
