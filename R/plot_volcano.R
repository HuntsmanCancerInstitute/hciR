#' Highchart version of volcano plot
#'
#' Plot fold changes and adjusted p-values in an interactive volcano plot
#'
#' @param res Annotated DESeq results table from results_all
#' @param  padj p-value cutoff for labeling points, default 0.05
#' @param log2FoldChange absolute value of log2 fold change cutoff for labeling points, default 2
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

plot_volcano <- function(res, padj = 0.05 , log2FoldChange = 2, ...){
   if(!is_tibble(res)){
      if(is.list(res)){
        message("Plotting the first table in the list")
        res <- res[[1]]
      }else{
         stop("Results should be a tibble")
      }
   }
   x <- filter(res, !is.na(padj))
   ### Grouping column for enableMouseTracking
   x$sig = ifelse( x$padj  < padj | abs(x$log2FoldChange)> log2FoldChange, "Y", "N")
   n <- sum(x$sig == "Y")
   if(n ==0) stop("No points above cutoffs, need to fix hchart to plot this case")
    message("Adding mouseover labels to ", n, " genes (",  round( n/nrow(x)*100, 1), "%)")
    ## use Ensembl ID if gene_name is missing ?
     n <- is.na(x[["gene_name"]]) | x[["gene_name"]] == ""
     if(sum(n)>0 ){
         message("Missing ", sum(n), " gene names, using Ensembl IDs instead")
         x[["gene_name"]][n] <- x[["id"]][n]
      }
    ## center at zero...
     fc <- max(abs(x$log2FoldChange), na.rm=TRUE)
    hchart(x, "scatter", hcaes(log2FoldChange,  -log10(padj), group = sig, value = gene_name), color = 'rgba(0,0,255, 0.3)',
             enableMouseTracking = c(FALSE, TRUE), showInLegend=FALSE) %>%
        hc_tooltip( pointFormat = "{point.value}", headerFormat = "") %>%
         hc_xAxis(title = list(text = "Log2 Fold Change"), gridLineWidth = 1, tickLength = 0, startOnTick = "true", endOnTick = "true" , min= -fc, max=fc) %>%
         hc_yAxis(title = list(text = "-Log10 Adjusted P-value")) %>%
         hc_chart(zoomType = "xy", ...) %>%
         hc_exporting(enabled=TRUE, filename = "volcano")
}
