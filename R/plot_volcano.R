#' Volcano plot
#'
#' Plot fold changes and adjusted p-values in an interactive volcano plot or
#' ggplot
#'
#' For ggplot, the results should not be sorted by p-value or fold change to
#' avoid stacking close overlapping points.  Labels are added using
#' \code{ggrepel}, so avoid labeling too many points (200 is the limit).
#'
#' @param res Annotated DESeq results table from results_all.
#' @param pvalue_cutoff either a single p-value cutoff on the y-axis or two
#' p-values to label down- and up-regulated genes. The default is 1.3
#' (corresponding to padj = 0.05) for highcharts and no labels for ggplot.
#' @param foldchange_cutoff either the absolute value of the log2 fold change
#' cutoff or negative and positive fold changes to label genes, default is 2 for
#' highcharts and no labels for ggplot.
#' @param max_pvalue y-axis limit, maximum value on a -10 log10 y-axis scale,
#' the default is 200 (padj < 1e-200), so genes below this cutoff are assigned
#' the maximum p-value.
#' @param radius point size, default 3
#' @param ggplot plot ggplot version, default TRUE
#' @param palette RColorBrewer palette name, vector of colors, or "RdGn" for
#' ggplot
#' @param \dots other options like width passed to \code{hc_chart}
#'
#' @return A highchart or ggplot.
#'
#' @author Chris Stubben
#'
#' @examples
#' plot_volcano(pasilla$results, pvalue=c(35,25), foldchange=2.5)
#' @export

plot_volcano<- function(res, pvalue_cutoff, foldchange_cutoff, max_pvalue = 200,
 radius=3, ggplot=TRUE, palette="RdBu", ...){
   if(!tibble::is_tibble(res)){
      if(is.list(res)){
        message("Plotting the first table in the list")
        res <- res[[1]]
      }else{
         stop("Results should be a tibble")
      }
   }
   ## limma top_tibble?
   n <- match(c("logFC", "adj.P.Val"), names(res))
   if( !any(is.na(n))){
      names(res)[n] <- c("log2FoldChange", "padj")
      # Symbol, Gene.Symbol, others
      n <- which( names(res) %in% c("Symbol", "Gene.Symbol"))
      if(length(n)==1){
         names(res)[n] <- "gene_name"
         res$gene_name <- gsub(", .*", "", res$gene_name)
      }
   }
   x <- dplyr::filter(res, !is.na(padj))
   ## center at zero...
   fc <- max(abs(x$log2FoldChange), na.rm=TRUE)

   ## fix points with zero
   n <- x$padj ==0
   if(sum(n)>0){
      message( "Setting ", sum(n), " zero p-values to 1e-", max_pvalue)
      x$padj[n] <-  1/10^max_pvalue
   }

   ## fix points with really low -p-values?
   n <- x$padj < 1/10^max_pvalue
   if(sum(n)>0){
      message( "Setting ", sum(n), " low p-values to 1e-", max_pvalue)
      x$padj[ n ] <-  1/10^max_pvalue
   }
   ## ggplot ------------------------------
   if(ggplot){
      clrs <- palette
      if(length(clrs)==1) clrs <- palette255(clrs, ramp=FALSE)
      p1 <- ggplot2::ggplot(data=x,
              ggplot2::aes(x=log2FoldChange, y= -log10(padj))) +
         ggplot2::geom_point(ggplot2::aes(fill = log2FoldChange),
             color="gray20", shape = 21, size=radius) +
         ggplot2::xlim( -fc, fc) + ggplot2::theme_light() +
         ggplot2::xlab("Log2 Fold Change") +
         ggplot2::ylab("-Log10 Adjusted P-value") +
         ggplot2::scale_fill_gradientn(colors=clrs, limits=c(-fc, fc),
             guide=FALSE)
		if( missing(pvalue_cutoff) & missing(foldchange_cutoff)){
			p1
		}else{
		   ## add labels, pvalue, fc or BOTH
		   ## avoid ggrepel warning if missing name
		   y <- filter(x, !is.na(gene_name))
		   y1 <- NULL
		   y2 <- NULL
           if(!missing(pvalue_cutoff)){
             if(length(pvalue_cutoff) == 2){
                 y1 <- filter(y,
					   padj < 1/10^pvalue_cutoff[1] & log2FoldChange < 0 |
					   padj < 1/10^pvalue_cutoff[2] & log2FoldChange > 0 )
            }else{
                 y1 <- filter(y, padj < 1/10^pvalue_cutoff )
            }
		  }
		  ##
		  if(!missing(foldchange_cutoff)){
            if(length(foldchange_cutoff)==2){
                 y2 <- filter(y, log2FoldChange < foldchange_cutoff[1] |
                                 log2FoldChange > foldchange_cutoff[2])
            }else{
                 y2 <- filter(y, abs(log2FoldChange) > foldchange_cutoff)
            }
		  }
		  # combine cutoffs
		  y <- dplyr::bind_rows(y1, y2) %>% unique()
          if(nrow(y) > 0){
            if(nrow(y) > 200){
               message("Too many points for ggrepel to label (", nrow(y),
                   "), check pvalue and fold change cutoffs")
               p1
            }else{
                p1 + ggrepel::geom_text_repel(data=y, ggplot2::aes(label=
                        gene_name), cex=3, box.padding=.1, point.padding=.1)
            }
		 }
	  }
   # Highcharts ------------------------------
   }else{
      ### Grouping column for enableMouseTracking
      if(missing(pvalue_cutoff)) pvalue_cutoff <-  -log10(0.05)
      if(missing(foldchange_cutoff))  foldchange_cutoff <- 2
      if(length(foldchange_cutoff)==2){
         x$sig <- ifelse( x$padj  < 1/10^pvalue_cutoff |
               x$log2FoldChange < foldchange_cutoff[1] |
               x$log2FoldChange > foldchange_cutoff[2],  "Y", "N")
      }else{
         x$sig <- ifelse( x$padj  < 1/10^pvalue_cutoff |
              abs(x$log2FoldChange) > foldchange_cutoff, "Y", "N")
      }
      n <- sum(x$sig == "Y")
      if(n ==0){
         # message("No points above cutoffs")
         x$sig[which.min(x$padj)] <- "Y"  #need one Y for enableMousetracking
         n <- 1
      }
      message("Adding mouseover labels to ", n, " genes (",
                 round( n/nrow(x)*100, 1), "%)")
      ## use Ensembl ID if gene_name is missing ?
      n <- x$sig=="Y" & (is.na(x[["gene_name"]]) | x[["gene_name"]] == "")
      if(sum(n)>0 ){
         message("Missing ", sum(n), " gene names, using Ensembl ID instead")
         x[["gene_name"]][n] <- x[["id"]][n]
      }
      highcharter::hchart(x, "scatter", highcharter::hcaes(log2FoldChange,
         -log10(padj), group=sig, value=gene_name), color= 'rgba(0,0,255,0.3)',
          enableMouseTracking = c(FALSE, TRUE), showInLegend=FALSE,
          marker = list(radius = radius, lineColor="blue")) %>%
      highcharter::hc_tooltip(pointFormat="{point.value}", headerFormat="") %>%
      highcharter::hc_xAxis(title = list(text = "Log2 Fold Change"),
         gridLineWidth = 1, tickLength = 0, startOnTick = "true",
         endOnTick = "true" , min= -fc, max=fc) %>%
      highcharter::hc_yAxis(title = list(text = "-Log10 Adjusted P-value")) %>%
      highcharter::hc_chart(zoomType = "xy", ...) %>%
      highcharter::hc_exporting(enabled=TRUE, filename = "volcano")
   }
}
