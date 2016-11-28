#' Highchart version of plotPCA in DESeq2
#'
#' Highchart version of sample PCA plot for transformed data
#'
#' @param object a DESeqTransform object
#' @param intgroup a character vector of names in colData(x) for grouping, default 'condition'
#' @param tooltip a character vector of names in colData(x) for tooltip display,
#'       default displays the object column names'
#' @param ntop number of top variable genes to use for principal components
#'
#' @return A highchart
#'
#' @author Chris Stubben
#'
#' @examples
#' dds <- makeExampleDESeqDataSet(betaSD=1)
#' rld <- rlog(dds)
#' plotPCA(rld)
#' hc_plotPCA(rld)
#' colData(rld)$sample <- 1:12
#' colData(rld)$gender <- sample(c("M", "F"), 12, TRUE)
#' hc_plotPCA(rld, tooltip=c("sample", "gender"))
#' }
#' @export

hc_plotPCA <- function(object, intgroup="condition", tooltip, ntop = 500){
   n  <- apply(assay(rld), 1, var)
   x  <-  head(assay(rld)[ order(n, decreasing=TRUE),], ntop)
   pca <- prcomp(t( x ))
   percentVar <- round(pca$sdev^2/sum(pca$sdev^2)*100, 1)

   if(!all(intgroup %in% names(colData(object)))) stop("intgroup should match columns of colData(object)")
   group <-  apply( as.data.frame(colData(object)[, intgroup, drop=FALSE]), 1, paste, collapse=": ")
   d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], INTGRP=group, COLNAMES = colnames(object), colData(object) )

   # if tooltip is missing use column names
   if(missing(tooltip)){
      tooltipJS <- "this.point.COLNAMES"
    }else{
       if(!all(tooltip %in% names(colData(object)))) stop("tooltip should match columns of colData(object)")
       ## if tooltip = ID, patient, gender  then tooltipJS =
      ##  "'ID: ' + this.point.ID + '<br>patient: ' + this.point.patient + '<br>gender: ' + this.point.gender
      tooltipJS <-  paste0("'", paste( tooltip, ": ' + this.point.", tooltip, sep="", collapse = " + '<br>"))
   }

   hchart(d, "scatter", x=PC1, y=PC2, group=INTGRP ) %>%
    hc_tooltip(formatter = JS( paste0("function(){ return (", tooltipJS, ")}"))) %>%
     hc_xAxis(title = list(text = paste0("PC1: ", percentVar[1], "% variance")),
          gridLineWidth=1, tickLength=0, startOnTick="true", endOnTick="true") %>%
      hc_yAxis(title = list(text = paste0("PC2: ", percentVar[2], "% variance"))) %>%
       hc_chart(zoomType = "xy")  %>%
        hc_exporting(enabled=TRUE, filename = "pca")
}
