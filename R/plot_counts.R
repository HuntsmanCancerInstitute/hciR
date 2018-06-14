#' Plot counts for a single gene
#'
#' @param rld rlog or other counts in DESeqTransform object
#' @param gene gene id, should match rownames in rld
#' @param intgroups one or two names in colData(rld) to use for grouping
#' @param ylab y-axis label
#' @param title ggplot title, use NULL for no title
#'
#' @note See \code{link{plot_interactions}} to plot many genes.
#'
#' @return A ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  DESeq2::plotCounts(dds, "ENSG00000075624",  "time")
#'  plot_counts(rld, "ENSG00000075624",  "time")
#'  # interaction plot
#'  plot_counts(rld, "ENSG00000075624",  c("time", "trt"), title = "ACTB")
#' }
#' @export

plot_counts <- function(rld, gene, intgroups, ylab="count", title){
   if(!class(rld) == "DESeqTransform") stop("rld shoud be a DESeqTransform")
   if(length(gene) > 1) gene <- gene[1]
   if(!gene %in% rownames(rld)) stop(gene, " is not found in rld, check rownames(rld) for valid input")
   if(missing(title)) title <- gene
   x <- data.frame( colData(rld)[, intgroups, drop=FALSE])
   x$gene <- as.vector( assay( rld[gene,]))
   if(length(intgroups) == 1){
       ggplot(x, aes_string(x=intgroups, y="gene")) + geom_point() +
        stat_summary(aes(y=gene, group=1), fun.y="mean", geom="line") +
         ylab(ylab) + ggtitle(title)
    }else{
       ggplot(x, aes_string(x=intgroups[1], y="gene", color=intgroups[2])) + geom_point() +
        stat_summary(aes_string(y="gene", group=intgroups[2]), fun.y="mean", geom="line") +
         ylab(ylab) + ggtitle(title)
   }
}
