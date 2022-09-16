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
#'  plot_counts(pasilla$rlog, "FBgn0033635",  c("type","condition"), title = "Prip")
#' @export

plot_counts <- function(rld, gene, intgroups, ylab="count", title){
   if(!inherits(rld, "DESeqTransform")) stop("rld shoud be a DESeqTransform")
   if(length(gene) > 1) gene <- gene[1]
   if(!gene %in% rownames(rld)){
       stop(gene, " is not found in rld, check rownames(rld) for valid input")
   }
   if(missing(title)) title <- gene
   x <- data.frame( SummarizedExperiment::colData(rld)[, intgroups, drop=FALSE])
   x$gene <- as.vector( SummarizedExperiment::assay( rld[gene,]))
   if(length(intgroups) == 1){
      ggplot2::ggplot(x, ggplot2::aes_string(x=intgroups, y="gene")) +
       ggplot2::geom_point() +
        ggplot2::stat_summary(ggplot2::aes(y=gene, group=1), fun.y="mean",
          geom="line") + ggplot2::ylab(ylab) +
           ggplot2::ggtitle(title)
    }else{
       ggplot2::ggplot(x, ggplot2::aes_string(x=intgroups[1], y="gene",
        color=intgroups[2])) + ggplot2::geom_point() +
         ggplot2::stat_summary(ggplot2::aes_string(y="gene",
                   group=intgroups[2]), fun.y="mean", geom="line") +
          ggplot2::ylab(ylab) + ggplot2::ggtitle(title)
   }
}
