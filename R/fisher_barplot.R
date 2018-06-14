#' Gene set enrichment using Fisher's test
#'
#' Gene set enrichment using \code{fisher.test}.
 #'
#' @param x  output from \code{\link{fisher_all}}
#' @param top number of top genes
#' @param ylab y-axis label
#'
#' @return a ggplot
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' data(msig)
#' fisher_all(res, msig$KEGG)
#' }
#' @export

fisher_barplot <- function(x, top =10, ylab = "Significant genes in pathway"){
    y <- dplyr::top_n(x, top, -pvalue)
    z <- dplyr::mutate(y, `Up-regulated` = up/total, `Down-regulated` =down/total,
       signif = signif/total) %>%
     dplyr::select(c(1:3,8:9)) %>%
       tidyr::gather("Match", "Percent", -c(term,total, signif))  %>%
        dplyr::mutate(n=row_number())
   z$Match <- factor(z$Match,  c( "Up-regulated", "Down-regulated"))
   ## order by percent

   ggplot2::ggplot(z, ggplot2::aes(x=stats::reorder(term, Percent), y=Percent, fill=Match )) +
   ggplot2::geom_bar(stat="identity") +
   ggplot2::geom_text( ggplot2::aes(label=total, y= signif), size=2.5, hjust=-0.1) +
     ggplot2::coord_flip() + ggplot2::scale_fill_manual(values=c( "red", "green")) +
     ggplot2::scale_y_continuous(labels=scales::percent, expand = ggplot2::expand_scale(mult = c(0, 0.07))) +
      ggplot2::xlab("") + ggplot2::ylab(ylab) +
      ggplot2::theme_classic() + ggplot2::guides(fill = guide_legend(reverse=TRUE, title=NULL))
}
