#' Gene set enrichment using Fisher's test
#'
#' Gene set enrichment using \code{fisher.test}.
 #'
#' @param res  A list of DESeq results
#' @param gsets Gene sets
#' @param padj Adjusted p-value cutoff for significant genes
#'
#' @return a tibble with gene sets intersections and pvalue from \code{fisher.test}
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
    y <- top_n(x, top, -pvalue)
    z <- mutate(y, `Up-regulated` = up/total, `Down-regulated` =down/total,
       signif = signif/total) %>%
     dplyr::select(c(1:3,8:9)) %>%
       gather("Match", "Percent", -c(term,total, signif))  %>% mutate(n=row_number())
   z$Match <- factor(z$Match,  c( "Up-regulated", "Down-regulated"))
   ## order by percent

   ggplot(z, aes(x=reorder(term, Percent), y=Percent, fill=Match )) +
   geom_bar(stat="identity") +
   geom_text(aes(label=total, y= signif), size=2.5, hjust=-0.1) +
     coord_flip() + scale_fill_manual(values=c( "red", "green")) + xlab("") +
     scale_y_continuous(labels=scales::percent, expand = expand_scale(mult = c(0, 0.07))) +
      ylab(ylab) +
      theme_classic() + guides(fill = guide_legend(reverse=TRUE, title=NULL))
}
