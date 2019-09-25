#' Plot reads assigned to biotypes
#'
#' @param x count matrix from featureCounts with `-g gene_biotype`
#' @param n number of top features to plot, default 12
#' @param group plot all features in 12 groups, default FALSE
#' @param stack, plot percent or total reads (percent or normal)
#' @param barwidth bar width, default null
#' @param \dots additional options passed to \code{hc_size}
#'
#' @return A highchart
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  x <- read_tsv("biotypes.txt")
#'  plot_biotypes(x)
#' }
#' @export

plot_biotypes <- function(x, n = 12, group = FALSE, stack = "percent",
 barwidth = NULL, ...){
   names(x)[1] <- "feature"
   if(stack == "percent"){
       ylab <- "Percent"
   }
   else{
      stack <- "normal"
      ylab <- "Total Reads"
   }
   if(!group){
      m1 <- round(rowMeans(as_matrix(x)),0)
      m2 <- tibble::tibble(name = names(m1), total = m1)
      y <- dplyr::top_n(m2, n, total)
      z <- dplyr::filter(x, feature %in% y$name)
      df <- tidyr::gather(z, "sample", "count", -feature)
   }else{
      biotypes <- list(
      antisense = c('antisense', 'antisense_RNA'),
      `protein coding`='protein_coding',
      `IG/TR gene`=c( 'IG_C_gene', 'IG_D_gene',  'IG_J_gene', 'IG_LV_gene', 'IG_V_gene',
         'TR_C_gene', 'TR_J_gene', 'TR_V_gene',  'TR_D_gene'),
      `processed transcript`='processed_transcript',
      `MT rRNA`='Mt_rRNA',
      `MT tRNA`='Mt_tRNA',
      rRNA='rRNA',
      `small RNA`=c('miRNA', 'scRNA', 'snRNA', 'snoRNA', 'sRNA', 'scaRNA',
         'vaultRNA'),
      `other ncRNA`=c('misc_RNA', 'ribozyme', 'non_coding', 'sense_intronic',
         'sense_overlapping', 'known_ncrna'),
      pseudogene=c('pseudogene', 'processed_pseudogene',
      'polymorphic_pseudogene', 'retrotransposed',
      'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene',
      'transcribed_unitary_pseudogene', 'translated_processed_pseudogene',
      'translated_unprocessed_pseudogene', 'unitary_pseudogene',
      'unprocessed_pseudogene', 'IG_pseudogene', 'IG_C_pseudogene',
      'IG_J_pseudogene', 'IG_V_pseudogene', 'TR_V_pseudogene',
      'TR_J_pseudogene', 'Mt_tRNA_pseudogene', 'tRNA_pseudogene',
      'snoRNA_pseudogene', 'snRNA_pseudogene', 'scRNA_pseudogene',
      'rRNA_pseudogene', 'misc_RNA_pseudogene', 'miRNA_pseudogene'),
      lincRNA=c('lincRNA', 'macro_lncRNA', 'bidirectional_promoter_lncRNA'),
      other=c('TEC', 'artifact', '3prime_overlapping_ncRNA',
      'nonsense_mediated_decay', 'non_stop_decay', 'retained_intron',
      'disrupted_domain', 'ambiguous_orf'))
      y <- dplyr::bind_rows(lapply(biotypes,
                 function(x) tibble::tibble(name=x)), .id="feature")
      ## combine with biotypes and sum
      z <- dplyr::inner_join(y, x, by=c(name="feature")) %>%
            dplyr::select(-name) %>%
             dplyr::group_by(feature) %>%
              dplyr::summarise_all(sum)
      df <- tidyr::gather(z, "sample", "count", -feature)
   }

   # PLOT percent or total
   highcharter::hchart(df, "bar",
      highcharter::hcaes(x=sample, y= count, group = feature)) %>%
   highcharter::hc_plotOptions(series = list(stacking = stack),
      bar=list(pointWidth=barwidth)) %>%
   highcharter::hc_yAxis(title=list(text=ylab), reversedStacks=FALSE) %>%
   highcharter::hc_xAxis(title=list(text="")) %>%
   highcharter::hc_size(...) %>%
   highcharter::hc_exporting(enabled=TRUE, filename = "biotypes")

}
