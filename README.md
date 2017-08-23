## RNA-seq workflows at HCI

This package includes functions to  [RNA-seq workflows].

basic workflow

intended to work with the [tidyverse].

1.  read featureCounts, RSEM or other output files into a count tibble
2.  read samples into a tibble  (read_tsv or other readr command)
3.  run deseq using sample and count tibbles  (deseq)
4.  plot interactive pca using highcharter (plot_pca)
5.  plot sample distances with annotation bars (plot_dist)
6.  results_all
7.  join top genes and normalized ocunts (rlog or vst)
8.  plot genes pheatmaps with annotation bars



1 functions for [RNA-seq workflows]    (read featureCounts or RSEM output files)
2.  intended to work with the [tidyverse].    (run deseq using sample and count tibbles)
3. results_all    (run many contrasts and add annotations to result tables)
3.  create interactive graphics     ( plot_pca using highcharter)
4.  heatmaps (plot_dist,  plot_genes)

2.  simiplify the code blocks in R Markdown reprots
to generate [R Markdown] reprots
3
It imports but does not load [DESeq2] and other [Bioconductor] packages.

Use `devtools` to install the package from GitHub.

```r
library(devtools)
install_github("HuntsmanCancerInstitute/hciR")
```

A few [R Markdown] sample files are in the [inst/Rmd] directory and rendered
output in the [docs] directory, which are best viewed from the Github pages.

1. [pasilla_DESeq.html] - Load pasilla count matrix and samples, run DESeq2, plot PCA and heatmaps.
2. [pasilla_flex.html] - Browse linked MA plot, volcano plot, and result table in a [Flex dashboard] using [Crosstalk].
3. [GSE81784.html] - Load `featureCounts` and run DESeq2 for GSE81784 at HCI.



[RNA-seq workflows]: http://www.bioconductor.org/help/workflows/rnaseqGene/
[tidyverse]: http://r4ds.had.co.nz/
[DESeq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[Bioconductor]: https://bioconductor.org/
[inst/Rmd]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/inst/Rmd
[pasilla_flex.html]: https://huntsmancancerinstitute.github.io/hciR/pasilla_flex.html
[pasilla_DESeq.html]: https://huntsmancancerinstitute.github.io/hciR/pasilla_DESeq.html
[GSE81784.html]: https://huntsmancancerinstitute.github.io/hciR/GSE81784.html
[inst/Rmd]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/inst/Rmd
[docs]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/docs
[R Markdown]: http://rmarkdown.rstudio.com/
[Flex dashboard]: http://rmarkdown.rstudio.com/flexdashboard/
[Crosstalk]: https://rstudio.github.io/crosstalk/
