## RNA-seq workflows at HCI

Use `devtools` to install the package from GitHub.

```r
library(devtools)
install_github("HuntsmanCancerInstitute/hciR")
```

The package includes functions to read output files from [featureCounts] and [RSEM],
display sample visualization using interactive [highcharts] (see [plot_pca] for the highcharts version of `plotPCA` in [DESeq2]),
run all possible contrasts and add [Biomart] annotations to results tables, create heatmaps using [d3heatmap]
or [pheatmap] and other RNA-seq workflow tasks.   The package is intended to work with the [tidyverse] and
imports but does not load [DESeq2] and other [Bioconductor] packages.

A few [R Markdown] sample files are in the [inst/Rmd] directory and rendered output in the [docs] directory,
which are best viewed from the Github pages.

1. [pasilla_DESeq.html] - Load pasilla count matrix and samples, run DESeq2, plot PCA and heatmaps.
2. [pasilla_flex.html] - Browse linked MA plot, volcano plot, and result table in a [Flex dashboard] using [Crosstalk].
3. [GSE81784.html] - Load `featureCounts` and run DESeq2 for GSE81784 at HCI.


[featureCounts]: http://bioinf.wehi.edu.au/featureCounts/
[RSEM]: http://deweylab.github.io/RSEM/
[Biomart]: http://uswest.ensembl.org/biomart/martview
[d3heatmap]: http://www.htmlwidgets.org/showcase_d3heatmap.html
[pheatmap]: https://cran.r-project.org/web/packages/pheatmap/index.html
[highcharts]: http://jkunst.com/highcharter/
[plot_pca]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/R/plot_pca.R

[DESeq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[Bioconductor]: https://bioconductor.org/
[inst/Rmd]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/inst/Rmd
[pasilla_flex.html]: https://huntsmancancerinstitute.github.io/hciR/pasilla_flex.html
[pasilla_DESeq.html]: https://huntsmancancerinstitute.github.io/hciR/pasilla_DESeq.html
[GSE81784.html]: https://huntsmancancerinstitute.github.io/hciR/GSE81784.html
[inst/Rmd]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/inst/Rmd
[docs]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/docs
[tidyverse]: http://r4ds.had.co.nz/
[R Markdown]: http://rmarkdown.rstudio.com/
[Flex dashboard]: http://rmarkdown.rstudio.com/flexdashboard/
[Crosstalk]: https://rstudio.github.io/crosstalk/
