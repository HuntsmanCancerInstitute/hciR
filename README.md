## RNA-seq workflows at [HCI]

This package is intended to work with the [tidyverse] and therefore imports but does not
load [DESeq2] and other [Bioconductor] packages.  Use `devtools` to install the package from GitHub.

```r
library(devtools)
install_github("HuntsmanCancerInstitute/hciR")
```

To learn more about the package, check the [Markdown] files in the [inst/Rmd] directory.
The rendered output is in the [docs] directory, which are best viewed from the Github pages.

1. [pasilla_DESeq.html] - Load pasilla count matrix and samples, run DESeq2, plot PCA and heatmaps.
2. [pasilla_flex.html] - Browse linked MA plot, volcano plot, and result table in a [Flex dashboard] using [Crosstalk].
3. [GSE81784.html] - Load `featureCounts` and run DESeq2 for GSE81784 at HCI.

A few functions are also included in [inst/Rscript] like `read_featureCounts.R` to run on the command line.  See the
[R scripts for RNA-Seq at HCI] file for details.


[HCI]: http://healthcare.utah.edu/huntsmancancerinstitute/
[RNA-seq workflows]: http://www.bioconductor.org/help/workflows/rnaseqGene/
[tidyverse]: http://r4ds.had.co.nz/
[DESeq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[Bioconductor]: https://bioconductor.org/
[inst/Rmd]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/inst/Rmd
[inst/Rscript]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/inst/Rscript
[pasilla_flex.html]: https://huntsmancancerinstitute.github.io/hciR/pasilla_flex.html
[pasilla_DESeq.html]: https://huntsmancancerinstitute.github.io/hciR/pasilla_DESeq.html
[GSE81784.html]: https://huntsmancancerinstitute.github.io/hciR/GSE81784.html
[inst/Rmd]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/inst/Rmd
[docs]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/docs
[Markdown]: http://rmarkdown.rstudio.com/
[Flex dashboard]: http://rmarkdown.rstudio.com/flexdashboard/
[Crosstalk]: https://rstudio.github.io/crosstalk/
[R scripts for RNA-Seq at HCI]: https://huntsmancancerinstitute.github.io/hciR/hciR_scripts.html
