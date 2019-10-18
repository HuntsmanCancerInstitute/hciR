## RNA-Seq workflows at [HCI]

Use `devtools` to install the package from GitHub.

```r
library(devtools)
install_github("HuntsmanCancerInstitute/hciR")
```

The `hciR` package works best with [tidyverse] packages (readr, dplyr, tibble,
etc.) and simplifies the code in a differential expression analysis.  The
package includes functions to run DESeq2 using sample and count tibbles as
input, get annotated DESeq results for all pairwise combination and create
interactive plots and other visualizations.

The basic workflow for a mouse experiment with three groups in a `trt` column is
listed below.

```r
library(hciR)
samples <- read_tsv("samples.txt")
counts <- read_tsv("counts.txt")
counts <- filter_counts(counts, n = 5)
dds <- deseq_from_tibble(counts, samples, design = ~ trt)
rld <- r_log(dds)
plot_pca(rld,  "trt", tooltip= c("id", "name"))
plot_dist(rld, "trt", na_col="white")
library(hciRdata)
res <- results_all(dds, mouse98)
plot_volcano(res[[1]])
x <- top_counts(res[[1]], rld, top=40)
plot_genes(x, "trt", scale ="row", annotation_names_col=FALSE)
write_deseq(res, dds, rld, mouse98)
```


Check the vignettes directory to learn more about the package.  The [Pasilla] vignette
runs through an analysis with a single contrast and [Liver] includes an interaction
model and gene set enrichment.  The [Ensembl] file has details on loading annotations.

The [hciRscripts] package wraps functions like `read_featureCounts` to run
on the command line.  See the [hciR scripts] file for more details.


[tidyverse]: http://tidyverse.org/
[Ensembl]: https://github.com/HuntsmanCancerInstitute/hciR/vignettes/Ensembl.md
[Pasilla]: https://github.com/HuntsmanCancerInstitute/hciR/vignettes/Pasilla.md
[hciR  scripts]: https://huntsmancancerinstitute.github.io/hciRscripts/hciR_scripts.html
[hciRscripts]: https://github.com/HuntsmanCancerInstitute/hciRscripts
