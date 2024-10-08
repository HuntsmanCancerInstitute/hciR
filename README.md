## RNA-Seq workflows at HCI

The `hciR` package works best with [tidyverse] packages (readr, dplyr, tibble,
etc.) and simplifies the code in a differential expression analysis.  The
package includes functions to run [DESeq2] using sample and count tibbles as
input, get annotated DESeq results for all pairwise comparisons and create
interactive plots and other visualizations.

***NOTE***: Starting in hciR version 1.5, the control group should be listed first
and then pairwise comparisons are selected in reverse order, so A, B, C levels
will return results for C vs. B, C vs. A, and B vs. A from `results_all`.

Use `devtools` to install `hciR` and the [hciRdata] package
with Ensembl annotations.

```r
library(devtools)
install_github("HuntsmanCancerInstitute/hciR")
install_github("HuntsmanCancerInstitute/hciRdata")
```


The basic workflow for a mouse experiment with three groups in a `trt` column is
listed below.

```r
library(hciR)
samples <- read_tsv("samples.txt")
samples$trt <- factor(samples$trt, levels = c("WT", "OE", "KO"))
counts <- read_tsv("counts.txt")
counts <- filter_counts(counts, n = 5)
dds <- deseq_from_tibble(counts, samples, design = ~ trt)
rld <- r_log(dds)
plot_pca(rld, "trt", tooltip= c("id", "name"))
plot_dist(rld, "trt", na_col="white")
library(hciRdata)
res <- results_all(dds, mouse104)
plot_volcano(res[[1]])
x <- top_counts(res[[1]], rld, top=40)
plot_genes(x, "trt", scale ="row", annotation_names_col=FALSE)
write_deseq(res, dds, rld, mouse104)
```


Check the vignettes directory to learn more about the package.  The [Pasilla] vignette
runs through an analysis with a single contrast and [Liver] includes an interaction
model and gene set enrichment.  The [Ensembl] file has details on loading annotations.

The [hciRscripts] package wraps functions like `read_featureCounts` to run
on the command line.  See the [hciR scripts] file for more details.

[hciRdata]: https://github.com/HuntsmanCancerInstitute/hciRdata
[DESeq2]: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
[tidyverse]: http://tidyverse.org/
[Ensembl]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/vignettes/Ensembl.md
[Liver]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/vignettes/Liver.md
[Pasilla]: https://github.com/HuntsmanCancerInstitute/hciR/blob/master/vignettes/Pasilla.md
[hciR scripts]: https://huntsmancancerinstitute.github.io/hciRscripts/hciR_scripts.html
[hciRscripts]: https://github.com/HuntsmanCancerInstitute/hciRscripts
