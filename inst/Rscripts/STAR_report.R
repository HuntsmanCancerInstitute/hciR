#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("hciR"))

opts <-  list(
   make_option(c("-g", "--genome"),  default="human",
      help="Genome name, either human, mouse or rat"),
    make_option(c("-o", "--output"), default="README.html",
       help="Output html file name"),
    make_option(c("-r", "--run"),
       help="Run ID for report title"),
    make_option(c("-a", "--analysis"),
       help="Analysis ID for GNomEx file links")
)
opt = parse_args(OptionParser(option_list=opts))

## path to Rmd files
x <- system.file("Rmd", package="hciR")
if(!opt$genome %in% c("human", "mouse", "rat")) stop("Genome name should be human, mouse or rat")
input <- paste0(x, "/STAR_", opt$genome, ".Rmd")

render(input, output_file=opt$output, output_dir=".", params=list( run = opt$run, analysis = opt$analysis ) )
