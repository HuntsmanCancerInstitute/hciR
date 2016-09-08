#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rmarkdown"))

opts <-  list(
    make_option(c("-o", "--output"), default="out",
       help="CollectRnaSeqMetrics output directory, default out/"),
    make_option(c("-p", "--pattern"),   default="txt$",
       help="Load file names matching regular expression, default txt$"),
    make_option(c("-r", "--recursive"), default=TRUE, action="store_true",
       help="Recurse subdirectories in output directory, default TRUE"),
    make_option(c("-i", "--input"),     default="RnaSeqMetrics.Rmd", metavar="Rmd",
       help="R markdown input file, default RnaSeqMetrics.Rmd")
)

opt = parse_args(OptionParser(option_list=opts))

render(opt$input, params=list( output = opt$output, pattern = opt$pattern, recursive=opt$recursive ) )
