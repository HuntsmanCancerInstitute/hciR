#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

opts <-  list(
   make_option(c("-s", "--samples"), default="samples.txt",
         help="Tab-delimited file with ids in column 1 matching count column names
                and treatment column for contrasts, default samples.txt"),
   make_option(c("-c", "--counts"), default="counts.txt",
         help="Tab-delimited count matrix file, default counts.txt"),
   make_option(c("-d", "--database"), default="human",
      help="Annotation database, either human, mouse, fly or rat, default human"),
    make_option(c("-f", "--filter"), default=5,
       help="Low count cutoff to filter counts, default 5"),
    make_option(c("-p", "--padj"), default=0.05,
       help="Adjusted p-value cutoff, default 0.05"),
    make_option(c("-t", "--trt"), default="trt",
       help="Name of treatment column in sample table, default trt or treatment"),
    make_option(c("-r", "--relevel"), default="NA",
       help="A comma-separated list to reorder treatments.  By default treatments are sorted
                alphabetically, so use 'C,B,A' to compare C vs A, C vs B and B vs A"),
    make_option(c("-m", "--mouseover"), default="NA",
       help="A comma-separated list of sample column names for tooltips in PCA plot, default is id column")
)

parser <- OptionParser(option_list=opts, description = "
Loads a sample and count table and runs DESeq2 using all possible contrasts.  The results
are saved to an Excel file and html report with sample visualizations.
")
opt <- parse_args(parser)

## use Rmd file in hciR package if missing
if( !file.exists("DESeq.Rmd"))  invisible( file.copy( system.file("Rmd/DESeq.Rmd", package="hciR"), "DESeq.Rmd"))

if( !file.exists( opt$samples) ){
   print_help(parser)
   quit(status=1)
}
#if( !opt$database %in% c("human", "mouse", "fly", "rat")) stop("Database name should be human, mouse, rat, or fly")


suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("hciR"))
suppressPackageStartupMessages(library("DESeq2"))


render("DESeq.Rmd",  params=list( samples = opt$samples, counts = opt$counts,
           database = opt$database, filter = opt$filter, padj = opt$padj,
           trt = opt$trt, relevel = opt$relevel, mouseover = opt$mouseover  ) )
