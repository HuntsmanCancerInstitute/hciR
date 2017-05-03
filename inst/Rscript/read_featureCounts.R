#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("hciR"))

opts <-  list(
   make_option(c("-d", "--directory"), default="NA",
      help="Directory with featuerCounts output files"),
    make_option(c("-o", "--output"), default="counts.txt",
       help="Output file name, default counts.txt")
)


parser <- OptionParser(option_list=opts, description = "
Combine featureCount output files into a single count matix")

 opt <- parse_args(parser)

  if( opt$directory == "NA" ){
     print_help(parser)
     quit(status=1)
  }


counts <- read_featureCounts(opt$directory)
message("Saved ", ncol(counts) - 1, " samples to ", opt$output)
write_tsv(counts, opt$output )
