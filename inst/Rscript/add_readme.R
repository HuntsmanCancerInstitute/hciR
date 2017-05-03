#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("rmarkdown"))
suppressPackageStartupMessages(library("hciR"))

   ## using -g returns WARNING: unknown gui ...  so no --genome option

opts <-  list(
   make_option(c("-r", "--run"), default="NA",
         help="Run ID for report title"),
   make_option(c("-a", "--analysis"), default="NA",
         help="Analysis ID for GNomEx file links"),
   make_option(c("-d", "--database"), default="human",
      help="Reference database, either human, mouse, fly or rat, default human"),
    make_option(c("-l", "--length"), default=50,
       help="Read length, default 50")
)

parser <- OptionParser(option_list=opts, description = "
Creates a README.html file for the  STAR RNA-seq workflow at HCI using defaults
from the cmd.txt file in setup_jobs.R.
")
 opt <- parse_args(parser)

## use Rmd file in hciR package if missing
if( !file.exists("README.Rmd"))  invisible(file.copy( system.file("Rmd/README.Rmd", package="hciR"), "README.Rmd"))

if( "NA"  %in% c(opt$run, opt$analysis  )){
   print_help(parser)
   quit(status=1)
}
if( !opt$database %in% c("human", "mouse", "fly", "rat")) stop("Database name should be human, mouse, rat, or fly")

render("README.Rmd", quiet=TRUE, params=list( run = opt$run, analysis = opt$analysis,
    database = opt$database, length = opt$length ) )

message("Added README.html")
