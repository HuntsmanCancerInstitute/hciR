#' Read MSigDB GMT file into a nested list of databases and genes sets
#'
#' Reads MSigDB GMT files from \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
#'
#' @param gmt a GMT file from MSigDB
#'
#' @return A nested list of databases and genes sets (if only one database in GMT file,
#' then a list of gene sets only)
#'
#' @author Chris Stubben
#'
#' @examples
#'  gmt <- system.file("extdata", "c2.cp.v6.0.symbols.gmt", package = "hciR")
#'  msig <- read_msigdb(gmt)
#'  # a list with 9 pathways dbs
#'  sapply(msig, length)
#'  msig$KEGG[1:3]
#' @export

read_msigdb <- function(gmt){
   x <- readLines(gmt)
   msig <- strsplit(x, "\t")
   # Get source before first _   (KEGG from KEGG_GLYCOLYSIS)
   n <-  gsub("\\_.+", "", sapply(msig, "[", 1))
   ## Drop source from gene set name (GLYCOLYSIS)
   names(msig) <- gsub("^[^_]+_", "",  sapply(msig, "[", 1) )
   # Gene sets start at 3rd element
   msig <- lapply( msig, function(x) x[3: length(x)])
   ## split gene sets by source
   msig <- split(msig, n)
   ## if only 1 source, then list of gene sets only
   if(length(msig) == 1) msig <- msig[[1]]
   msig
}
