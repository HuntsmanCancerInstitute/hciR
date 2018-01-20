#' Read biomart annotations
#'
#' Read Ensembl annotations using \code{biomaRt} library
#'
#' @param dataset The first letter of genus and full species name like btaurus, ecaballus, sscrofa from \code{listEnsembl}.
#' A few common names are accepted for human, mouse, rat, fruitfly, yeast and zebrafish.
#' @param attributes vector of column names to pass to \code{getBM}, default ensembl_gene_id,
#'  external_gene_name, gene_biotype, chromosome_name, start_position, end_position, strand,
#'  description, transcript_count, and entrezgene
#' @param host URL from \code{listEnsemblArchives()} to download previous versions
#' @param list return a list of either datasets, attributes or filters only.
#' @param \dots additional options passed to \code{getBM} or \code{listAttributes}
#'
#' @note The version option in \code{useEnsembl} no longer works, so this was replaced with host on Jan 2, 2018.
#  Many attributes like entrezgene have a 0 to many relationship with ensembl_gene_id
#' causing duplicate ensembl ids to be added.  The default returns a 10 column table where
#' entrezgene is grouped into a comma-separated list so ensembl id is unique.
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' x1 <- read_biomart( list= "datasets")
#' x1
#' # Download latest version
#'  mmu <- read_biomart("mouse")
#'  # Download mouse genes and human homologs
#'  x1 <- read_biomart("mouse", list= "attributes")
#'  dplyr::count(x1, page)
#'  filter(x1, grepl("hsapiens_homolog", name))
#'  mmu_homologs <- read_biomart("mouse",  attributes = c("ensembl_gene_id",
#'    "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name",
#'    "hsapiens_homolog_perc_id"))
#' # Download version 87
#' listEnsemblArchives()
#  mmu87 <- read_biomart("mouse", host = "http://Dec2016.archive.ensembl.org")
#'  # Human genes with SignalP
#'  x2 <- read_biomart("human", list= "filters")
#'  filter(x2, grepl("signal", name))
#'  hsa_signalp  <- read_biomart("human", attributes = c("ensembl_transcript_id",
#'        "ensembl_gene_id",  "external_gene_name", "signalp_start",  "signalp_end"),
#'        filter="with_signalp", values=TRUE)
#' }
#' @export


read_biomart <- function( dataset="human", attributes, host="www.ensembl.org",
   list = NULL, ...){
   common <- c(human = "hsapiens",
               mouse = "mmusculus",
               rat = "rnorvegicus",
               zebrafish = "drerio",
               fly = "dmelanogaster",
               fruitfly = "dmelanogaster",
               pig = "sscrofa",
               yeast = "scerevisiae")
   if( tolower(dataset) %in% names(common))  dataset <- common[[tolower(dataset)]]
   if( !grepl("gene_ensembl$", dataset) )  dataset <- paste0(dataset, "_gene_ensembl")

  ## LIST
   if(!is.null(list)){
      if(tolower(list) == "datasets"){
         ensembl <- biomaRt::useEnsembl(biomart="ensembl", host=host)
         bm <- biomaRt::listDatasets(ensembl)
         ## remove AsIs class
         for(i in 1:3) class(bm[,i]) <- "character"
         message("Downloaded ", nrow(bm), " datasets")
      }else{
         ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset, host=host)
         if(tolower(list) == "filters"){
             bm <- biomaRt::listFilters(ensembl, ...)
             message("Downloaded ", nrow(bm), " filters")
         }else{
             bm <- biomaRt::listAttributes(ensembl, ...)
             message("Downloaded ", nrow(bm), " attributes")
         }
      }
   }else{
      ## SEARCH
      ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset, host=host)
      # default search
      if(missing(attributes)){
          bm <- biomaRt::getBM(attributes=c('ensembl_gene_id','external_gene_name', 'gene_biotype',
                  'chromosome_name', 'start_position', 'end_position','strand', 'description',
                  'transcript_count', 'entrezgene'), mart = ensembl, ...)
           # replace long names like ensembl_gene_id
           names(bm)[1:6] <- c("id", "gene_name", "biotype", "chromosome", "start", "end")
           n <-   length(unique(bm$id))
           message("Downloaded ", nrow(bm), " results, grouping into ", n, " rows with a unique Ensembl ID" )
           ## group by first 9 columns and paste entrezgene into comma-separated list so ensembl id is unique
           bm <- dplyr::group_by_( bm , .dots= colnames(bm)[-ncol(bm)]) %>%
             dplyr::summarize( entrez_gene = paste(entrezgene, collapse=","))
           # drop source from description  [Source:MGI Symbol;Acc:MGI:102478]
         bm$description <- gsub(" \\[.*\\]$", "" , bm$description)
      }else{
         bm <- biomaRt::getBM(attributes= attributes, mart = ensembl, ...)
         message("Downloaded ", nrow(bm), " features")
      }
   }
   # will also drop the grouped_df class
   bm <- tibble::as_data_frame(bm)
   attr(bm, "downloaded") <- Sys.Date()
   attr(bm, "host") <- host
   bm
}
