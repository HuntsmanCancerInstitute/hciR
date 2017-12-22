#' Read biomart annotations
#'
#' Read Ensembl annotations using \code{biomaRt} library
#'
#' @param dataset The first letter of genus and full species name like btaurus, ecaballus, sscrofa from \code{listEnsembl}.
#' A few common names are accepted for human, mouse, rat, fruitfly, yeast and zebrafish.
#' @param attributes vector of column names to pass to \code{getBM}, default ensembl_gene_id,
#'  external_gene_name, gene_biotype, chromosome_name, start_position, end_position, strand,
#'  description, transcript_count, and entrezgene
#' @param version Ensembl version number for \code{useEnsembl}. Default is current version.
#' @param list return a list of either datasets, attributes or filters only.
#' @param \dots additional options passed to \code{getBM} or \code{listAttributes}
#'
#' @note Many attributes like entrezgene have a 0 to many relationship with ensembl_gene_id
#' causing duplicate ensembl ids to be added.  The default returns a 10 column table where
#' entrezgene is grouped into a comma-separated list so ensembl id is unique.
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'  mmu <- read_biomart("mouse")
#'  # Mouse and Human homologs
#'  x1 <- read_biomart("mouse", list= "attributes")
#'  count(x1, page)
#'  filter(x1, grepl("hsapiens_homolog", name)
#'  mmu_homologs <- read_biomart("mouse",  attributes = c("ensembl_gene_id",
#'    "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name",
#'    "hsapiens_homolog_perc_id"))
#'  # genes with SignalP
#'  x2 <- read_biomart("human", list= "filters")
#'  filter(x2, grepl("signal", name))
#'    mmu_signalp  <- read_biomart("human", attributes = c("ensembl_transcript_id",
#'         "ensembl_gene_id",  "external_gene_name", "signalp_start",  "signalp_end"),
#'          filter="with_signalp", values=TRUE)
#' }
#' @export

read_biomart <- function( dataset="human" , attributes, version = NULL, list = NULL, ...){
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

   ## SAVE version as attribute
   if(is.null(version)){
      x <- biomaRt::listEnsembl()
      version2 <- x$version[x$biomart=="ensembl"]
      version2 <- gsub("Ensembl Genes ", "", version2)
   }
   message("Ensembl Release ", version2)

  ## LIST
   if(!is.null(list)){
      if(tolower(list) == "datasets"){
         ensembl <- biomaRt::useEnsembl(biomart="ensembl")
         bm <- biomaRt::listDatasets(ensembl)
         message("Downloaded ", nrow(bm), " datasets")
      }else{
         ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset, version = version)
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
      ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset, version = version)
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
   attr(bm, "version") <- version2
   bm
}
