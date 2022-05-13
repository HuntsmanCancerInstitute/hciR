#' Read biomart annotations
#'
#' Read Ensembl annotations using \code{biomaRt} library
#'
#' @param dataset The first letter of genus and full species name like btaurus,
#' ecaballus, sscrofa from \code{listEnsembl}. A few common names are accepted
#' for human, mouse, rat, fruitfly, yeast and zebrafish.
#' @param attributes vector of column names to pass to \code{getBM}, default
#' ensembl_gene_id, external_gene_name, gene_biotype, chromosome_name,
#' start_position, end_position, strand, description and transcript_count
#' @param host host for connection, default https://www.ensembl.org
#' @param version Ensembl version for previous releases
#' @param patch Keep features on patches starting with CHR_, default FALSE
#' @param list return a list of either datasets, attributes or filters only.
#' @param \dots additional options passed to \code{getBM} or
#' \code{listAttributes}
#'
#' @note Many attributes like entrezgene have a 0 to many relationship with
#' ensembl_gene_id causing duplicate ensembl ids to be added.
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
#'  # Human genes with SignalP
#'  x2 <- read_biomart("human", list= "filters")
#'  filter(x2, grepl("signal", name))
#'  hsa_signalp <- read_biomart("human", attributes = c("ensembl_transcript_id",
#'   "ensembl_gene_id",  "external_gene_name", "signalp_start",  "signalp_end"),
#'   filter="with_signalp", values=TRUE)
#' }
#' @export
read_biomart <- function(dataset="human", attributes, host="https://www.ensembl.org",
   version=NULL, patch=FALSE, list=NULL, ...){
      common <- c(human     = "hsapiens",
                  mouse     = "mmusculus",
                  rat       = "rnorvegicus",
                  zebrafish = "drerio",
                  elephant  = "lafricana",
                  fly       = "dmelanogaster",
                  fruitfly  = "dmelanogaster",
                  pig       = "sscrofa",
                  sheep     = "oaries",
                  rabbit    = "ocuniculus",
                  roundworm = "celegans",
                  vervet    = "csabaeus",
                  worm      = "celegans",
                  yeast     = "scerevisiae")
   if( tolower(dataset) %in% names(common))  dataset <- common[[tolower(dataset)]]
   if( !grepl("gene_ensembl$", dataset) )  dataset <- paste0(dataset, "_gene_ensembl")
   release <- version
   if (is.null(version)) {
        x <- biomaRt::listEnsembl()
        release <- x$version[x$biomart == "ensembl"]
        release <- gsub("Ensembl Genes ", "", release)
        message("Using Ensembl release ", release)
    }

  ## LIST
   if(!is.null(list)){
      if(tolower(list) == "datasets"){
         ensembl <- biomaRt::useEnsembl(biomart="ensembl", host=host, version=version)
         bm <- biomaRt::listDatasets(ensembl)
         ## remove AsIs class
         for(i in 1:3) class(bm[,i]) <- "character"
         message("Downloaded ", nrow(bm), " datasets")
      }else{
         ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset, host=host, version=version)
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
      ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset, host=host, version=version)
      # default search
      if(missing(attributes)){
          bm <- biomaRt::getBM(attributes=c('ensembl_gene_id','external_gene_name', 'gene_biotype',
                  'chromosome_name', 'start_position', 'end_position','strand', 'description',
                  'transcript_count'), mart = ensembl, ...)
           # replace long names like ensembl_gene_id
           names(bm)[1:6] <- c("id", "gene_name", "biotype", "chromosome", "start", "end")
           n <-   length(unique(bm$id))
           # drop source from description  [Source:MGI Symbol;Acc:MGI:102478]
           bm$description <- gsub(" \\[.*\\]$", "" , bm$description)
           # white space in version 92
          bm$description <- trimws(bm$description)
          bm <- dplyr::arrange(bm, id)
          if(!patch){
              n1 <- nrow(bm)
              bm <- dplyr::filter(bm, substr(chromosome,1,4) != "CHR_")
              if(n1 != nrow(bm)) message("Removed ", n1 - nrow(bm), " features on patch CHR_*")
         }
      }else{
         bm <- biomaRt::getBM(attributes= attributes, mart = ensembl, ...)
         if(!patch){
            if("chromosome_name" %in% colnames(bm)){
              n1 <- nrow(bm)
              bm <- dplyr::filter(bm, substr(chromosome_name,1,4) != "CHR_")
              if(n1 != nrow(bm)) message("Removed ", n1 - nrow(bm), " features on patch CHR_*")
          }
        }
      }
      message("Downloaded ", nrow(bm), " features")
   }
   # will also drop the grouped_df class
   bm <- tibble::as_tibble(bm)
   attr(bm, "downloaded") <- Sys.Date()
   attr(bm, "version") <- release
   bm
}
