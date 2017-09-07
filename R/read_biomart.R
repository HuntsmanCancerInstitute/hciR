#' Read biomart annotations
#'
#' Read Ensembl annotations using \code{biomaRt} library
#'
#' @param dataset The first letter of genus and full species name like btaurus, ecaballus, sscrofa from \code{listEnsembl}.
#' A few common names are accepted for human, mouse, rat, fruitfly, yeast and zebrafish.
#' @param version Ensembl version number for \code{useEnsembl}. Default is current version.
#' @param attributes vector of column names to pass to \code{getBM}, default ensembl_gene_id,
#'  external_gene_name, gene_biotype, chromosome_name, start_position, end_position, strand,
#'  description, transcript_count, and entrezgene
#' @param fragments drop fragments like GL456211.1, CHR_MG132_PATCH (any chromosome
#' name with 3 or more characters), default FALSE
#' @param \dots additional options passed to \code{getBM}
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
#'    mm <- read_biomart("mouse")
#' }
#' @export

read_biomart <- function( dataset="human" , version = NULL, attributes, fragments = FALSE, ...){

   common <- c(human = "hsapiens",
               mouse = "mmusculus",
               rat = "rnorvegicus",
               zebrafish = "drerio",
               fly = "dmelanogaster",
               fruitfly = "dmelanogaster",
               pig = "sscrofa",
               scerevisiae = "yeast")

   if( tolower(dataset) %in% names(common))  dataset <- common[[tolower(dataset)]]
   if( !grepl("gene_ensembl$", dataset) )  dataset <- paste0(dataset, "_gene_ensembl")
   ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset, version = version)

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
      # drop fragments GL456211.1, CHR_MG132_PATCH
      if(fragments){
           bm <- subset(bm , nchar(chromosome) < 3)
           message("Dropped ", n-nrow(bm),  " rows with chromosome fragments, ", nrow(bm), " total rows")
       }
   }else{
       bm <- biomaRt::getBM(attributes= attributes, mart = ensembl, ...)
       message("Downloaded ", nrow(bm), " features")
   }
   # drop source from description  [Source:MGI Symbol;Acc:MGI:102478]
   bm$description <- gsub(" \\[.*\\]$", "" , bm$description)
   # drop the grouped_df class
   bm <- tibble::as_data_frame(bm)
   attr(bm, "downloaded") <- Sys.Date()
   bm
}
