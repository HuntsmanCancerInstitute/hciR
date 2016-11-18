#' Read biomart annotations
#'
#' Read Ensembl annotations using \code{biomaRt} library
#'
#' @param dataset ensembl dataset name from \code{listEnsembl}
#'
#' @return A tibble with id, biotype, description, chromosome, start, end, strand
#'  gene name and transcript count
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#'    mm <- read_biomart("mouse")
#' }
#' @export

read_biomart <- function( dataset="human" , fragments = TRUE){
   if(dataset == "mouse") dataset="mmusculus_gene_ensembl"
   if(dataset == "human") dataset="hsapiens_gene_ensembl"

   ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset=dataset)
   ## extra columns will usually create duplicate ensembl_gene_id rows
   bm <- biomaRt::getBM(attributes=c('ensembl_gene_id','gene_biotype', 'description', 'chromosome_name',
            'start_position','end_position','strand', 'external_gene_name', 'transcript_count'),
             mart = ensembl)
   names(bm)[c(1,2,4:6, 8)] <- c("id", "biotype", "chromosome", "start", "end",  "gene_name")
   # drop fragments GL456211.1, CHR_MG132_PATCH
   if(!fragments)  bm <- subset(bm , nchar(chromosome) < 3)
   # drop source from description  [Source:MGI Symbol;Acc:MGI:102478]
   bm$description <- gsub(" \\[.*\\]$", "" , bm$description)
   bm <- dplyr::tbl_df(bm)
   attr(bm, "downloaded") <- Sys.Date()
   bm
}
