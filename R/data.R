#' Human and Mouse Homologs
#'
#' Human and Mouse Homology with phenotype annotations from MGI
#'
#' @format A tibble with 18,463 rows and 6 variables:
#' \describe{
#'   \item{Human}{ Human Marker Symbol }
#'   \item{EntrezGene}{ Human Entrez Gene ID }
#'   \item{HomoloGene }{ HomoloGene ID }
#'   \item{Mouse}{ Mouse Marker Symbol }
#'   \item{MGI}{ MGI Marker Accession ID }
#'   \item{PhenotypeId}{ High-level Mammalian Phenotype ID  }
#' }
#' @source \url{http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt}
#' @examples
#' data(mgi)
#' mgi
#' table( mgi$Human == toupper(mgi$Mouse) )  # 85.7% are the same
"mgi"

#'  Mouse gene annotations from Ensembl
#'
#' @format A tibble with 50143 rows and 10 columns: id, gene_name, biotype, chromosome,
#' start, end, strand, description, transcript_count, entrez_gene
#' @source \code{read_biomart("mouse")}
#' @examples
#' data(mmu)
#' mmu
#' attr(mmu, "downloaded")
#' group_by(mmu, biotype) %>% summarize(n=n()) %>% arrange(desc(n))
#' # most ensembl ids have 1 entrez id
#' table(nchar(gsub("[^,]", "", mmu[["entrez_gene"]]))+1)
"mmu"

#'  Human gene annotations from Ensembl
#'
#' @format A tibble with 63305 rows and 10 columns: id, gene_name, biotype, chromosome,
#' start, end, strand, description, transcript_count, entrez_gene
#' @source \code{read_biomart("human")}
#' @examples
#' data(hsa)
#' hsa
#' attr(hsa, "downloaded")
#' group_by(hsa, biotype) %>% summarize(n=n()) %>% arrange(desc(n))
#' # most ensembl ids have 1 entrez id
#' table(nchar(gsub("[^,]", "", hsa[["entrez_gene"]]))+1)
"hsa"

#' Mouse pathway gene sets from KEGG
#'
#' @format A named list with 287 pathways with Entrez ID vectors
#' @note Signaling, metabolism and disease pathway indices are saved as attributes
#' @source \code{kegg.gsets} in \code{gage} package
#' @examples
#' data(kegg_mmu)
#' length(kegg_mmu)
#' kegg_mmu[1:3]
#' names(attributes(kegg_mmu))
#' length(kegg_mmu[attr(kegg_mmu, "sigmet.idx")] )
"kegg_mmu"

#' Human pathway gene sets from KEGG
#'
#' @format A named list with 291 pathways with Entrez ID vectors
#' @note Signaling, metabolism and disease pathway indices are saved as attributes
#' @source \code{kegg.gsets} in \code{gage} package
#' @examples
#' data(kegg)
#' kegg_hsa[1:3]
"kegg_hsa"


#' DESeqDataSet object with Pasilla knock-downs
#'
#' @format A DESeqDataSet object
#' @source \code{pasilla} package and Brooks et al. 2010. Conservation of an RNA regulatory map
#' between Drosophila and mammals.
#' @examples
#' data(pasilla_dds)
#' counts(pasilla_dds)[1:4, ]
#' \dontrun{
#' library("pasilla")
#' count_tbl   <- read_tsv(system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla") )
#' samples <- read_csv(system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla" ))
#' # Need sample data column matching count column names
#' samples$file <- gsub("fb$", "", samples$file )
#' ## remove 2240 features with 0 reads and 721 with only 1 read
#' count_tbl  <- filter_counts( count_tbl, sum=TRUE)
#' pasilla_dds <- deseq_from_tibble(count_tbl, samples,  design = ~ condition)
#' pasilla_dds <- DESeq(pasilla_dds)
#' }
"pasilla_dds"
