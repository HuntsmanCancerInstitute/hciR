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

#' Human gene annotations from Ensembl
#'
#' @format A tibble with 63,305 rows and 10 columns: id, gene_name, biotype, chromosome,
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

#' Mouse gene annotations from Ensembl
#'
#' @format A tibble with 50,143 rows and 10 columns
#' @source \code{read_biomart("mouse")}
#' @examples
#' data(mmu)
#' mmu
#' group_by(mmu, biotype) %>% summarize(n=n()) %>% arrange(desc(n))
"mmu"

#' Fruitfly gene annotations from Ensembl
#'
#' @format A tibble with 17,559 rows and 10 columns
#' @source \code{read_biomart("fly")}
#' @examples
#' data(mmu)
#' mmu
#' group_by(mmu, biotype) %>% summarize(n=n()) %>% arrange(desc(n))
"fly"

#' Human pathway gene sets from KEGG
#'
#' @format A tibble with 25,059 rows and 291 pathways
#' @note Signaling, metabolism and disease pathway indices are saved as attributes
#' @source \code{kegg.gsets} in \code{gage} package
#' @examples
#' data(kegg_hsa)
#' kegg_hsa
#' names(attributes(kegg_hsa))
#' filter(kegg_hsa, id %in% attr(kegg_hsa, "dise.idx")) %>%
#'   group_by(pathway) %>% summarize(n = n()) %>% arrange(desc(n))
#' \dontrun{
#' library(gage)
#' x <- kegg.gsets(species = "hsa")
#' n <- sapply(x$kg.sets, length)
#' y <- names(x$kg.sets)
#' kegg_hsa <- tibble( id = rep(1:length(y), n),  entry = rep(substr(y, 1,8), n),
#'   pathway = rep(substring(y, 10),n) , entrez = unlist(x$kg.sets) )
#' for(i in 2:5) attr(kegg_hsa, names(x1)[i]) <- x1[[i]]
#' }
"kegg_hsa"

#' DESeq objects and results for Pasilla knock-downs
#'
#' @format A list with DESeqDataSet (dds), DESeqTransform (rlog), and tibble (results)
#' @source \code{pasilla} package and Brooks et al. 2010. Conservation of an RNA regulatory map
#' between Drosophila and mammals.
#' @examples
#' data(pasilla)
#' pasilla
#' \dontrun{
#' library("pasilla")
#' counts   <- read_tsv(system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla") )
#' samples <- read_csv(system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla" ))
#' # Need sample data column matching count column names
#' samples$file <- gsub("fb$", "", samples$file )
#' ## remove 2240 features with 0 reads and 721 with only 1 read
#' counts  <- filter_counts( counts, sum=TRUE)
#' dds <- deseq_from_tibble(counts, samples,  design = ~ condition)
#' rld <- rlog(dds)
#' fly <- read_biomart("fly")
#' res <- results_all(dds, fly)
#' pasilla <- list( dds = dds, rlog = rld, results = res)
#' }
"pasilla"
