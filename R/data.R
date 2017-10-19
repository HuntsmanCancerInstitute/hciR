#' Canonical pathways in the MSigDB curated gene sets.
#'
#' Nine canonical pathways from the MSigDB C2 curated gene sets.  See \code{\link{read_msigdb}} to load file.
#'
#' @format A list with 9 pathways with 1329 gene sets
#'
#' @source \url{http://software.broadinstitute.org/gsea/msigdb/collections.jsp}
#' @examples
#' data(msig)
#' sapply(msig, length)
#' msig$KEGG[1:3]
"msig"

#' Human and Mouse Homologs
#'
#' Human and Mouse Homology with phenotype annotations from MGI
#'
#' @format A tibble with 18,501 rows and 6 variables:
#' \describe{
#'   \item{Human}{ Human Marker Symbol }
#'   \item{EntrezGene}{ Human Entrez Gene ID }
#'   \item{HomoloGene }{ HomoloGene ID }
#'   \item{Mouse}{ Mouse Marker Symbol }
#'   \item{MGI}{ MGI Marker Accession ID }
#'   \item{PhenotypeId}{ High-level Mammalian Phenotype ID  }
#' }
#' @source \url{http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt}.
#' @notes See the \code{data-raw} directory for code to download.
#' @examples
#' data(mgi)
#' mgi
#' table( mgi$Human == toupper(mgi$Mouse) )  # 86% are the same
#' table(duplicated(mgi$Mouse))
#' attr(mgi, "downloaded")
"mgi"

#' Human gene annotations from Ensembl release 90
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

#' Mouse gene annotations from Ensembl release 90
#'
#' @format A tibble with 50,143 rows and 11 columns
#' @source \code{read_biomart("mouse")}
#' @note Human homologs added from \code{\link{mgi}}
#' @examples
#' data(mmu)
#' mmu
#' group_by(mmu, biotype) %>%
#'  summarize(n=n(), human_homologs = sum(!is.na(human_homolog))) %>%
#'   arrange(desc(n))
"mmu"

#' Fruitfly gene annotations from Ensembl release 90
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
"kegg_hsa"

#' Mouse pathway gene sets from KEGG
#'
#' @format A tibble with 26,255 rows and 288 pathways
#' @note Signaling, metabolism and disease pathway indices are saved as attributes
#' @source \code{kegg.gsets} in \code{gage} package
#' @examples
#' data(kegg_mmu)
#' kegg_mmu
#' names(attributes(kegg_mmu))
#' filter(kegg_mmu, id %in% attr(kegg_mmu, "dise.idx")) %>%
#'   group_by(pathway) %>% summarize(n = n()) %>% arrange(desc(n))
"kegg_mmu"

#' DESeq objects and results for Pasilla knock-downs
#'
#' @format A list with DESeqDataSet (dds), DESeqTransform (rlog), and tibble (results)
#' @source \code{pasilla} package and Brooks et al. 2010. Conservation of an RNA regulatory map
#' between Drosophila and mammals.
#' @examples
#' data(pasilla)
#' pasilla
"pasilla"
