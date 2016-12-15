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

#' KEGG pathway gene sets for mouse
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

#' KEGG pathway gene sets for human
#'
#' @format A named list with 291 pathways with Entrez ID vectors
#' @note Signaling, metabolism and disease pathway indices are saved as attributes
#' @source \code{kegg.gsets} in \code{gage} package
#' @examples
#' data(kegg)
#' kegg_hsa[1:3]
"kegg_hsa"
