#' Run enricher from clusterProfiler on all DESeq2 result tables
#'
#' @param res A list of DESeq results
#' @param gsets A list of gene sets (or two column dataframe with pathway and gene)
#' @param maxGSSize Maximum set size, default 2000
#' @param expressed_universe Use all expressed genes as universe, default TRUE
#' @param deseq.padj Adjusted p-value cutoff for significant genes, default 0.05
#' @param logFC absolute fold change cutoff
#' @param protein_coding compare protein coding genes only, default TRUE
#' @param quiet suppress messages except total significant
#' @param \dots Additional options like minGSSize and pvalueCutoff passed to \code{enricher}
#'
#' @return A list of tibbles
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(hciRdata)
#' enricher_all(res, msig_pathways$KEGG)
#' }
#' @export

enricher_all <- function(res, gsets, maxGSSize = 2000, expressed_universe=TRUE, deseq.padj = 0.05,
	logFC, protein_coding = TRUE, quiet=FALSE, ...){
   ## convert to 2 column data.frame
   if(class(gsets)[1] == "list"){
      term2gene <- dplyr::bind_rows(
		 lapply(gsets, tibble::enframe, value="gene"), .id="pathway") %>%
         dplyr::select(1,3)
  }else{
	  term2gene <- gsets
  }
   ##  if results are a tibble (since simplify=TRUE by default)
   if( class(res)[1] != "list"){
         n <- attr(res, "contrast")
         if(is.null(n)) stop("Missing contrast name")
         res <- list(res)
         names(res) <- n
   }
   n <- length(res)
   genes_res <- vector("list", n)
   vs <- names(res)
   ## padded for message
   vs1 <- sprintf(paste0("%-", max(nchar(vs))+2, "s"), paste0(vs, ":") )
   names(genes_res) <-  vs
   pc <- " "
   if(protein_coding) pc <- " protein coding "
   if(!missing(logFC)){
      if(!quiet) message( "Using padj < ", deseq.padj, " AND abs(log2FoldChange) > ", logFC, " to find significant", pc, "genes")
   }else{
      if(!quiet) message( "Using padj < ", deseq.padj, " to find significant", pc, "genes")
   }

   for(i in 1:n){
       r1 <- res[[i]]
       if(protein_coding){
           r1 <- dplyr::filter(r1, biotype == "protein_coding")
       }
       ## use first gene name in comma-separated lists of human homologs?
       if( "human_homolog" %in% colnames(r1)){
           # drop genes without homolog???
           r1 <- dplyr::filter(r1, !is.na(human_homolog))
           r1$human_homolog <- gsub(",.*", "", r1$human_homolog)
       }
       sig <- dplyr::filter( r1, padj <= deseq.padj)
       if( !missing(logFC)){
           sig <- dplyr::filter(sig, abs(log2FoldChange) >= logFC)
       }
       if(nrow(sig)==0){
		   message( i, ". ", vs1[[i]], " No significant genes")
		   genes_res[[i]]  <- tibble::tibble()
      }else{
          # get unique genes??
          if( "human_homolog" %in% colnames(r1)){
            sig_genes     <- unique(sig$human_homolog)
			universe     <- unique(r1$human_homolog)
          }else{
            sig_genes     <- unique(sig$gene_name)
			universe     <- unique(r1$gene_name)
          }
		  # use NULL to swith universe to all genes in term2gene
		  if(!expressed_universe) universe <- NULL
	      em1 <- clusterProfiler::enricher(sig_genes, TERM2GENE=term2gene, universe = universe, maxGSSize=maxGSSize, ...)
	      # ID and Description are the same
	      z <- tibble::as_tibble(em1)[, -2]
	      names(z)[1] <- "pathway"
          genes_res[[i]]  <- z
		  message( i, ". ", vs1[[i]], " ", nrow(z), " enriched sets")
	  }
   }
   genes_res <- genes_res[sapply(genes_res, nrow) > 0]
   if(length(genes_res) == 1) genes_res <- genes_res[[1]]
   genes_res
}
