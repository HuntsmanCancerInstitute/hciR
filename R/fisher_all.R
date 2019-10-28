#' Run Fisher's test on all DESeq2 result tables
#'
#' Compare significant genes to gene sets using a one-sided Fisher's test.
#'
#' @param res  A list of DESeq results
#' @param gsets Gene sets with human gene names.
#' @param deseq.padj Adjusted p-value cutoff for significant genes, default 0.05
#' @param logFC absolute fold change cutoff
#' @param min_set Minimum number of overlaps, default 2
#' @param protein_coding compare protein coding genes only, default TRUE
#' @param quiet suppress messages except total significant
#'
#' @return a tibble with gene sets intersections and pvalue from \code{fisher.test}.
#' Sets with no intersecting genes are dropped.
#'
#' @note Assumes gene sets are human genes and will intersect using human_homolog column if present.
#' Creates a 2 x 2 contingency table with significant and non-significant genes in each gene set and
#' the p-value is returned from a one-sided \code{fisher.test}
#'
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(hciRdata)
#' fisher_all(res, msig_pathways$KEGG)
#' }
#' @export

fisher_all <- function(res, gsets, deseq.padj = 0.05, logFC, min_set =2,
   protein_coding = TRUE, quiet=FALSE){
   if( any(sapply(gsets, function(x) any(is.na(x)))) ) stop( "Please remove NAs from gene sets")
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
   if(!missing(logFC)){
      if(!quiet) message( "Using padj < ", deseq.padj, " AND abs(log2FoldChange) > ", logFC, " to find significant genes")
   }else{
      if(!quiet) message( "Using padj < ", deseq.padj, " to find significant genes")
   }
   if(!quiet) message("Note: Sets with < ", min_set, " overlapping genes are removed")

   for(i in 1:n){

       r1 <- res[[i]]
       if(protein_coding){
           if(i ==1 & !quiet) message("Dropping ", sum(r1$biotype!= "protein_coding", na.rm=TRUE), " non-coding genes")
           r1 <- filter(r1, biotype == "protein_coding")
       }
       ## use first gene name in comma-separated lists of human homologs?
       if( "human_homolog" %in% colnames(r1)){
           # drop genes without homolog???
           r1 <- filter(r1, !is.na(human_homolog))
           r1$human_homolog <- gsub(",.*", "", r1$human_homolog)
       }
       sig <- filter( r1, padj <= deseq.padj)
       not_sig <- filter( r1, padj > deseq.padj)
       if( !missing(logFC)){
           not_sig <- dplyr::bind_rows(not_sig, filter(sig, abs(log2FoldChange) < logFC) )
           sig <- filter(sig, abs(log2FoldChange) >= logFC)
       }
       if(nrow(sig)==0) stop("No significant genes")

       up_reg <- filter( sig, log2FoldChange > 0)

      # get unique genes??
       if( "human_homolog" %in% colnames(r1)){
            sig_genes     <- unique(sig$human_homolog)
            up_reg_genes  <- unique(up_reg$human_homolog)
            not_sig_genes <- unique(not_sig$human_homolog)
       }else{
            sig_genes     <- unique(sig$gene_name)
            up_reg_genes  <- unique(up_reg$gene_name)
            not_sig_genes <- unique(not_sig$gene_name)
       }

      ## CREATE contigency table.   This is slower than intersect
      # tbls <- lapply(gsets, function(x) rbind( table(factor(sig_genes %in% x, levels=c(FALSE, TRUE))),
      #                                     table(factor(not_sig_genes %in% x, levels=c(FALSE, TRUE)) )))
      # intersect will ignore duplicates and NAs
      # intersect(c(1,1,2, NA), 1:4)
      ## OR c(1,1,2, NA) %in% 1:4
   #   n1 <- sapply(gsets, function(x) length( dplyr::intersect(sig_genes, x) ))
       n1 <- sapply(gsets, function(x) sum( sig_genes %in% x) )

      dropN <- n1 < min_set
   #   if(sum(dropN)>0) message("    Dropping ", sum(dropN), " gsets with < ", min_set, " hits")
      n1 <- n1[!dropN]
      if(length(n1) > 0){
         gsets1 <- gsets[!dropN]
         n2 <- sapply(gsets1, function(x) sum( not_sig_genes %in% x) )
         x <- cbind( n1, n2, x1 = length(sig_genes) - n1,  x2 = length(not_sig_genes) - n2)
         tbls <- lapply( split(x, 1:nrow(x)), matrix, ncol=2)
         ## Create a matrix with total = number of expressed genes
            ## set others
         # sig   49  2787
         #  ns   74 11467
         pvalue <- sapply(tbls, function(x) stats::fisher.test(x, alt="greater")$p.value)
         # count up-regulated genes in set for barplots...
         up <- sapply(gsets1, function(x) sum( up_reg_genes %in% x) )
         ### set size or sig+ns?
         set_size <- sapply(gsets1, length)
         x <- tibble::tibble(term = names(gsets1), size= set_size, overlap = n1+n2, signif = n1, up = up, down = n1-up, pvalue = pvalue)
         x <- dplyr::mutate(x, percent = round( signif/overlap*100,1))
         ## arrange by p-value?
         x <-  dplyr::arrange(x, pvalue)
         genes_res[[i]]  <- x
         #   message( i, ". ", vs1[[i]], sum(x$pvalue < 0.05), " signficant sets (", nrow(x), " total, ", length(sig_genes), " genes)")
         message( i, ". ", vs1[[i]], sum(x$pvalue < 0.05), " signficant sets (", nrow(x), " total)")
      }else{
        message( i, ". ", vs1[[i]], "No sets with minimum number of overlaps")
        genes_res[[i]]  <- tibble::tibble()
      }
   }
   genes_res <- genes_res[sapply(genes_res, nrow) > 0]
   if(length(genes_res) == 1) genes_res <- genes_res[[1]]
   genes_res
}
