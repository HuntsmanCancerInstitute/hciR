#' Run Fisher's test on all DESeq2 result tables
#'
#' Compare significant genes to gene sets using a one-sided Fisher's test.
#'
#' @param res  A list of DESeq results
#' @param gsets Gene sets
#' @param deseq.padj Adjusted p-value cutoff for significant genes, default 0.05
#' @param min Minimum number of overlaps, default 2
#'
#' @return a tibble with gene sets intersections and pvalue from \code{fisher.test}.
#' Sets with no intersecting genes are dropped.
#'
#' @note This functions removes non-coding genes in results and creates a 2 x 2
#' contingency table with significant and non-significant genes in each gene set.
#' A one-sided \code{fisher.test} returns p-values for all sets with at least
#' the minimum number of overlapping genes.
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' data(msig)
#' fisher_all(res, msig$KEGG)
#' }
#' @export

fisher_all <- function(res, gsets, deseq.padj = 0.05, min =2){
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
   names(genes_res) <-  vs
   message( "Using padj < ", deseq.padj, " to find significant genes")
   for(i in 1:n){
       r1 <- filter(res[[i]], biotype == "protein_coding")
       sig <- filter( r1, padj <= deseq.padj)
       if(nrow(sig)==0) stop("No significant genes")
       not_sig <-  filter( r1, padj > deseq.padj)
       up_reg <- filter( sig, log2FoldChange>0)

       if( "human_homolog" %in% colnames(sig)){
            ## some human homologs are comma-separated
            sig_genes     <- unique(unlist(strsplit(sig$human_homolog, ",")) )
            not_sig_genes <- unique(unlist(strsplit(not_sig$human_homolog, ",")) )
            up_reg_genes  <- unique(unlist(strsplit(up_reg$human_homolog, ",")) )
       }else{
            sig_genes     <- unique(sig$gene_name)
            not_sig_genes <- unique(not_sig$gene_name)
            up_reg_genes  <- unique(up_reg$gene_name)
       }
      message( i, ". ", vs[[i]])
      message("   Found ", length(sig_genes), " significant genes (", length(not_sig_genes), " not signficant)")
      ## CREATE contigency table.   This is slower than intersect
      # tbls <- lapply(gsets, function(x) rbind( table(factor(sig_genes %in% x, levels=c(FALSE, TRUE))),
      #                                     table(factor(not_sig_genes %in% x, levels=c(FALSE, TRUE)) )))
      n1 <- sapply(gsets, function(x) length( dplyr::intersect(sig_genes, x) ))
      dropN <- n1 < min
      message("    Dropping ", sum(dropN), " gsets with < ", min, " hits")
      n1 <- n1[!dropN]
      gsets1 <- gsets[!dropN]

      n2 <- sapply(gsets1, function(x) length( dplyr::intersect(not_sig_genes, x) ))
      x <- cbind( x1 = length(sig_genes) - n1,  x2 = length(not_sig_genes) - n2, n1, n2)
      tbls <- lapply( split(x, 1:nrow(x)), matrix, ncol=2)

      pvalue <- sapply(tbls, function(x) stats::fisher.test(x, alt="less")$p.value)
      # count up-regulated genes in set for barplots...
      up <- sapply(gsets1, function(x) length( dplyr::intersect(up_reg_genes, x) ))
      ### set size or sig+ns?
      set_size <- sapply(gsets1, length)
      x <- tibble::data_frame(term = names(gsets1), total= n1+n2, signif = n1, up = up, down = n1-up, pvalue = pvalue)
      x <- dplyr::mutate(x, percent = round( signif/total*100,1))
      ## drop gsets with no overlapping genes?
      x <- dplyr::filter(x, signif > 1) %>% dplyr::arrange(pvalue)
      genes_res[[i]]  <- x
      message( "   ", nrow(x), " gsets with ", min, " or more genes (", sum(x$pvalue < 0.05), " with Fisher p-value < 0.05)" )

   }
   if(length(genes_res) == 1) genes_res <- genes_res[[1]]
   genes_res
}
