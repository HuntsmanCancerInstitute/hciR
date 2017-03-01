#' Extract and annotate all results from a DESeq analysis
#'
#' Extract all possible contrasts and annotate result tables from a DESeq object.
#' Currently supports simple designs with a single variable.
#'
#' @param object a DESeqDataSet
#' @param biomart annotations from \code{read_biomart} with column 1 matching row names in results
#' @param add a vector of biomart columns to add to result table, default gene_name, biotype and description
#' @param vs1 either compare all vs all (default) or a specific treatment vs all
#' @param alpha the significance cutoff for the adjusted p-value cutoff (FDR)
#' @param simplify return a tibble if only 1 contrast present
#' @param \dots additional options passed to \code{results}
#'
#' @return A list of tibbles for each contrast
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' data(pasilla)
#' res <- results_all(dds, fly)
#' res
#' # Set factor levels in the DESeq object to change contrast order
#' #  since results_all uses combn on levels
#' dds$trt <- factor(dds$trt, levels=c("heart", "lung", "control"))
#' apply( combn(c("heart", "lung", "control"), 2), 2, paste, collapse= " vs. ")
#' }
#' @export

results_all <- function( object, biomart, add, vs1= "all", alpha = 0.05, simplify=TRUE,  ...){
   message("Using adjusted p-value < ", alpha)
   n <- as.character( DESeq2::design(object))
   ## [1] "~"   "condition"
   if(grepl(" + ", n[2], fixed=TRUE)){
      message("The design has multiple variables and only the first variable will be used")
      n[2] <- gsub(" \\+.*", "", n[2])
   }
   trt <- n[2]
   n <- levels( object[[trt]] )
   # add option to re-level ?

   if(vs1 == "all"){
      contrast <- utils::combn(n, 2)
    }else{
       if(!vs1 %in% n) stop("No level in ", trt, " named ", vs1)
       contrast <- rbind( vs1, n[n!=vs1])
    }
      vs <- apply(contrast, 2, paste, collapse = " vs. ")
      res <- vector("list", length(vs))
      names(res) <- vs

   ## padded for message
   vs1 <- sprintf(paste0("%-", max(nchar(vs))+2, "s"), paste0(vs, ":") )

   for(i in seq_along( vs )){
       res1 <- DESeq2::results(object, contrast = c( trt, contrast[1,i], contrast[2,i] ), alpha = alpha, ...)
       x <- summary_deseq(res1)
       message(i, ". ", vs1[i], x[1,2], " up and ", x[2,2], " down regulated" )
       if(!missing(biomart)){
          if(missing(add)) add <- c("gene_name", "biotype", "description")
           # suppress messages  like 70 rows in results are missing from biomart table and print once
             res1 <- suppressMessages( annotate_results( res1, biomart, add) )
             attr(res1, "contrast") <- vs[i]
            attr(res1, "alpha") <- alpha
       }
       res[[i]] <- res1
   }
   n1 <-is.na(res[[1]][[3]])
   if(any(n1)) message("Note: ", sum(n1), " rows in results are missing from biomart table")
   if(simplify & length(res) == 1)  res <- res[[1]]
   res
}
