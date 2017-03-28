#' Extract and annotate all results from a DESeq analysis
#'
#' Extract all possible contrasts and annotate result tables from a DESeq object.
#' Currently supports simple designs with a single variable.
#'
#' @param object a DESeqDataSet
#' @param biomart annotations from \code{read_biomart} with column 1 matching row names in results
#' @param vs1 either compare all vs all (default) or a specific treatment vs all
#' @param alpha the significance cutoff for the adjusted p-value cutoff (FDR)
#' @param add_columns a vector of biomart columns to add to result table, default
#'        gene_name, biotype, chromosome, start and description
#' @param other_columns a vector of additional columns in biomart table
#' @param simplify return a tibble if only 1 contrast present
#' @param \dots additional options passed to \code{results}
#'
#' @note If you combine factors of interest into a single group following section 3.3 in the DESeq2 vignette,
#'  you can set vs1 = "combined1", "combined2", or "combined" to limit the comparisons.
#'  For example, if you combine 3 cell types and 2 time points, then the default
#' returns 18 contrasts,  "combined1" returns 3 contrasts comparing time within cell types, combined2 returns 6 contrasts
#' comparing cell types at the same time point or "combined" to return both (9 total).
#' Factors should be separated by spaces in the combined treatment group for parsing.
#'
#' @return A list of tibbles for each contrast
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' data(pasilla)
#' data(fly)
#' res <- results_all(pasilla$dds, fly)
#' res
#' # Set factor levels in the DESeq object to change contrast order
#' #  since results_all uses combn on levels
#' dds$trt <- factor(dds$trt, levels=c("heart", "lung", "control"))
#' apply( combn(c("heart", "lung", "control"), 2), 2, paste, collapse= " vs. ")
#' }
#' @export

results_all <- function( object, biomart,  vs1= "all", alpha = 0.05, add_columns, other_columns , simplify=TRUE,  ...){
   message("Using adjusted p-value < ", alpha)
   n <- as.character( DESeq2::design(object))
   ## [1] "~"   "condition"
   if(grepl(" + ", n[2], fixed=TRUE)){
      # message("The design has multiple variables and only the first variable will be used")
      n[2] <- gsub(" \\+.*", "", n[2])
   }
   trt <- n[2]
   n <- levels( object[[trt]] )
   # add option to re-level ?
   contrast <- utils::combn(n, 2)

   if( vs1 == "combined1"){
      ## if two columns are combined into a single trt group, compare first group
      n1 <- apply(contrast, 2, function(x) length(unique( gsub(" [^ ]+", "", x)))==1)
      contrast <- contrast[, n1]
   }else if( vs1 == "combined2"){
      ## or second group
      n2 <- apply(contrast, 2, function(x) length(unique( gsub("[^ ]+ ", "", x)))==1)
      contrast <- contrast[, n2]
   }else if( vs1 == "combined"){
      ## or both groups
      n1 <- apply(contrast, 2, function(x) length(unique( gsub(" [^ ]+", "", x)))==1)
      n2 <- apply(contrast, 2, function(x) length(unique( gsub("[^ ]+ ", "", x)))==1)
      contrast <- contrast[, n1 | n2]
    }else if( vs1 %in% n){
       contrast <- rbind( vs1, n[n!=vs1])
    }
      vs <- apply(contrast, 2, paste, collapse = " vs. ")
if(length(vs)==0) stop("No contrasts found")
      res <- vector("list", length(vs))
      names(res) <- vs

   ## padded for message
   vs1 <- sprintf(paste0("%-", max(nchar(vs))+2, "s"), paste0(vs, ":") )

   if(missing(add_columns)) add_columns <- c("gene_name", "biotype", "chromosome", "start", "description")
   # add one extra to defaults...
   if(!missing(other_columns)) add_columns <- c(add_columns, other_columns)

   for(i in seq_along( vs )){
       res1 <- DESeq2::results(object, contrast = c( trt, contrast[1,i], contrast[2,i] ), alpha = alpha, ...)
       x <- suppressMessages( summary_deseq(res1) )
       message(i, ". ", vs1[i], x[1,2], " up and ", x[2,2], " down regulated" )
       if(!missing(biomart)){
           # suppress messages  like 70 rows in results are missing from biomart table and print once
           res1 <- suppressMessages( annotate_results( res1, biomart, add_columns) )
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
