#' Run all results from a DESeq analysis using apeglm shrinkage
#'
#' @param object a DESeqDataSet
#' @param biomart annotations from \code{read_biomart} with column 1 matching row names in results
#' @param alpha the significance cutoff for the adjusted p-value cutoff (FDR)
#' @param add_columns a vector of biomart columns to add to result table, default
#'        gene_name, biotype, chromosome, description and human_homolog if present
#' @param trt Compare groups within trt group, default is first term in the design formula
#' @param simplify return a tibble if only 1 contrast present
#' @param \dots additional options passed to \code{results}
#'
#' @return A list of tibbles for each contrast
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' library(hciRdata)
#' res <- results_apeglm(pasilla$dds, fly98)
#' res
#' }
#' @export

results_apeglm <- function( object, biomart, alpha = 0.05, add_columns, trt, simplify=TRUE,  ...){
   message("Using adjusted p-value < ", alpha)

   if(missing(trt)){
      n <- as.character( DESeq2::design(object))
      ## [1] "~"   "trt"
      # drop intercept if "0 + trt"?
      n[2] <- gsub("^0 \\+ ", "", n[2])
      ## multiple terms in design formula  ~ trt + mouse
      if(grepl(" + ", n[2], fixed=TRUE)){
         n[2] <- gsub(" \\+.*", "", n[2])
      }
      trt <- n[2]
   }
   ## sample groups
   grps <- levels(object[[trt]])
   n <- length(grps)
   contrast <- utils::combn(grps, 2)[2:1,, drop=FALSE]
  # contrast names
   vs <- apply(contrast, 2, paste, collapse = " vs ")
   if(length(vs) == 0) stop("No contrasts found")
   ## padded for message
   vs1 <- sprintf(paste0("%-", max(nchar(vs))+2, "s"), paste0(vs, ":") )

   ## resultNames for contrasts
   resNames <- paste0( trt, "_", gsub(" ", "_", vs))
   ## number contrasts for each reference level 1 1 1 2 2 3
   n2 <- rep(1:(n-1), (n-1):1)

   res <- vector("list", length(vs))
   names(res) <- vs
   if(missing(add_columns)){
        add_columns <- c("gene_name", "biotype", "chromosome", "description")
        if(!missing(biomart)){
           if("human_homolog" %in% names(biomart)) add_columns <- c(add_columns, "human_homolog")
        }
   }
   message("Adding apeglm's shrunken estimates to log2FoldChange, saving unshrunken in MLE_log2FC")
   for(i in 1:(n-1)){
      if(i > 1){
       #  message("Setting ", grps[i], " as reference level")
        ## create new levels...
        object[[trt]] <- stats::relevel(object[[trt]], ref = grps[i])
        object <- suppressMessages( DESeq2::nbinomWaldTest(object))
      }
      for(j in which(n2 %in% i)){
        res1 <- DESeq2::results(object, name=resNames[j] , alpha = alpha, ...)
          mle_log2FC <- res1$log2FoldChange
          res1 <-  DESeq2::lfcShrink(object, resNames[j], res=res1, type = "apeglm", quiet=TRUE)
          res1$MLE_log2FC <- mle_log2FC

          ft <- S4Vectors::metadata(res1)$filterThreshold
          x <- suppressMessages( summary_deseq(res1) )
          message(j, ". ", vs1[j], x[1,2], " up and ", x[2,2], " down regulated" )
          if(!missing(biomart)){
              # suppress messages like 70 rows in results are missing from biomart table and print once
              res1 <- suppressMessages( annotate_results( res1, biomart, add_columns) )
          }else{
             ## gene_name or id?
             res1 <- tibble::as_tibble(tibble::rownames_to_column(data.frame(res1), var="id"))
          }
          attr(res1, "contrast") <- vs[j]
          attr(res1, "alpha") <- alpha
          attr(res1, "filterThreshold") <- ft
          res[[j]] <- res1
      }
   }
   n1 <-is.na(res[[1]][[3]])
   if(any(n1)) message("Note: ", sum(n1), " rows in results are missing from biomart table")
   if(simplify & length(res) == 1)  res <- res[[1]]
   res
}
