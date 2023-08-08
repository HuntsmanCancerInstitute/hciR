#' GDC column data
#'
#' Removes empty and duplicated columns in the \code{colData} slot from \code{GDCprepare}
#' and returns a tibble.
#'
#' @param gdc a SummarizedExperiment object from \code{GDCprepare} in the \code{TCGAbiolinks} package.
#'
#' @return A tibble
#'
#' @author Chris Stubben
#'
#' @examples
#' \dontrun{
#' brca <- GDCprepare(query)
#' meta <- gdc_coldata(brca)
#' meta
#' x <- summary_tibble(meta)
#' x
#' filter(x, unique==2)
#' }
#' @export

gdc_coldata <- function(gdc){
   x <- as.data.frame(SummarizedExperiment::colData(gdc))
   # treatments is a list of dataframes
   x <- x[, colnames(x) != "treatments"]
   # Fix other list columns... disease_type, primary_site
   ## table( sapply(x, class))
   n <- which(sapply(x, class) %in% c("list", "AsIs"))
   #for(i in n) x[[i]] <- unlist(x[[i]])
   for(i in n) x[[i]] <- sapply( x[[i]], paste, collapse=",")
   ## some columns with all values =  Not reported
   x[x == "not reported"] <- NA
   x <- tibble::as_tibble(x)  %>% dplyr::mutate_if(is.factor, as.character)
   ## Drop duplicate column names, state.x, state.y,  (all live?), updated.datetime
   ## dplyr::select(meta,  starts_with("state"))
     x <- dplyr::select(x,  -dplyr::starts_with("state"), -dplyr::starts_with("updated"))
   ## some samples have NAs
   if(all(x$name == x$disease_type, na.rm=TRUE)) x <- x[, colnames(x) != "name"]
   if(all(x$tissue_or_organ_of_origin == x$primary_diagnosis, na.rm=TRUE)) x <- x[, colnames(x) != "tissue_or_organ_of_origin"]
   if(all(x$site_of_resection_or_biopsy == x$primary_diagnosis, na.rm=TRUE)) x <- x[, colnames(x) != "site_of_resection_or_biopsy"]
   if(all(x$released, na.rm=TRUE)) x <- x[, colnames(x) != "released"]  ## all TRUE
   ## not currently exported - need to FIX
   x  <- suppressMessages( drop_empty_columns(x) )
   x
}
