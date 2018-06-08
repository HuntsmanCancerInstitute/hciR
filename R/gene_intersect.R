#' Compare genes to sets
#'
#' Intersect genes and a list of gene sets
#'
#' @param genes a vector of genes
#' @param sets a list of gene sets
#' @param min Minimum number of overlaps, default 2
#'
#' @return a tibble with gene set names, size and number of overlapping genes
#'
#' @author Chris Stubben
#'
#' @examples
#' gene_intersect(letters[5:8], list(x=letters[3:6], y=letters[6:14], z=letters[10:15]))
#' @export

gene_intersect <- function(genes, sets, min=2){
message("Comparing ", length(genes), " genes to ", length(sets), " sets")
   n <- sapply(sets, function(x) length( dplyr::intersect(genes, x) ))
   y <- tibble(term = names(sets), size= sapply(sets, length), overlap=n)
   y <- mutate(y, percent = round( overlap/size*100,1))
   filter(y, overlap>=2) %>% arrange(desc(percent))
}
