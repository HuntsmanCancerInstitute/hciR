#' Expand RColorBrewer palettes
#'
#' Expands palettes into 255 colors  for \code{plot_genes} and \code{plot_volcano}
#'
#' @param palette RColorBrewer palette name or "RdGn" for Red-Green color scale.
#' @param ramp use \code{colorRampPalette}, default returns 255 colors
#'
#' @return A vector of 255 colors, diverging palettes are reversed
#'
#' @author Chris Stubben
#'
#' @examples
#' x <- palette255("RdGn")
#' plot(c(1,256),c(0,3), type = "n", axes = FALSE,  xlab="", ylab="")
#' i <- 1:255
#' rect(i, 1, i+1, 2, col=x, border=NA)
#' @export

palette255 <- function(palette, ramp=TRUE){
   if( palette %in% c("RdGn", "RdGr")){
         clrs <- c(rev(RColorBrewer::brewer.pal(7,"Greens")), "white", RColorBrewer::brewer.pal(7,"Reds"))
    # OR reverse divergent color palette
    }else if(palette %in% c("BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral")){
       clrs <- rev( RColorBrewer::brewer.pal(11, palette))
    }else{
       clrs <-  RColorBrewer::brewer.pal(9, palette)
    }
    if(ramp) clrs <- grDevices::colorRampPalette(clrs)(255)
    clrs
}
