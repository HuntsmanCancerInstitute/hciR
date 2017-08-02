#' Read IPA Expression Analysis output file
#'
#' Read an IPA Expression Analysis output text file and optionally
#' write to an  Excel file with each table in a separate tab.
#'
#' @param file IPA output file (created using Export All on the main Summary tab)
#' @param excel write to Excel file
#' @param mylists include optional My Pathways and My Lists
#'
#' @return a list of tibbles or Excel file
#'
#' @author Chris Stubben
#'
#' @examples
#'  \dontrun{
#'  x <- read_ipa("RPR2_vs_RPM.txt")
#'  names(x)
#'  read_ipa("RPR2_vs_RPM.txt",  excel=TRUE)
#'  }
#' @export

read_ipa <- function(file,  excel=FALSE, mylists=FALSE){
   x <- readLines(file)
   ## Find start of tables
   n <- grep("for My Projects", x)
   # save results
   n1 <- length(n)
   z <- vector("list", n1)
   ## get names for Excel tabs
   TABS <-  gsub( " for My .*", "", x[n])
   names(z) <- TABS
   ## add length of x
   n <- c(n, length(x)+1)
   # Loop through tables...
   for(i in 1:n1){
      start <- n[i] + 1
      end <- n[i + 1] - 2
      if(!mylists) if( TABS[i]  %in% c("My Lists", "My Pathways") ) next
      if(end - start == 0) next
      message("Loading ", TABS[i])
      y <- x[  start:end  ]
      if(i == 1){
           ## fix parsing problems with 1 or 3 columns to avoid warnings
           ## firt tried suppressWarnings, but not working within function, ok in Rstudio
           y[12] <- gsub("\tand ", " and ", y[12])
           y <- y[!y %in% c("", "Analysis Settings", "Filter Summary\t")]
           n2 <- grep("\t", y, invert=TRUE)
           y[n2] <- paste0(y[n2], "\t\t")  # 1 \t gets trimmed
           y <- c("Analysis Settings\tValue", y)
      }
      z[[i]] <- readr::read_tsv( paste( gsub("\t$", "", y), collapse ="\n"))
    }
    ## sort pathways by p-value
    z[["Canonical Pathways"]] <- arrange(z[["Canonical Pathways"]], desc(`-log(p-value)`))
    ## zscores, either #NUM! or lots of decimals in Excel
    z[["Canonical Pathways"]]$zScore  <-  sprintf("%.3f", z[["Canonical Pathways"]]$zScore)
    # drop Analysis name and empty ID
    z[["Upstream Regulators"]] <- z[["Upstream Regulators"]][,-(1:2)]
    z[["Networks"]] <-  z[["Networks"]][,-2]
    # drop missing tables , usually My lists
    z <- z[!sapply(z, is.null)]
    if(excel){
      z <- lapply(z, data.frame, check.names=FALSE)
      outfile <- gsub(".txt", ".xlsx", file)
      message("Saved to ", outfile)
      openxlsx::write.xlsx(z, file = outfile )
    }else{
      z
    }
}
