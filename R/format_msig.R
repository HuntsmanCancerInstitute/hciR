#' Format MsigDB names
#'
#' Formats MsigDB names with all upper case and underlines for spaces
#'
#' @param name a vector of names
#' @author Chris Stubben
#' @examples
#'  format_msig("THIS_IS_AN_MSIGDB_NAME")
#' @export

format_msig <- function( name ){
   name <- stringr::str_to_title(gsub("_", " ", name))
   CAPs <- c(Dna= "DNA", Rna = "RNA", Trna="tRNA",  Mrna="mRNA", Tca="TCA", ` Als` = " ALS", Abc = "ABC",
    ` And `=" and ", ` In ` = " in ", ` Of ` = " of ", ` The `= " the ", ` By ` = " by ",
      ` To `= " to ", ` For ` = " for ", ` From `= " from ", ` Or `= " or ", ` An `=" an ",
      ` Is ` = " is ", ` Ii `= " II ",  Iii = "III", ` Ri ` = " RI "  )
   stringr::str_replace_all(name, CAPs)
}
