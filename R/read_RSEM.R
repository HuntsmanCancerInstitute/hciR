#' Read RSEM output files
#'
#' Reads RSEM counts or stats files and optionally reshape into wide format.
#'
#' @param path the path to RSEM output files, the default corresponds to the working directory.
#' @param pattern regular expression for count file name matching, deafult .genes.results
#' @param reshape reshape into wide format with samples in rows (a count matrix).
#' @param stats read stat files, default counts
#'
#' @note The cnt and model files in the stats directory vary depending on RSEM options
#' and the parser may fail.   Reshape uses only expected counts or alignment stats
#'
#' @return A tibble in long or wide format if reshape=TRUE
#'
#' @author Chris Stubben
#' @references \url{https://github.com/deweylab/RSEM/blob/master/cnt_file_description.txt} and
#' \url{https://github.com/deweylab/RSEM/blob/master/model_file_description.txt} for output format
#'
#' @examples
#' \dontrun{
#'  # count matrix
#' rsem_counts <- read_RSEM( ".")
#' rsem_counts
#' # need to round expected counts for some functions
#'  rlog_values <- rlog( round( as.matrix(rsem_counts)))
#' rsem_all <- read_RSEM( ".", , reshape=TRUE)
#' rsem_all
#'  # create TPM or FPKM matrix
#' tpm <- dplyr::select( x, sample, gene_id, tpm) %>% tidyr::spread(sample, tpm)
#'
#'  # reshape uses alignments stats only (rows 1-3 in *.cnt files)
#' read_RSEM( ".", stats=TRUE)
#' # read all cnt and model stats into long format
#' rsem <- read_RSEM( ".", stats=TRUE, reshape=FALSE)
#'  plot Fragment length, Read length, Read start position, or Reads per alignment
#' x <- filter(rsem, stat=="Fragment length", n < 350)
#' hchart(x, "line", x=n ,  y= value, group=sample)
#'
#' xCol <- c("Unique", "Multiple", "Unaligned", "Filtered")
#' x <- filter(rsem, stat %in% xCol ) %>%
#'  mutate(stat = factor(stat, levels=xCol))
#' hchart(x, "bar", x=sample, y= value, group=stat) %>%
#'   hc_plotOptions(bar = list(stacking = "normal"))
#' }
#' @export

read_RSEM <- function(path = ".", pattern = "genes.results$", reshape = TRUE, stats = FALSE){
   if(!stats){
      res1 <- read_sample_files(path, pattern)
      if(reshape)  res1 <- dplyr::select( res1, sample, gene_id, expected_count) %>% spread(sample, expected_count)
   }else{
      cntF <- list.files(path, "\\.cnt$", recursive=TRUE, full.names=TRUE)
      if(length(cntF)==0) stop("No *.cnt files found")
      samples <- sample_names(cntF)
      out1 <- vector("list", length(cntF))
      for(i in seq_along(cntF)){
         message("Reading count and model files in ", gsub("/[^/]+$", "", cntF[i]))
         x1 <- read.table(file = cntF[i], nrow=3 , fill=TRUE)
         x1 <- as.matrix(x1)  # to select rows without unlisting
         ## Reads per alignment starting in row 4
         x2 <- read.table(file = cntF[i], skip=3, sep="\t")
         cnt <- bind_rows(
           tibble(row=1, stat=c("Unaligned", "Aligned",  "Filtered", "Total" ), n=1:4,  value=x1[1,1:4] ),
           tibble(row=2, stat=c("Unique",    "Multiple", "Uncertain" ),         n=1:3,  value=x1[2,1:3] ),
           tibble(row=3, stat=c("Hits"),                                        n=1,    value=x1[3,1] ),
           tibble(row=1:nrow(x2)+3, stat="Reads per alignment",                 n= x2[,1], value=x2[,2] )
                        )
         # add sample and file ending
         cnt <- tibble::add_column(cnt, sample=samples[i], file = "count", .before=1)
         ## read model
         x3 <- readLines(gsub("cnt$", "model", cntF[i]), n=14)
         vec1 <- as.numeric( strsplit(x3[5], " ")[[1]])
         vec2 <- as.numeric( strsplit(x3[8], " ")[[1]])
         ## Check if has_optional_length_dist =0 and then adjust coordinates for RSPD (add 1 or 2?)
         if(x3[8] == "0"){
            message("  Missing Fragment length distribution")
            rld <- tibble(row=6, stat= "Read length", n=(vec1[1] + 1):vec1[2],  value=as.numeric(strsplit(x3[6], " ")[[1]]) )
            fld <- tibble()
            x3 <- c("add extra element for rspd parsing", x3)
         }else{
            fld <- tibble(row=6, stat= "Fragment length", n=(vec1[1] + 1):vec1[2], value=as.numeric(strsplit(x3[6], " ")[[1]]) )
            rld <- tibble(row=9, stat= "Read length",     n=(vec2[1] + 1):vec2[2], value=as.numeric(strsplit(x3[9], " ")[[1]]) )
         }
         ## Read start position is off by default.
         if(x3[11] == 0 & x3[12] == ""){
            message("  Missing Read Start position, add --estimate-rspd to rsem-calculate-expression")
            rspd <- tibble()
         }else{
            rspd <- tibble(row=13, stat= "Read start position",  n=1:x3[12], value=as.numeric(strsplit(x3[13], " ")[[1]]) )
         }
         # TO DO get Quality scores
         model <- bind_rows(fld, rld, rspd)
         model <- tibble::add_column(model, sample=samples[i], file = "model", .before=1)
         out1[[i]] <- bind_rows(cnt, model)
      }
      res1 <- bind_rows(out1)
      if(reshape){
         res1 <-  filter(res1, file=="count", row <= 3 ) %>%
              dplyr::select(sample, stat, value) %>% spread(stat, value)
      }
  }
  res1
}
