# KEGG pathway gene sets  for human
library(gage)

x <- kegg.gsets(species = "hsa")
n <- sapply(x$kg.sets, length)
y <- names(x$kg.sets)
kegg_hsa <- tibble( id = rep(1:length(y), n),  entry = rep(substr(y, 1,8), n),
  pathway = rep(substring(y, 10),n) , entrez = unlist(x$kg.sets) )
for(i in 2:5) attr(kegg_hsa, names(x1)[i]) <- x1[[i]]
