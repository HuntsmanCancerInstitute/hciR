#  Human and Mouse Homology with phenotype annotations from MGI

url <- "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"
# add empty column 8 to avoid parsing failures
mgi <- read_tsv(url, col_names=c("Human", "EntrezGene", "HomoloGene", "x1",  "Mouse", "MGI", "PhenotypeId", "x2"))
## 408 duplicated mouse genes
table(duplicated(mgi$Mouse))
# not sure about column 4?   There are 41=no  18460=yes
table(mgi$x1)
mgi <- mgi[, -c(4,8)]
attr(mgi, "downloaded") <- Sys.Date()
#save(mgi, file="mgi.rda")
