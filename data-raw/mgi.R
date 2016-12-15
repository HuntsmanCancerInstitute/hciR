

#  Human and Mouse Homology with phenotype annotations from MGI

url <- "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"

# add empty column 8 to avoid paring failures
mgi <- read_tsv(url, col_names=c("Human", "EntrezGene", "HomoloGene", "x1",  "Mouse", "MGI", "PhenotypeId", "x2"))

# not sure about column 4?   There are 18412 = yes and 51= no
mgi <- mgi[, -c(4,8)]

#save(mgi, file="mgi.rda")
