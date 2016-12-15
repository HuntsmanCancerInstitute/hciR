
library(gage)

ks.mmu <- kegg.gsets(species = "mmu", id.type = "kegg")
## list of  1 list and 4 vectors
names(ks.mmu)
#[1] "kg.sets"    "sigmet.idx" "sig.idx"    "met.idx"    "dise.idx"

# add index vectors as attributes and save the list of entrez ids
kegg_mmu <- ks.mmu$kg.sets
for(i in 2:5) attr( kegg_mmu, names(ks.mmu)[i]) <- ks.mmu[[i]]

#---

ks.hsa <- kegg.gsets(species = "hsa", id.type = "kegg")
kegg_hsa <- ks.hsa$kg.sets
for(i in 2:5) attr( kegg_hsa, names(ks.hsa)[i]) <- ks.hsa[[i]]
