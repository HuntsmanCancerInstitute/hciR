### ADD human homologs to mmu

mmu <- read_biomart("mouse")

mmu_homologs <- read_biomart("mouse",  attributes = c("ensembl_gene_id", "external_gene_name",
  "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name",
 "hsapiens_homolog_perc_id"))

colnames(mmu_homologs) <- c("id", "gene_name", "homolog_id", "homolog", "perc_id")
mmu_homologs <- filter(mmu_homologs, homolog !="")

## GET homolog with max percent identity and combine into a comma-delimited list

x2 <- group_by(mmu_homologs, id) %>%
      filter(perc_id == max(perc_id)) %>%
       group_by(id) %>%
        summarize(n=n_distinct(homolog), human_homolog = paste(unique(homolog), collapse=","))

## 266 with more than 1 top homolog
table(x2$n>1)

## add to mmu

n <- match(mmu$id, x2$id)
mmu$human_homolog <- x2$human_homolog[n]
