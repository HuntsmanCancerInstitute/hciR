# pasilla dataset

library("pasilla")
counts   <- read_tsv(system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla") )
samples <- read_csv(system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla" ))
# Need sample data column matching count column names
samples$file <- gsub("fb$", "", samples$file )
## remove 2240 features with 0 reads and 721 with only 1 read
counts  <- filter_counts( counts, sum=TRUE)
dds <- deseq_from_tibble(counts, samples,  design = ~ condition)
rld <- rlog(dds)
fly <- read_biomart("fly")
res <- results_all(dds, fly)
pasilla <- list( dds = dds, rlog = rld, results = res)
