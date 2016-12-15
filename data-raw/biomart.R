# Ensembl annotations from biomart

library(biomaRt)

 ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

## ensembl gene id is primary key

mm <- getBM(attributes=c('ensembl_gene_id','gene_biotype', 'description', 'chromosome_name',
 'start_position','end_position','strand', 'external_gene_name', 'transcript_count'),
  mart = ensembl)
names(mm)[c(1,2,4:6, 8)] <- c("id", "biotype", "chromosome", "start", "end",  "gene_name")


## or with entrez gene id as primary key for GSEA and pathway mapping
## most ensembl IDs are 1:1 mapping, but 1414 ids are 1: many


mm_entrez <- getBM(attributes=c('entrezgene', 'ensembl_gene_id','gene_biotype', 'description', 'chromosome_name',
 'start_position','end_position','strand', 'external_gene_name', 'transcript_count'),
  mart = ensembl)
names(mm_entrez)[c(1,2,4:6, 8)] <- c("id", "biotype", "chromosome", "start", "end",  "gene_name")


##drop fragments GL456211.1, CHR_MG132_PATCH


table(mm$chromosome[nchar(chromosome) < 3])

   1   10   11   12   13   14   15   16   17   18   19    2    3    4    5    6    7    8    9   MT    X
3428 1639 2998 1537 1522 1802 1270 1191 1723  897 1079 3776 2893 2867 3252 3133 4721 2312 2180   37 2613
   Y
1570



 # drop source from description ?  ADD as new columns
  [Source:MGI Symbol;Acc:MGI:102478]

  table(gsub(".*\\[Source:([^ ;]+).*", "\\1" , mm$description))

             EntrezGene       HGNC        MGI     RefSeq
         869          9          2      47559          1


mm <- subset(mm , nchar(chromosome) < 3)
 mm$description <- gsub(" \\[.*\\]$", "" , mm$description)
 mm<- tbl_df(mm)


--
 mm_entrez <- subset(mm_entrez , nchar(chromosome) < 3)
  mm_entrez$description <- gsub(" \\[.*\\]$", "" , mm_entrez$description)
  mm_entrez<- tbl_df(mm_entrez)


save(mm, mm_entrez, file="~/Documents/R/hci/datasets/mouse_genes_ensembl.rda")

## 206!
mm4 <- subset(mm, chromosome_name=="4" & start_position > 77.8*1e6  & end_position < 94.7*1e6)



#------------------------------------------------
#  HUMAN

Homo sapiens genes (GRCh38.p7)

p8 (patch release 8) on 2016/06/30

library(biomaRt)
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

listFilters(ensembl)

tbl_df(listFilters(ensembl))
# A tibble: 310 x 2
             name     description
            <chr>           <chr>
1 chromosome_name Chromosome name
2           start Gene Start (bp)
3             end   Gene End (bp)
4      band_start      Band Start
5        band_end        Band End
6    marker_start    Marker Start
7      marker_end      Marker End
8   encode_region   Encode region
9          strand          Strand
# ... with 301 more rows


tbl_df(listAttributes(ensembl) )

tbl_df(listAttributes(ensembl) )
# A tibble: 1,464 x 3
                    name           description         page
                   <chr>                 <chr>        <chr>
1        ensembl_gene_id       Ensembl Gene ID feature_page
2  ensembl_transcript_id Ensembl Transcript ID feature_page
3     ensembl_peptide_id    Ensembl Protein ID feature_page
4        ensembl_exon_id       Ensembl Exon ID feature_page
5            description           Description feature_page
6        chromosome_name       Chromosome Name feature_page
7         start_position       Gene Start (bp) feature_page
8           end_position         Gene End (bp) feature_page
9                 strand                Strand feature_page
10                  band                  Band feature_page
# ... with 1,454 more rows



## 20 have more than one symbol

hgnc <- getBM(attributes=c('ensembl_gene_id','gene_biotype', 'description', 'chromosome_name',
 'start_position','end_position','strand', 'hgnc_symbol', 'external_gene_name', 'transcript_count'),
  mart = ensembl)

    names(hgnc)[c(1,2,4:6, 8:9)] <- c("id", "biotype", "chromosome", "start", "end",  "symbol", "gene_name")
    hgnc$description <- gsub(" \\[.*\\]$", "" , hgnc$description)
    hgnc<- tbl_df(hgnc)

    subset(hgnc, symbol!=gene_name & symbol!="") %>% dplyr::select(id, biotype, symbol, gene_name)

---


  hg <- getBM(attributes=c('ensembl_gene_id','gene_biotype', 'description', 'chromosome_name',
   'start_position','end_position','strand', 'external_gene_name', 'transcript_count'),
    mart = ensembl)

  names(hg)[c(1,2,4:6, 8)] <- c("id", "biotype", "chromosome", "start", "end",  "gene_name")


# drop source from description ?
RNA, U6 small nuclear 280, pseudogene [Source:HGNC Symbol;Acc:HGNC:47243]

hg$description <- gsub(" \\[.*\\]$", "" , hg$description)

class(hg)

# Grange

hg<- tbl_df(hg)

hg
# A tibble: 63,325 × 10
                id                biotype                           description       chromosome     start
             <chr>                  <chr>                                 <chr>            <chr>     <int>
1  ENSG00000252303                  snRNA RNA, U6 small nuclear 280, pseudogene CHR_HG2128_PATCH  67546651
2  ENSG00000281771               misc_RNA                                 Y RNA CHR_HG2128_PATCH  67631019
3  ENSG00000281256                lincRNA                                       CHR_HG2191_PATCH  74823667
4  ENSG00000283272                   sRNA                   Clostridiales-1 RNA CHR_HG2022_PATCH  91356877
5  ENSG00000280864                lincRNA                                        CHR_HG126_PATCH  72505075
6  ENSG00000280792                lincRNA                                       CHR_HG2233_PATCH 239734956
7  ENSG00000282878                lincRNA                                        CHR_HG986_PATCH  41242373
8  ENSG00000283276 unprocessed_pseudogene                                       CHR_HG2022_PATCH  91374331
9  ENSG00000281822                  snRNA  RNA, U1 small nuclear 62, pseudogene  CHR_HG126_PATCH  72577040
10 ENSG00000281384                lincRNA                                       CHR_HG2233_PATCH 239762860


table2(hg$biotype)
# A tibble: 45 x 2
                     name     n
                    <chr> <int>
1          protein_coding 22285
2    processed_pseudogene 10815
3                 lincRNA  7845
4               antisense  5703
5  unprocessed_pseudogene  3297
6                misc_RNA  2387
7                   snRNA  2047
8                   miRNA  1580
9                     TEC  1073
10                 snoRNA  1006


save(hg, file="~/Documents/R/hci/datasets/human_genes_ensembl.rda")
