Download Ensembl annotations
================
Chris Stubben
October 16, 2019

The [hciRdata](https://github.com/HuntsmanCancerInstitute/hciRdata)
package includes feature annotations from the latest Ensembl release 98
for twelve reference genomes (human, mouse, elephant, fly, pig, rabbit,
rat, sheep, vervet, worm, yeast and zebrafish), plus earlier even
numbered releases for mouse and human (human90, human92, etc). Load the
package and select an annotation table.

``` r
library(hciRdata)
dplyr::select(human98, 1:4,8)
#  # A tibble: 60,623 x 5
#     id              gene_name biotype       chromosome description                                             
#     <chr>           <chr>     <chr>         <chr>      <chr>                                                   
#   1 ENSG00000000003 TSPAN6    protein_codi… X          tetraspanin 6                                           
#   2 ENSG00000000005 TNMD      protein_codi… X          tenomodulin                                             
#   3 ENSG00000000419 DPM1      protein_codi… 20         dolichyl-phosphate mannosyltransferase subunit 1, catal…
#   4 ENSG00000000457 SCYL3     protein_codi… 1          SCY1 like pseudokinase 3                                
#   5 ENSG00000000460 C1orf112  protein_codi… 1          chromosome 1 open reading frame 112                     
#   6 ENSG00000000938 FGR       protein_codi… 1          FGR proto-oncogene, Src family tyrosine kinase          
#   7 ENSG00000000971 CFH       protein_codi… 1          complement factor H                                     
#   8 ENSG00000001036 FUCA2     protein_codi… 6          alpha-L-fucosidase 2                                    
#   9 ENSG00000001084 GCLC      protein_codi… 6          glutamate-cysteine ligase catalytic subunit             
#  10 ENSG00000001167 NFYA      protein_codi… 6          nuclear transcription factor Y subunit alpha            
#  # … with 60,613 more rows
```

To download additional species, use `read_biomart` to view a list of all
available datasets.

``` r
x <- read_biomart(list="datasets")
Using Ensembl release 98
Downloaded 212 datasets
> x
# # A tibble: 212 x 3
#   dataset                      description                                  version
#   <chr>                        <chr>                                        <chr>
# 1 abrachyrhynchus_gene_ensembl Pink-footed goose genes (ASM259213v1)        ASM259213v1
# 2 acalliptera_gene_ensembl     Eastern happy genes (fAstCal1.2)             fAstCal1.2
# 3 acarolinensis_gene_ensembl   Anole lizard genes (AnoCar2.0)               AnoCar2.0
# 4 acitrinellus_gene_ensembl    Midas cichlid genes (Midas_v5)               Midas_v5
# 5 ahaastii_gene_ensembl        Great spotted kiwi genes (aptHaa1)           aptHaa1
# 6 amelanoleuca_gene_ensembl    Panda genes (ailMel1)                        ailMel1
# 7 amexicanus_gene_ensembl      Mexican tetra genes (Astyanax_mexicanus-2.0) Astyanax_mexicanus-2.0
# # … with 205 more rows
```

Use the first letter of the genus and species name in the first column
to download features (common names are supported for the reference
genomes with existing release 98 tables).

``` r
panda98 <- read_biomart("amelanoleuca")
# Using Ensembl release 98
# Downloaded 23262 features
```

To download earlier versions, check the Ensembl archives for available
releases. Archived versions are available for 54, 67, and all versions
after 77.

``` r
biomaRt::listEnsemblArchives()
#              name     date                                url version current_release
# 1  Ensembl GRCh37 Feb 2014          http://grch37.ensembl.org  GRCh37
# 2      Ensembl 98 Sep 2019 http://sep2019.archive.ensembl.org      98               *
# 3      Ensembl 97 Jul 2019 http://jul2019.archive.ensembl.org      97
# ...
# 22     Ensembl 78 Dec 2014 http://dec2014.archive.ensembl.org      78
# 23     Ensembl 77 Oct 2014 http://oct2014.archive.ensembl.org      77
# 24     Ensembl 75 Feb 2014 http://feb2014.archive.ensembl.org      75
# 25     Ensembl 67 May 2012 http://may2012.archive.ensembl.org      67
# 26     Ensembl 54 May 2009 http://may2009.archive.ensembl.org      54
```

Add the version number to download earlier releases.

``` r
zb87 <- read_biomart("zebrafish", version=87)
# Downloaded 32266 features
```

In some cases, the default attributes may be missing. For example, the
`external_gene_name` is missing from early fly releases. Use the list
option to view the available attributes.

``` r
fly67 <- read_biomart("fly", version=67)
# Error in biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",  :
#   Invalid attribute(s): external_gene_name
# Please use the function 'listAttributes' to get valid attribute names
fly67_columns <- read_biomart("fly", version=67, list="attributes")
fly67_columns
#    name                   description                page
#  1 ensembl_gene_id        Ensembl Gene ID            feature_page
#  2 ensembl_transcript_id  Ensembl Transcript ID      feature_page
#  3 ensembl_peptide_id     Ensembl Protein ID         feature_page
#  4 description            Description                feature_page
# ...
# 13 external_gene_id       Associated Gene Name       feature_page
# ...
```

Replace the default `external_gene_name` with `external_gene_id` by
passing a vector of attributes and then follow the code in
`read_biomart` to update column names and
descriptions.

``` r
fly67 <- read_biomart("fly", version=67, attributes = c('ensembl_gene_id',
    'external_gene_id', 'gene_biotype', 'chromosome_name', 'start_position',
     'end_position','strand', 'description', 'transcript_count'))
names(fly67)[1:6] <- c("id", "gene_name", "biotype", "chromosome", "start", "end")
fly67$description <- gsub(" \\[Source.*", "", fly67$description)
```

<br>
