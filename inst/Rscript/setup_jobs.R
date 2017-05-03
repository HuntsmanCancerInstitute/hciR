#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

opts <-  list(
   make_option(c("-e", "--email"), default="NA",
         help="#e email or hci user name for pysano directive, required"),
   make_option(c("-c", "--cluster"), default="kingspeak",
         help="#c cluster name with at least 64 GB of RAM (not ember), default kingspeak"),
   make_option(c("-a", "--analysis"), default="NA",
         help="#a Copy results to this Analysis directory"),
   make_option(c("-r", "--run"), default="NA",
         help="Run ID with Fastq files to link"),
   make_option(c("-d", "--database"), default="human",
         help="Reference database, either human, mouse, fly or rat, default human"),
    make_option(c("-l", "--length"), default="50",
       help="Read length, default 50")
)

parser <- OptionParser(option_list=opts, description = "
Creates a cmd.txt file and sample directories with Fastq file links in order to run STAR,
featureCounts and quality metrics on the CHPC clusters using pysano.  The default is to align
to the human reference configured for 50 bp reads in /tomato/dev/data/Human/GRCh38/star50")
  opt <- parse_args(parser)

  if( opt$email == "NA" ){
     print_help(parser)
     quit(status=1)
  }


if(file.exists( "cmd.txt")){
   message("Note: cmd.txt file already exists")
}else{

   if( !grepl("@", opt$email )) opt$email <- paste0( opt$email, "@hci.utah.edu")
   release <- 87
if( grepl( "mouse|mus", opt$database, ignore.case = TRUE) ) {
    assembly <- "GRCm38"
    name <- "Mus_musculus"
    db <- "Mouse"
}else if(grepl( "human|homo", opt$database, ignore.case = TRUE)){
   assembly <- "GRCh38"
   name <- "Homo_sapiens"
   db <- "Human"
}else if(grepl( "fly|dros", opt$database, ignore.case = TRUE) ){
   assembly <- "BDGP6"
   name <- "Drosophila_melanogaster"
   db <- "Drosophila"
}else if(grepl( "rat", opt$database, ignore.case = TRUE) ){
   assembly <- "Rnor_6.0"
   name <- "Rattus_norvegicus"
   db <- "Rat"
}else{
   stop("Database should be human, mouse, rat, or fly")
}

## Configure header
cat(
 paste("#e", opt$email),
 paste("#c", opt$cluster),
 if( opt$analysis != "NA") paste("#a", opt$analysis),
 "",
 paste0("DB=/tomato/dev/data/", db, "/", assembly),
 paste0("INDEX=$DB/star", opt$length),
 paste0("GTF=$DB/", name, ".", assembly, ".", release, ".gtf"),
       "CHROM=$DB/chrom.sizes",
 paste0("REFFLAT=$DB/", name, ".", assembly, ".", release, ".refflat"),
 paste0("RIBOINT=$DB/", name, ".", assembly, ".rRNA.interval"),
   sep="\n",  file = "cmd.txt")

cat("
## Paths
APP=/tomato/dev/app
FASTQC=$APP/FastQC/0.11.5/fastqc
STAR=$APP/STAR/2.5.2b/STAR
FEATCOUNT=$APP/Subread/1.5.1/bin/featureCounts
SAMTOOLS=$APP/samtools/1.3/samtools
BIGWIG=$APP/UCSC/bedGraphToBigWig
RnaSeqMetrics=$APP/picard/1.87/CollectRnaSeqMetrics.jar

## GET sample prefix from Fastq file name
GZ=`echo *.txt.gz`
OUT=`echo ${GZ%%_*}`

# rename for multiqc
mv $GZ $OUT.fq.gz
gunzip $OUT.fq.gz

# FASTQ v 0.11.5
$FASTQC -f fastq $OUT.fq

#  STAR v 2.5.2b
$STAR --genomeDir $INDEX \\
--readFilesIn $OUT.fq \\
--runThreadN $NCPU \\
--twopassMode Basic \\
--outSAMtype BAM SortedByCoordinate \\
--outWigType bedGraph \\
--outWigStrand Unstranded \\
--clip3pAdapterSeq AGATCGGAAGAGCACACGTCTGAACTCCAGTCA

# rename for multiqc ID parsing
mv Aligned.sortedByCoord.out.bam $OUT.bam
mv Log.final.out $OUT.Log.final.out

# featureCounts v 1.5.1
$FEATCOUNT -T 16 -s 2 -a $GTF -o $OUT.counts $OUT.bam

# Samtools v 1.3.1
$SAMTOOLS index $OUT.bam
$SAMTOOLS idxstats $OUT.bam | sort -V > $OUT.idxstats

## RnaSeq metrics - v 1.87
java -Xmx20G -jar $RnaSeqMetrics REF_FLAT=$REFFLAT \\
STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$RIBOINT \\
I=$OUT.bam  O=$OUT.rna_metrics

# bedGraphToBigWig v 4
$BIGWIG Signal.Unique.str1.out.bg  $CHROM $OUT.unique.bw
$BIGWIG Signal.UniqueMultiple.str1.out.bg $CHROM $OUT.multiple.bw

rm $OUT.fq
rm Signal.Unique.str1.out.bg
rm Signal.UniqueMultiple.str1.out.bg
## STAR twopassMode
rm -rf _STARgenome
rm -rf _STARpass1
",  file="cmd.txt",  append=TRUE)

message("Created cmd.txt file for ", db, " star", opt$length)
}


if( opt$run == "NA" ){
   message("Add Run ID to link Fastq files")
}else{
   repo_found <- FALSE
   for(year in 2017:2012){
      repo <- paste("/Repository/MicroarrayData", year, opt$run, sep="/")
      if( file.exists (repo)){
        repo_found <- TRUE
        break
     }
   }
  if(!repo_found) stop("No match to ", opt$run, " in /Repository/MicroarrayData")
  fastq <- list.files(paste0(repo, "/Fastq"), pattern = "txt.gz$")
  if(length(fastq) == 0) stop("No Fastq files found in ", opt$run)
  fastq_full <- list.files(paste0(repo, "/Fastq"), pattern = "txt.gz$", full.names =TRUE)
  ids <- gsub("([^_]+)_.*", "\\1", fastq)

   for(i in seq_along(fastq)){
     if( !dir.exists(ids[i])) dir.create(ids[i])
     file.symlink( fastq_full[i], ids[i])
     file.link( "cmd.txt", paste0(ids[i], "/cmd.txt"))
   }
   message("Linked ",  length(fastq), " Fastq files in ", repo)
}
