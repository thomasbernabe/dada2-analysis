library(dada2);packageVersion("dada2")
library(Biostrings);packageVersion("Biostrings")
library(ShortRead);packageVersion("ShortRead")
library(ggplot2);packageVersion("ggplot2")
library(reshape2);packageVersion("reshape2")
library(RColorBrewer);packageVersion("RColorBrewer")

path <- "/Users/thomasbernabe/Documents/callahan-lab/DATA-ZymoATCCMix-32plex-Revio"
path.out <- "Figures/"
path.rds <- "RDS/"

fn <- c(
  "m84036_230702_205216_s2.MAS16S_Fwd_01--MAS16S_Rev_13.hifi_reads.fastq.gz",
  "m84036_230702_205216_s2.MAS16S_Fwd_02--MAS16S_Rev_14.hifi_reads.fastq.gz",
  "m84036_230702_205216_s2.MAS16S_Fwd_03--MAS16S_Rev_15.hifi_reads.fastq.gz",
  "m84036_230702_205216_s2.MAS16S_Fwd_04--MAS16S_Rev_16.hifi_reads.fastq.gz",
  "m84036_230702_205216_s2.MAS16S_Fwd_05--MAS16S_Rev_17.hifi_reads.fastq.gz",
  "m84036_230702_205216_s2.MAS16S_Fwd_06--MAS16S_Rev_18.hifi_reads.fastq.gz",
  "m84036_230702_205216_s2.MAS16S_Fwd_07--MAS16S_Rev_19.hifi_reads.fastq.gz",
  "m84036_230702_205216_s2.MAS16S_Fwd_08--MAS16S_Rev_20.hifi_reads.fastq.gz"
)

fn <- file.path(path, fn)

sample.names <- paste0("ATCC_Sample_", sprintf("%02d", 1:8))