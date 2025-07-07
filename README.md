# dada2-analysis

## Overview
This repository contains scripts and results for comparing three different error estimation models in DADA2 when analyzing PacBio Kinnex full-length 16S rRNA sequencing data. The study evaluates how differrent error correction approaches affect the accuracy of microbiome analysis. 

## Dataset Description
- The PacBio Kinnex dataset is [DATA-ZymoATCCMix-32plex-Revio](https://downloads.pacbcloud.com/public/dataset/Kinnex-16S/DATA-ZymoATCCMix-32plex-Revio/)
- Download the following ATCC even strain files
<img width="296" alt="32-plex barcodes for ATCC" src="https://github.com/user-attachments/assets/86186606-8cab-42df-9116-295940ede1e2" />

- Place FASTQ files in data/ directory 

## File Structure 
Raw data from PacBio Kinnex dataset: 
- m84036_230702_205216_s2.MAS16S_Fwd_01--MAS16S_Rev_13.hifi_reads.fastq.gz
- m84036_230702_205216_s2.MAS16S_Fwd_02--MAS16S_Rev_14.hifi_reads.fastq.gz
- ... (8 total files)
  
Each file contains full-length 16S amplicons from teh same mock community

## Error Model Comparison 
1. Default Loess Model
   - DADA2's standard smooth fitting across quality scores
2. Binned Quality Score Model
   - Piecewise linear fitting for discrete quality bins
   - Uses makeBinnedQualErrfun() with data detected bins  
3. Custom Modified Loess Model
   - Based on Gulliem Salazar's solution (GitHub issue #938)
   - Manual tuning for improved performance 
  
## The order to run files 
1. Setup.R
   - Installs required packages and establishes the working directory
3. QualityFiltering.R
   - Filters reads based on quality scores
   - Preserves full-length reads 
4. CheckBinnedQ.R
   - Analyze quality score distribution
   - Creates appropriate bins for binned error function
5. ErrorModels.R
   - Learns with the three different error models
   - Generates error rate plots for comparison 
6. ASVAnalysis.R
   - Runs denoising with each error model
   - Compares ASV counts and characteristics 
