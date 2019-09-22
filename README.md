# LDProject
Linkage Disequilibrium Project - Bachelor's Thesis - ICS FORTH

# Background

The goal of this Bachelor's Thesis is to optimize an existing tool that detects genomic loci that have been introgressed between Homo Sapiens and Neanderthal. The existing tool determines if specific genomic loci has been introgressed using Linkage Disequilibrium and windows to specify each genomic loci. SNPs are being compared using Pattern Matching and PopCount.

In this thesis the concept is to check another approach at SNP storing level. More specifically, binary Tries have been introduced as they are perfect for retrieving data faster, by approaching SNPs as words in a dictionary and storing each 'letter' hierarchically, avoiding also duplicates. This practice could help lower the computation time needed for traversing SNPs, which is very important as millions of traverses and comparisons are needed during the run of the algorithm.

# How to Run

1. Download a chromosome from 1000genomes database.
   
   Example for chromosome 22:
   
   wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
   
   gunzip ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
   
   (PS. To run for other chromosome you need to change file name in run.sh file)
   
2. Execute ./run.sh

# Metrics

This tool runs by default for 10000 patterns. To change number of patterns, change PATTERNS_MAX in the header file.

10000 patterns occupy 3GB of RAM and the execution time is around a minute.

The more the patterns, the more RAM is used, so be gentle changing it. Use htop to monitor changes in RAM used.
