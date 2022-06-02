#### Microbiome Taxonomy Analysis

## Author: Shalvi Chirmade
## BINF*6999 Research Project Summer 2022

## Van Raay Lab at the University of Guelph

## Primary Advisor: Dr. Terry Van Raay
## Secondary Advisor: Dr. Lewis Lukens

## Description:
# Microbiome analysis of ten patient fecal samples for the identification of diet-derived metabolites and their effects on mental health.

# Please see the README file for the project proposal.

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#### 1. Data Import and Cleaning ----

# Data has been pre-processed through command-line on Compute Canada using:
      # FastQC 0.11.9
      # Trimmomatic 0.39
      # fastq-join 1.3.1

# Taxonomic analysis has been carried out on command-line using Bowtie2 and MetaPhlAn 3

# The taxonomic profile files for each sequence has been combined into one file called merged_abundance_table.txt

