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

#### 1. Load Required Packages ----

# CRAN

#install.packages("tidyverse")
library(tidyverse)

# Bioconductor

#BiocManager::install(c("MicrobiotaProcess"))
library(MicrobiotaProcess)


# GitHub



#### 2. Data Import and Cleaning ----

# Data has been pre-processed through command-line on Compute Canada using:
      # FastQC 0.11.9
      # Trimmomatic 0.39
      # fastq-join 1.3.1

# Taxonomic analysis has been carried out on command-line using Bowtie2 and MetaPhlAn 3

# The taxonomic profile files for each sequence has been combined into one file called merged_abundance_table.txt

# Import the merged taxa file.
dfData <- read.delim("merged_abundance_table.txt", header = F)
View(dfData)

# Row 1 contains the name of the database used for the taxonomy analysis - will be deleted.
dfData <- dfData[-1,]

# Now Row 1 contains the column names, will make those the names of the columns of the dataframe and delete the row.
colnames(dfData) <- dfData[1,]
colnames(dfData)
dfData <- dfData[-1,]

# Clean the column names to display only the patient number
colnames(dfData)[3:12] <- str_extract(colnames(dfData)[3:12], "Patient_[:alnum:]+")
colnames(dfData)
colnames(dfData)[3] <- str_extract(colnames(dfData)[3], "Patient_[:alpha:]+")
colnames(dfData)

# Create a new dataframe with the clade information separated into different columns.
dfTidy <- dfData
dfTidy <- dfTidy %>% 
  separate(clade_name, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|", 
           extra = "merge")

# Clean all the rank names - i.e. remove "x__".
for (i in 1:7){
  dfTidy[,i] <- sapply(str_split(dfTidy[,i], "__"), "[", 2)
} 






# The number of NAs in each taxonomic rank column
sapply(dfTidy[,1:7], function(x) sum(is.na(x)))
# Kingdom  Phylum   Class   Order  Family   Genus  Species 
# 0        3        10      24     44       81     154 





