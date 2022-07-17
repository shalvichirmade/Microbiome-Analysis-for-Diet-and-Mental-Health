#### Script for Microbiome Taxonomy Analysis

## Author: Shalvi Chirmade
## BINF*6999 Research Project Summer 2022

## Van Raay Lab at the University of Guelph

## Primary Advisor: Dr. Terry Van Raay
## Secondary Advisor: Dr. Lewis Lukens

## Description: This is a reproducible script for conducting microbiome taxonomy analysis. It allows for the creation of a taxonomy stacked barplot, PCA analysis and hierarchical clustering. Use CTRL+F to find the word "TODO" and edit the code accordingly to fit your use.

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#### 1. Load Required Packages ----

# CRAN

#install.packages("cluster")
library(cluster)
#install.packages("factoextra")
library(factoextra)
#install.packages("forcats")
library(forcats)
#install.packages("ggfortify")
library(ggfortify)
#install.packages("matrixStats")
library(matrixStats)
#install.packages("phyloseq")
library(phyloseq)
#install.packages("plotly")
library(plotly)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("vegan")
library(vegan)



#### 2. Data Import and Cleaning ----

# Data has been pre-processed through command-line on Compute Canada using:
# FastQC 0.11.9
# Trimmomatic 0.39
# fastq-join 1.3.1

# Taxonomic analysis has been carried out on command-line using Bowtie2 and MetaPhlAn 3

# The taxonomic profile files for each sequence has been combined into one file called merged_abundance_table.txt

# Import the merged taxa file.
#TODO Change the file name accordingly.
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

#TODO This was to remove the numbers after ASD (Patient_ASD977), can use these lines if you have any patient names with letters. Make sure to change the index of the column. If there are no patient names with letters, skip this step.
colnames(dfData)[3] <- str_extract(colnames(dfData)[3], "Patient_[:alpha:]+")
colnames(dfData)

# Store sample names in order.
#TODO Edit the order of columns based on your data.
samples <- sort(colnames(dfData)[3:12])
samples <- reorder(samples, c(1,9,2,3,4,5,6,7,8,10))

# Make rownames clean.
rownames(dfData) <- 1:(nrow(dfData))


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


# Remove genus name from species column.
dfTidy[,7] <- sapply(str_split(dfTidy[,7], "_", n = 2), "[", 2)


# The number of NAs in each taxonomic rank column
sapply(dfTidy[,1:7], function(x) sum(is.na(x)))


# Remove variables that are not required downstream.
rm(i)


# Export dfTidy as a .tsv file for cleaner raw data. Only needs to be run once; comment the code once completed to prevent overwriting.
write_tsv(dfTidy, file = "Tidy_Merged_Abundance_Table.tsv", col_names = T)



#### 3. Relative Abundance Visualization ----

## Create stacked bar chart by using ggplot and the Family taxa. First create a manipulated dataframe containing only the Family information. Due to the way our abundance table was created, we have to select the rows where the data is describing the abundance for the whole Family; the entries with Genus and Species specific s is encompassed in these Family abundance rows.
dfFamily <- dfTidy[,5:ncol(dfTidy)]

# Remove rows with NAs in Family
dfFamily <- dfFamily[complete.cases(dfFamily$Family),]

# Keep rows with NAs in both Genus and Species - we only want these rows as the abundance information encompasses the whole Family.
dfFamily <- dfFamily[!complete.cases(dfFamily$Genus),]
dfFamily <- dfFamily[,c(1, 5:ncol(dfFamily))]

# Convert abundance information from character to numerical values.
dfFamily[, 2:ncol(dfFamily)] <- lapply(dfFamily[, 2:ncol(dfFamily)], function(x) as.numeric(as.character(x)))
sapply(dfFamily, class)

# Sum the abundance for each sample.
round(sapply(dfFamily[, 2:ncol(dfFamily)], sum), 3)


# Make the columns into rows - result checked to make sure the transformation was accurate.
dfFamily <- pivot_longer(dfFamily, cols = 2:ncol(dfFamily), names_to = "Sample")
colnames(dfFamily)[3] <- "Abundance"
dfFamily <- dfFamily %>%
  filter(Abundance != 0)


# Ideas of cleaning by Riffomonas Project https://www.youtube.com/watch?v=w4X3o6MQjVA

# Show a summary fo the relative abundance of the Family taxa.
dfFamily %>%
  group_by(Family) %>%
  summarise(max = max(Abundance)) %>%
  arrange(desc(max))

dfFamily %>%
  group_by(Family) %>%
  summarise(mean = mean(Abundance)) %>%
  arrange(desc(mean))

# Create an "Other" category based on a threshold of abundance 5%
dfFamily_pool <- dfFamily %>%
  group_by(Family) %>%
  summarise(pool = max(Abundance) < 10, #TODO can edit value if needed
            mean = mean(Abundance),
            .groups = "drop")

dfFamily_Others <- inner_join(dfFamily, dfFamily_pool, by = "Family") %>%
  mutate(Family = if_else(pool, "Other", Family)) %>%
  group_by(Sample, Family) %>%
  summarise(Abundance = sum(Abundance), 
            mean = min(mean),
            .groups = "drop") %>%
  mutate(Family = factor(Family),
         Family = fct_reorder(Family, mean, .desc = T))


# Create stacked barplot using this new dataframe.
#TODO can edit colors accordingly
colors <- c("#debbcc", "#d4cce6", "#bad6ef", "#a8e1de", "#ebc9a7", "#bedcb8", "#dedede")

ggplot(dfFamily_Others, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors, 
                    breaks = c("Bifidobacteriaceae", #TODO edit names accordingly
                               "Coriobacteriaceae",
                               "Enterococcaceae",
                               "Lachnospiraceae",
                               "Ruminococcaceae",
                               "Streptococcaceae",
                               "Other")) +
  scale_x_discrete(limits = levels(samples)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Relative Abundance of Family") +
  ylab("Relative Abundance")

# Remove variables no longer required.
rm(dfFamily, dfFamily_Others, dfFamily_pool)


## Carry out same analysis using Phylum level.
dfPhylum <- dfTidy[,2:ncol(dfTidy)]

# Remove rows with NAs in Family
dfPhylum <- dfPhylum[complete.cases(dfPhylum$Phylum),]

# Keep rows with NAs in Class - we only want these rows as the abundance information encompasses the whole Family.
dfPhylum <- dfPhylum[!complete.cases(dfPhylum$Class),]
dfPhylum <- dfPhylum[,c(1, 8:ncol(dfPhylum))]

# Convert abundance information from character to numerical values.
dfPhylum[, 2:ncol(dfPhylum)] <- lapply(dfPhylum[, 2:ncol(dfPhylum)], function(x) as.numeric(as.character(x)))
sapply(dfPhylum, class)

# Sum the abundance for each sample.
round(sapply(dfPhylum[, 2:ncol(dfPhylum)], sum), 3)


# Make the columns into rows - result checked to make sure the transformation was accurate.
dfPhylum <- pivot_longer(dfPhylum, cols = 2:ncol(dfPhylum), names_to = "Sample")
colnames(dfPhylum)[3] <- "Abundance"
dfPhylum <- dfPhylum %>%
  filter(Abundance != 0)


# Show a summary fo the relative abundance of the Phylum taxa.
dfPhylum %>%
  group_by(Phylum) %>%
  summarise(max = max(Abundance)) %>%
  arrange(desc(max))

dfPhylum %>%
  group_by(Phylum) %>%
  summarise(mean = mean(Abundance)) %>%
  arrange(desc(mean))

# Create an "Other" category based on a threshold of abundance 5%
dfPhylum_pool <- dfPhylum %>%
  group_by(Phylum) %>%
  summarise(pool = max(Abundance) < 5, #TODO can edit value is needed
            mean = mean(Abundance),
            .groups = "drop")

dfPhylum_Others <- inner_join(dfPhylum, dfPhylum_pool, by = "Phylum") %>%
  mutate(Phylum = if_else(pool, "Other", Phylum)) %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance),
            mean = min(mean),
            .groups = "drop") %>%
  mutate(Phylum = factor(Phylum),
         Phylum = fct_reorder(Phylum, mean, .desc = T))


# Create stacked barplot using this new dataframe.
#TODO can edit colros accordingly
colors <- c("#debbcc", "#d4cce6", "#bad6ef", "#a8e1de", "#dedede")

ggplot(dfPhylum_Others, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors,
                    breaks = c("Actinobacteria", #TODO edit names accordingly
                               "Bacteroidetes",
                               "Firmicutes",
                               "Verrucomicrobia",
                               "Other")) +
  scale_x_discrete(limits = levels(samples)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Relative Abundance of Phylum") +
  ylab("Relative Abundance")

# Remove variables no longer required.
rm(dfPhylum, dfPhylum_Others, dfPhylum_pool)
















