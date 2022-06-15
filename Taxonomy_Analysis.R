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

#install.packages("cluster")
library(cluster)
#install.packages("factoextra")
library(factoextra)
#install.packages("forcats")
library(forcats)
#install.packages("ggfortify")
library(ggfortify)
#install.packages("phyloseq")
library(phyloseq)
#install.packages("plotly")
library(plotly)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("vegan")
library(vegan)

# Bioconductor

#BiocManager::install(c("MicrobiotaProcess"))
#library(MicrobiotaProcess)


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

# Store sample names in order.
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
# Kingdom  Phylum   Class   Order  Family   Genus  Species 
# 0        3        10      24     44       81     154 


# Remove variables that are not required downstream.
rm(i)


# Export dfTidy as a .tsv file for cleaner raw data. Code commented out as it does not need to be redone.
# write_tsv(dfTidy, file = "Tidy_Merged_Abundance_Table.tsv", col_names = T)



#### 3. Relative Abundance Visualization ----

## Create stacked bar chart by using ggplot and the Family taxa. First create a manipulated dataframe containing only the Family information. Due to the way our abundance table was created, we have to select the rows where the data is describing the abundance for the whole Family; the entries with Genus and Species specific s is encompassed in these Family abundance rows.
dfFamily <- dfTidy[,5:18]

# Remove rows with NAs in Family
dfFamily <- dfFamily[complete.cases(dfFamily$Family),]

# Keep rows with NAs in both Genus and Species - we only want these rows as the abundance information encompasses the whole Family.
dfFamily <- dfFamily[!complete.cases(dfFamily$Genus),]
dfFamily <- dfFamily[,c(1, 5:14)]

# Convert abundance information from character to numerical values.
dfFamily[, 2:11] <- lapply(dfFamily[, 2:11], function(x) as.numeric(as.character(x)))
sapply(dfFamily, class)

# Sum the abundance for each sample.
lapply(dfFamily[, 2:11], sum)


# Make the columns into rows - result checked to make sure the transformation was accurate.
dfFamily <- pivot_longer(dfFamily, cols = 2:11, names_to = "Sample")
colnames(dfFamily)[3] <- "Abundance"
dfFamily <- dfFamily %>%
  filter(Abundance != 0)


# Ideas of cleaning by Riffomonas Project https://www.youtube.com/watch?v=w4X3o6MQjVA

# Show a summary fo the relative abundance of the Family taxa.
dfFamily %>%
  group_by(Family) %>%
  summarise(max = max(Abundance)) %>%
  arrange(desc(max))

# Create an "Other" category based on a threshold of abundance 5%
dfFamily_pool <- dfFamily %>%
  group_by(Family) %>%
  summarise(pool = max(Abundance) < 10, 
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
colors <- RColorBrewer::brewer.pal(9,"Set3")
colors <- c("#debbcc", "#d4cce6", "#bad6ef", "#a8e1de", "#ebc9a7", "#bedcb8", "#dedede")

ggplot(dfFamily_Others, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors, 
                    breaks = c("Bifidobacteriaceae",
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
dfPhylum <- dfTidy[,2:18]

# Remove rows with NAs in Family
dfPhylum <- dfPhylum[complete.cases(dfPhylum$Phylum),]

# Keep rows with NAs in Class - we only want these rows as the abundance information encompasses the whole Family.
dfPhylum <- dfPhylum[!complete.cases(dfPhylum$Class),]
dfPhylum <- dfPhylum[,c(1, 8:17)]

# Convert abundance information from character to numerical values.
dfPhylum[, 2:11] <- lapply(dfPhylum[, 2:11], function(x) as.numeric(as.character(x)))
sapply(dfPhylum, class)

# Sum the abundance for each sample.
lapply(dfPhylum[, 2:11], sum)

# Save this dataframe to be used later to create a phyloseq object.
dfOTU <- dfPhylum


# Make the columns into rows - result checked to make sure the transformation was accurate.
dfPhylum <- pivot_longer(dfPhylum, cols = 2:11, names_to = "Sample")
colnames(dfPhylum)[3] <- "Abundance"
dfPhylum <- dfPhylum %>%
  filter(Abundance != 0)


# Show a summary fo the relative abundance of the Phylum taxa.
dfPhylum %>%
  group_by(Phylum) %>%
  summarise(max = max(Abundance)) %>%
  arrange(desc(max))

# Create an "Other" category based on a threshold of abundance 5%
dfPhylum_pool <- dfPhylum %>%
  group_by(Phylum) %>%
  summarise(pool = max(Abundance) < 5, 
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
colors <- RColorBrewer::brewer.pal(5,"Set3")
colors <- c("#debbcc", "#d4cce6", "#bad6ef", "#a8e1de", "#dedede")

ggplot(dfPhylum_Others, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = colors,
                    breaks = c("Actinobacteria",
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

#TODO can manipulate ggplot to have the top two taxa on either side of the stacked bar --> anchor for the top two, one at the top, and one at the bottom


#### 4. PCA ----

# Will be using the covariance matrix of this dataset as the units for all variables are the same. Have to transform dfData as prcomp() uses the columns as variables.

## Standardize the data as some taxa are found is high and low abundance, which may skew the analysis. Make the taxa the row names before transforming.
dfScaled <- dfData[,c(1, 3:12)]
rownames(dfScaled) <- dfScaled[,1]
dfScaled <- dfScaled[,-1]
dfScaled <- as.data.frame(t(dfScaled))

# Convert abundance data to numeric.
dfScaled[] <- lapply(dfScaled, function(x) as.numeric(as.character(x)))
sapply(dfScaled, class)

# Standardize data.
dfScaled <- as.data.frame(scale(dfScaled))

# Add column for samples to differentiate in visualizations.
dfScaled$Sample <- factor(rownames(dfScaled), levels = levels(samples))

# Perform PCA.
pca1 <- prcomp(dfScaled[,-326])

# Scree plot to visualize variance.
dfPCA1_var <- data.frame(PC = paste0("PC", 1:length(row.names(dfScaled))),
                          var = (pca1$sdev)^2 / sum((pca1$sdev)^2))

ggplot(dfPCA1_var[1:10,], aes(x = reorder(PC, -var), y = var)) +
  geom_bar(stat = "identity", fill = "darkslategray3") +
  ggtitle("Scree plot: PCA based on variance for taxonomy") +
  xlab("Principal Component") +
  ylab("Proportion of Variance") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0))


# Get scores for each PC.
summary(pca1)

# 2D PCA visualization
colors_pca1 <- c("#ef9e98", "#eaae7a", "#f7de92", "#aef9a1", "#96e3de",
                 "#9bd7fd", "#73adff", "#c8abe3", "#f68bf9", "#dedede")

autoplot(pca1, data = dfScaled,
         colour = "Sample",
         main = "PCA for Differentiation of Samples",
         size = 5) +
  theme_minimal() +
  scale_color_manual(values = colors_pca1) +
  geom_text(label = rownames(dfScaled),
            nudge_x = 0.07,
            cex = 3)
  
# 3D PCA visualization - to determine if the points grouped together are separated based on PC3 or are similar to one another. Create a dataframe for the PC scores used for this visualization.
dfPCA1_scores <- as.data.frame(pca1$x)
dfPCA1_scores$Sample <- factor(rownames(dfPCA1_scores), levels = levels(samples))

plot_ly(dfPCA1_scores, x = ~PC1, y = ~PC2, z = ~PC3, 
        color = dfPCA1_scores$Sample, 
        colors = colors_pca1,
        type = "scatter3d",
        mode = "markers+text") %>%
  layout(title = "3D PCA of Samples")


#### 5. Hierarchical Clustering ----

# Create numerical matrix from data.
matScaled <- as.matrix(dfScaled[,1:325])
matScaled <- matScaled[levels(samples),,drop = F]

# Crete dissimilarity matrix.
distDissim <- dist(matScaled)
matDissim <- as.matrix(distDissim)

# Visualize the dissimilarity matrix for the samples.
fviz_dist(distDissim,
          order = F,
          gradient = list(low = "darkslategray3",
                          mid = "white",
                          high = "salmon2")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Heatmap of dissimilarity between the samples") +
  xlab("") +
  ylab("") 

# Complete linkage clustering.
hc <- hclust(distDissim) # complete is the default

colors <- c("#f68bf9", "#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")
  
fviz_dend(hc,k = 5,
          cex = 0.8,
          k_colors = colors,
          rect = T,
          rect_fill = T,
          rect_border = colors,
          main = "Hierarchical Clustering of Samples using Complete Linkage")
          #ggtheme = theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())) # Add if y axis lines are required



#### 6. Statistical Analysis ----

# Create a matrix of the relative abundance data.
matData <- dfData[,c(1, 3:12)]
rownames(matData) <- matData[,1]
matData <- matData[,-1]
matData <- as.data.frame(matData)

# Convert abundance data to numeric.
matData[] <- lapply(matData, function(x) as.numeric(as.character(x)))
sapply(matData, class)
matData <- as.matrix(matData)



## Try to create Phyloseq object. Require three dataframes: OTU, taxonomy and samples table.

# OTU table needs to have the OTU names as a column and all the samples as columns. I will make the full taxonomy name as the OTU name (already created when making dfPhylum). Only using the Phylum information as this is the clade being used for analysis and the relative abundance for all samples at this stage equate to 100%.

# Changing the columns to go in order.
dfOTU <- relocate(dfOTU, levels(samples))
dfOTU <- relocate(dfOTU, Phylum, .before = levels(samples))
colnames(dfOTU)[1] <- "otu"

# Make the row names, the Phylum names. First, save the current row names into a varibale.
index <- rownames(dfOTU)
rownames(dfOTU) <- dfOTU$otu
matOTU <- as.matrix(dfOTU[,-1])
otu <- otu_table(matOTU, taxa_are_rows = T)


# Making the Taxonomy table, keeping the otu names column as the first (what is required).
dfTaxa <- tibble(otu = dfOTU$otu, dfTidy[index, 1:2])
rownames(dfTaxa) <- dfTaxa$otu
matTaxa <- as.matrix(dfTaxa)
tax <- tax_table(matTaxa)

# Making the Samples table. The frist column has to be the sample names.
dfSamples <- data.frame(sample = colnames(dfOTU)[2:11])
rownames(dfSamples) <- dfSamples$sample

# Create phyloseq object.
phy <- phyloseq(otu, tax, dfSamples)

# Make another phyloseq object with the full dataset.
dfOTU2 <- dfData[,c(1,3:12)]
dfOTU2 <- relocate(dfOTU2, levels(samples))
dfOTU2 <- relocate(dfOTU2, clade_name, .before = levels(samples))
colnames(dfOTU2)[1] <- "otu"
dfOTU2$otu <- paste0("OTU", 1:nrow(dfOTU2))
rownames(dfOTU2) <- dfOTU2$otu
matOTU2 <- as.matrix(dfOTU2[,-1])
storage.mode(matOTU2) <- "numeric"
otu2 <- otu_table(matOTU2, taxa_are_rows = T)

dfTaxa2 <- tibble(otu = dfOTU2$otu, dfTidy[,1:7])
rownames(dfTaxa2) <- dfTaxa2$otu
matTaxa2 <- as.matrix(dfTaxa2)
tax2 <- tax_table(matTaxa2)

phy2 <- phyloseq(otu2, tax2, dfSamples)


# Alpha diversity measure using Shannon index.
# "In our study, microbiota diversity was quantified using Shannon index. This diversity index is a quantitative indicator of the number of different bacteria that are present in a stool sample, taking into account the uniformity in the distribution of these bacteria in these species. Diversity index value increases both when the number of species increases and when evenness increases. The Shannon index is a well-known diversity index used in microecological studies. The higher the Shannon index value, the higher the community diversity ." Yin et al., 2019
plot_richness(phy, measures = "Shannon") +
  scale_x_discrete(limits = levels(samples)) +
  ggtitle("Alpha Diversity based on Phylum alone")

plot_richness(phy2, measures = "Shannon") +
  scale_x_discrete(limits = levels(samples)) +
  ggtitle("Alpha Diversity based on full data")

# There's only one point for each sample as we have not distinguished the samples by any variable.

















#### . REFERENCES ----

# Yin, L., Wan, Y. D., Pan, X. T., Zhou, C. Y., Lin, N., Ma, C. T., Yao, J., Su, Z., Wan, C., Yu, Y. W., & Zhu, R. X. (2019). Association Between Gut Bacterial Diversity and Mortality in Septic Shock Patients: A Cohort Study. Medical science monitor : international medical journal of experimental and clinical research, 25, 7376–7382. https://doi.org/10.12659/MSM.916808


