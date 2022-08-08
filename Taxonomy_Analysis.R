#### Bioinformatic Characterization of Microbiomes Related to Gut-Health and Mental Well-being

## Author: Shalvi Chirmade
## BINF*6999 Bioinformatics Masters Project - Summer 2022

## Van Raay Lab at the University of Guelph

## Primary Advisor: Dr. Terry Van Raay
## Secondary Advisor: Dr. Lewis Lukens

## Description:
# Microbiome analysis of ten patient fecal samples for the identification of diet-derived metabolites and their effects on mental health.

#https://github.com/shalvichirmade/Microbiome-Analysis-for-Diet-and-Mental-Health

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
#install.packages("ggeasy")
library(ggeasy)
#install.packages("ggfortify")
library(ggfortify)
#install.packages("ggVennDiagram")
library(ggVennDiagram)
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
round(sapply(dfFamily[, 2:11], sum), 3)


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

dfFamily %>%
  group_by(Family) %>%
  summarise(mean = mean(Abundance)) %>%
  arrange(desc(mean))

# Create an "Other" category based on a threshold of abundance 10%
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
round(sapply(dfPhylum[, 2:11], sum), 3)

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

dfPhylum %>%
  group_by(Phylum) %>%
  summarise(mean = mean(Abundance)) %>%
  arrange(desc(mean))

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

# Eigen vectors for this PCA.
summary(pca1)

# Remove variables no longer required.
rm(pca1, dfPCA1_var, dfPCA1_scores)


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
  
fviz_dend(hc, k = 5,
          cex = 0.8,
          k_colors = colors,
          rect = T,
          rect_fill = T,
          rect_border = colors,
          main = "Hierarchical Clustering of Samples using Complete Linkage")
          #ggtheme = theme_minimal() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank())) # Add if y axis lines are required

# As the value of k decreases, the groups on the right come together first, leaving samples 9,6 5 to be added last.

# plot(agnes(distDissim,
#            method = "average"),
#      main = "B. Complete Linkage using Agnes",
#      xlab = "Differentially Expressed Genes",
#      which.plot = 2,
#      cex = 0.5)

# Remove variables no longer required.
rm(hc)


#### 6. K-medoid Clustering ----

# Using the elbow method to find the optimal cluster (value of k).
# Takes a few minutes to run.
fviz_nbclust(t(matScaled), pam, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2) +
  labs(subtitle = "WSS Elbow method for Kmedoid (samples)")
# Optimal cluster at 3

# Create the k-medoid plot.
medoid3 <- pam(matScaled, 3)
medoid3
medoid3$clustering # Patients 6 and 9 are their own clusters
table(medoid3$clustering)

fviz_cluster(medoid3, data = matScaled, 
             main = "K-Medoid cluster plot", 
             repel = TRUE)

summary(medoid3)



#### 7. PCA Comparison to Publicly Available Data ----

# Dataset was obtained from MicrobiomeDB. It is the maturation of the human gut microbiome during the first 5 years of life in children living in Bangladesh. 
# https://microbiomedb.org/mbio/app/record/dataset/DS_01668ecdbf

# This dataset will be used to for comparison against the ten samples in this study. A PCA plot will be used to show the differences between the healthy children in the obtained dataset and the children of the ten samples.

# Input the data.
dfExtDetails <- read_tsv("Bangladesh_healthy_5yr.16s_DADA2.sample_details.tsv")
dfExtAbund <- read_tsv("Bangladesh_healthy_5yr.16s_DADA2.taxon_abundance.tsv")

# What is the oldest age of a child form this dataset?
max(dfExtDetails$`Age at sample collection (months)`) # 60.4 - 5rs of age

# Remove children below 2 of age as the lowest age of the ten samples is 2 years.
sapply(dfExtDetails, class)
dfExtDetails <- dfExtDetails %>%
  filter(`Age at sample collection (months)` > 24)

unique(dfExtDetails$`Host body product`) # all fecal samples

# Rename sample column header.
colnames(dfExtDetails)[1] <- "Sample"

# How many samples from each age group?
dfExtDetails %>%
  filter(`Age at sample collection (months)` > 24 & `Age at sample collection (months)` < 36.1) %>%
  count() # 438

dfExtDetails %>%
  filter(`Age at sample collection (months)` > 36 & `Age at sample collection (months)` < 48.1) %>%
  count() # 602

dfExtDetails %>%
  filter(`Age at sample collection (months)` > 48 & `Age at sample collection (months)` < 61) %>%
  count() # 543


# Randomly select 10 samples from each age group for further analysis.
set.seed(6999)

dfExtDetails_Subset <- dfExtDetails %>%
  filter(`Age at sample collection (months)` > 24 &
           `Age at sample collection (months)` < 36.1) %>%
  sample_n(10)

dfExtDetails_Subset <- rbind(dfExtDetails_Subset, dfExtDetails %>%
                               filter(`Age at sample collection (months)` > 36 &
                                        `Age at sample collection (months)` < 48.1) %>%
                               sample_n(10))

dfExtDetails_Subset <- rbind(dfExtDetails_Subset, dfExtDetails %>%
                               filter(`Age at sample collection (months)` > 48 &
                                        `Age at sample collection (months)` < 60.4) %>%
                               sample_n(10))


# Subset the Abundance table according to the samples selected.
dfExtAbund_Subset <- dfExtAbund %>% 
  select(...1 ,dfExtDetails_Subset$Sample)
colnames(dfExtAbund_Subset)[1] <- "Taxa"


# Separate the Taxa column into individual columns by clade.
dfExtAbund_Subset <- dfExtAbund_Subset %>%
  separate(Taxa, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";")

# Remove genus name from species column.
dfExtAbund_Subset[,7] <- sapply(dfExtAbund_Subset[,7], function(x) word(x, 2))


# Replace NA in count table with 0.
dfExtAbund_Subset[,8:37] <- sapply(dfExtAbund_Subset[,8:37], function(x) replace_na(x, 0))


# Make count data into relative abundance data. Create a function to change each value based on the sum of the vector (column of dataframe).
make_relative <- function(numvec){
  # Store the sum of the vector into a variable.
  sumvec <- sum(numvec)
  
  # Divide each element of the vector by the total sum and multiply by 100 - this creates relative abundance data.
  numvec <- sapply(numvec, function(x) (x/sumvec)*100)
  
  return(numvec)
}

dfExtAbund_Subset[,8:37] <- sapply(dfExtAbund_Subset[,8:37], make_relative)


# How many taxa match between dfTidy and dfExtAbund_Subset?
nrow(inner_join(dfTidy, dfExtAbund_Subset, by = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) #29

dfJoin <- inner_join(dfTidy, dfExtAbund_Subset, 
                     by = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# Remove the taxa columns as that information is not required for PCA.
dfJoin <- dfJoin[,9:48]

# Convert abundance data to numeric.
dfJoin[] <- lapply(dfJoin, function(x) as.numeric(as.character(x)))
sapply(dfJoin, class)

# Transpose the data so the samples are rows.
dfJoin <- t(dfJoin)

# Standardize data.
dfJoin <- as.data.frame(scale(dfJoin))

# Add column for samples to differentiate in visualizations.
dfJoin$Sample <- rep(c("Sample", "External"), c(10, 30))

# Perform PCA.
pca2 <- prcomp(dfJoin[,-30])

# Scree plot to visualize variance.
dfPCA2_var <- data.frame(PC = paste0("PC", 1:29),
                         var = (pca2$sdev)^2 / sum((pca2$sdev)^2))

ggplot(dfPCA2_var[1:10,], aes(x = reorder(PC, -var), y = var)) +
  geom_bar(stat = "identity", fill = "darkslategray3") +
  ggtitle("Scree plot: PCA based on variance for taxonomy") +
  xlab("Principal Component") +
  ylab("Proportion of Variance") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0))


# Get scores for each PC.
summary(pca2)

# 2D PCA visualization
colors_pca2 <- c("darkslategray3", "lightsalmon2")

autoplot(pca2, data = dfJoin,
         colour = "Sample",
         main = "PCA for Differentiation of Samples",
         size = 5) +
  theme_minimal() +
  scale_color_manual(values = colors_pca2) +
  geom_text(label = c(rownames(dfJoin)[1:10], rep("", 30)),
            nudge_x = 0.04,
            cex = 3) 

# Decided not to use this PCA as part of the analysis as the taxa used is a very small subset of the actual information available. Will be using the whole taxa information for analysis.

# Remove variables no longer required.
rm(dfJoin, dfPCA2_var, pca2, colors_pca2)


## Try again using the FULL taxonomic information.
dfTidy[,9:18] <- sapply(dfTidy[,9:18], function(x) as.numeric(as.character(x)))

#full_join from dplyr did not work --> it gave a lot of errors when converting NAs to 0's. Some were becoming 0's and some were becoming 0.0000's. The single )'s became NaN during scaling. Used baseR merge() instead to conduct a full-join.
dfFullJoin <- merge(dfTidy, dfExtAbund_Subset, 
                    by = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    all = T)
dfFullJoin <- dfFullJoin[,9:48]
dfFullJoin[is.na(dfFullJoin)] <- 0

# Delete rows with all 0's. Otherwise it produces errors during scaling.
dfFullJoin <- dfFullJoin[rowSums(dfFullJoin[]) > 0,]

# Transpose the data so the samples are rows.
dfFullJoin <- t(dfFullJoin)

# Standardize data.
dfFullJoin <- as.data.frame(scale(dfFullJoin))

# Add column for samples to differentiate in visualizations.
dfFullJoin$Sample <- rep(c("Sample", "External"), c(10, 30))

# Perform PCA.
pca3 <- prcomp(dfFullJoin[,-602])

# Scree plot to visualize variance.
dfPCA3_var <- data.frame(PC = paste0("PC", 1:40),
                         var = (pca3$sdev)^2 / sum((pca3$sdev)^2))

ggplot(dfPCA3_var[1:10,], aes(x = reorder(PC, -var), y = var)) +
  geom_bar(stat = "identity", fill = "darkslategray3") +
  ggtitle("Scree plot: PCA based on variance for taxonomy from both external dataset and the ten samples") +
  xlab("Principal Component") +
  ylab("Proportion of Variance") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0))


# Get scores for each PC.
summary(pca3)

# 2D PCA visualization
colors_pca3 <- c("darkslategray3", "lightsalmon2")

autoplot(pca3, data = dfFullJoin,
         colour = "Sample",
         main = "PCA for Differentiation of Samples",
         size = 5) +
  theme_minimal() +
  scale_color_manual(values = colors_pca3) +
  geom_text(label = c(rownames(dfFullJoin)[1:10], rep("", 30)),
            nudge_x = 0.04,
            cex = 3)+
  guides(color = guide_legend(override.aes = list(size = 5)),
         size = F)

# 3D PCA visualization - to determine if the points grouped together are separated based on PC3 or are similar to one another. Create a dataframe for the PC scores used for this visualization.
dfPCA3_scores <- as.data.frame(pca3$x)
dfPCA3_scores$Sample <- c(rownames(dfPCA3_scores)[1:10], rep("External", 30))
dfPCA3_scores$Sample <- factor(rownames(dfPCA3_scores), levels = levels(samples))
dfPCA3_scores$Sample <- addNA(dfPCA3_scores$Sample)
levels(dfPCA3_scores$Sample) <- c(levels(dfPCA3_scores$Sample), "External")
dfPCA3_scores$Sample[is.na(dfPCA3_scores$Sample)] <- "External"

plot_ly(dfPCA3_scores, x = ~PC1, y = ~PC2, z = ~PC3, 
        color = dfPCA3_scores$Sample, 
        colors = c(colors_pca1, "lightsalmon2"),
        type = "scatter3d",
        mode = "markers+text") %>%
  layout(title = "3D PCA of Samples")

# Remove variables no longer required.
rm(pca3, colors_pca1, colors_pca3, dfPCA3_var, dfPCA3_scores, 
   dfFullJoin, dfExtAbund, dfExtAbund_Subset, dfExtDetails, dfExtDetails_Subset)


#### 8. Evaluation of Pipeline ----

# The following are the terminal outputs regarding the pre-processing and MetaPhlAn steps of the sequence file being evaluated against the created pipeline.

## Trimmomatic terminal output:
  # Input Read Pairs: 4664804 Both Surviving: 4427724 (94.92%) Forward Only Surviving: 157515 (3.38%) Reverse Only Surviving: 50893 (1.09%) Dropped: 28672 (0.61%)
  # TrimmomaticPE: Completed successfully

## Fastq-join terminal output:
  # Total reads: 4427724
  # Total joined: 1703390
  # Average join len: 73.67
  # Stdev join len: 39.45
  # Version: 1.3.1

## MetaPhlAn terminal output:
  # WARNING: The metagenome profile contains clades that represent multiple species merged into a single representant.
  # An additional column listing the merged species is added to the MetaPhlAn output.



### Import and clean the taxonomy file from the evaluation, after MetaPhlAn.

# Import the merged taxa file.
dfEvaluate <- read.delim("N_S002_9_S57_merged_profile.txt", header = F)
View(dfEvaluate)

# Row 1, 2 and 3 contains the name of the database used for the taxonomy analysis - will be deleted.
dfEvaluate <- dfEvaluate[-c(1,2,3),]

# Now Row 1 contains the column names, will make those the names of the columns of the dataframe and delete the row.
colnames(dfEvaluate) <- dfEvaluate[1,]
colnames(dfEvaluate)
dfEvaluate <- dfEvaluate[-1,]

# Make rownames clean.
rownames(dfEvaluate) <- 1:(nrow(dfEvaluate))

# Clean column names.
colnames(dfEvaluate)[1] <- "clade_name"


# Create a new dataframe with the clade information separated into different columns.
dfEvaluate_Tidy <- dfEvaluate
dfEvaluate_Tidy <- dfEvaluate_Tidy %>% 
  separate(clade_name, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
           sep = "\\|", 
           extra = "merge")


# Clean all the rank names - i.e. remove "x__".
for (i in 1:7){
  dfEvaluate_Tidy[,i] <- sapply(str_split(dfEvaluate_Tidy[,i], "__"), "[", 2)
} 


# Remove genus name from species column.
dfEvaluate_Tidy[,7] <- sapply(str_split(dfEvaluate_Tidy[,7], "_", n = 2), "[", 2)


# The number of NAs in each taxonomic rank column
sapply(dfEvaluate_Tidy[,1:7], function(x) sum(is.na(x)))
# Kingdom  Phylum   Class   Order  Family   Genus  Species 
# 0        2        7       17     28       50     88 


# Remove variables that are not required downstream.
rm(i)


# Separate out the additional species column into the number of species present.
str_count(dfEvaluate_Tidy$additional_species, ",") # one row has 13 species

dfEvaluate_Tidy <- dfEvaluate_Tidy %>% 
  separate(additional_species, 
           into = paste0("Species", 1:13), 
           sep = ",", 
           extra = "merge")

# The additional species columns can be used later when comparing taxonomy found.

# Export dfEvaluate_Tidy as a .tsv file for cleaner raw data. Code commented out as it does not need to be redone.
#write_tsv(dfEvaluate_Tidy, file = "Evaluation_Abundance_Table.tsv", col_names = T)



### Import and clean the taxonomy file from the sequence source. https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&page_size=10&acc=ERR2744174&display=analysis

# Import the merged taxa file.
dfSource <- read_csv("bquxjob_bf4b274_181ed506074.csv")
View(dfSource)

# # Pivot table based on taxonomy rank.
# test <- dfSource
# test <- test %>%
#   pivot_wider(names_from = rank, values_from = name, values_fill = "NA")


# Create a separate dataframe for Family level data for the evaluation data.
dfEvaluate_Family <- dfEvaluate_Tidy[,5:9]

# Remove rows with NAs in Family
dfEvaluate_Family <- dfEvaluate_Family[complete.cases(dfEvaluate_Family$Family),]

# Keep rows with NAs in both Genus and Species - we only want these rows as the abundance information encompasses the whole Family.
dfEvaluate_Family <- dfEvaluate_Family[!complete.cases(dfEvaluate_Family$Genus),]
dfEvaluate_Family <- dfEvaluate_Family[,c(1, 5)]

# Convert abundance information from character to numerical values.
dfEvaluate_Family[,2] <- as.numeric(dfEvaluate_Family[,2])
sapply(dfEvaluate_Family, class)

# Sum the abundance for each sample.
round(sum(dfEvaluate_Family[, 2]), 3)

# Fix the rownames.
rownames(dfEvaluate_Family) <- 1:nrow(dfEvaluate_Family)


# Create a separate dataframe for Family level data for the source data.
dfSource %>%
  count(rank)
# 73 family level

dfSource_Family <- dfSource %>%
  filter(rank == "family")

# Number of unique Family names.
length(unique(dfSource_Family$name)) # all unique

# Convert count data to relative abundance.
dfSource_Family[, 5:6] <- sapply(dfSource_Family[, 5:6], make_relative)

dfSource_Family[, 5:6] <- sapply(dfSource_Family[, 5:6], function(x) round(x, 6))

# Use total count column for comparison with MetaPhlAn output.
dfSource_Family <- dfSource_Family[, 4:5]

# Reorder based on percentage.
dfSource_Family <- dfSource_Family %>%
  arrange(desc(total_count))


## Create a Venn Diagram showing the overlapping taxa detected. Data needs to be in a list form.
list_Family <- list(Source = dfSource_Family$name,
                    Pipeline = dfEvaluate_Family$Family)

ggVennDiagram(list_Family,
              label_alpha = 0,
              edge_size = 0.05) +
  scale_fill_gradient(low = "#cfebe9", high = "#2ad4c3") +
  ggtitle("Venn diagram for Family level comparison") +
  easy_center_title()
  
# Top five Family taxa found in each dataframe.
dfSource_Family[1:5,] ; dfEvaluate_Family[1:5,]

# Same top five.
#https://academic.oup.com/dnaresearch/article/26/5/391/5541856
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8317196/pdf/TSM2-4-174.pdf

# The intersecting taxa.
intersect(dfSource_Family$name, dfEvaluate_Family$Family)

# Additional taxa found in source dataframe.
diff_family <- setdiff(dfSource_Family$name, dfEvaluate_Family$Family)
dfDiff_Family <- dfSource_Family %>%
  filter(name %in% diff_family)

# All the additional families selected are less than 0.3% of relative abundance.



### Do the same for Phylum.
# Create a separate dataframe for Phylum level data for the evaluation data.
dfEvaluate_Phylum <- dfEvaluate_Tidy[,2:9]

# Remove rows with NAs in Phylum
dfEvaluate_Phylum <- dfEvaluate_Phylum[complete.cases(dfEvaluate_Phylum$Phylum),]

# Keep rows with NAs in Class - we only want these rows as the abundance information encompasses the whole Phylum.
dfEvaluate_Phylum <- dfEvaluate_Phylum[!complete.cases(dfEvaluate_Phylum$Class),]
dfEvaluate_Phylum <- dfEvaluate_Phylum[,c(1, 8)]

# Convert abundance information from character to numerical values.
dfEvaluate_Phylum[,2] <- as.numeric(dfEvaluate_Phylum[,2])
sapply(dfEvaluate_Phylum, class)

# Sum the abundance for each sample.
round(sum(dfEvaluate_Phylum[, 2]), 3)

# Fix the rownames.
rownames(dfEvaluate_Phylum) <- 1:nrow(dfEvaluate_Phylum)


# Create a separate dataframe for Phylum level data for the source data.
dfSource %>%
  count(rank)
# 12 Phylum level

dfSource_Phylum <- dfSource %>%
  filter(rank == "phylum")

# Number of unique Phylum names.
length(unique(dfSource_Phylum$name)) # all unique

# Convert count data to relative abundance.
dfSource_Phylum[, 5:6] <- sapply(dfSource_Phylum[, 5:6], make_relative)

dfSource_Phylum[, 5:6] <- sapply(dfSource_Phylum[, 5:6], function(x) round(x, 6))

# Use total count column for comparison with MetaPhlAn output.
dfSource_Phylum <- dfSource_Phylum[, 4:5]

# Reorder based on percentage.
dfSource_Phylum <- dfSource_Phylum %>%
  arrange(desc(total_count))


## Create a Venn Diagram showing the overlapping taxa detected. Data needs to be in a list form.
list_Phylum <- list(Source = dfSource_Phylum$name,
                    Pipeline = dfEvaluate_Phylum$Phylum)

ggVennDiagram(list_Phylum,
              label_alpha = 0,
              edge_size = 0.05) +
  scale_fill_gradient(low = "#cfebe9", high = "#2ad4c3") +
  ggtitle("Venn diagram for Phylum level comparison") +
  easy_center_title()

# Top five Phylum taxa found in each dataframe.
dfSource_Phylum[1:5,] ; dfEvaluate_Phylum[1:5,]

# Same top four; the sixth Phylum corresponds to the fifth in dfEvaluate_Phylum
#https://academic.oup.com/dnaresearch/article/26/5/391/5541856
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8317196/pdf/TSM2-4-174.pdf

# The intersecting taxa.
intersect(dfSource_Phylum$name, dfEvaluate_Phylum$Phylum)

# Additional taxa found in source dataframe.
diff_phylum <- setdiff(dfSource_Phylum$name, dfEvaluate_Phylum$Phylum)
dfDiff_Phylum<- dfSource_Phylum %>%
  filter(name %in% diff_phylum)

# All the additional families selected are less than 0.03% of relative abundance.


# I also used all the additional species column and conducted a Venn diagram again. Produced the same result so I deleted the code.


