#### Script for Microbiome Taxonomy Analysis

## Author: Shalvi Chirmade
## BINF*6999 Research Project Summer 2022

## Van Raay Lab at the University of Guelph

## Primary Advisor: Dr. Terry Van Raay
## Secondary Advisor: Dr. Lewis Lukens

## Description: This is a reproducible script for conducting microbiome taxonomy analysis. It allows for the creation of a taxonomy stacked barplot, PCA analysis and hierarchical clustering. 

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

















