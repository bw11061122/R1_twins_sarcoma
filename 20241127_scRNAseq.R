###################################################################################################################################
# SCRIPT 8

# Script to remove germline mutations in copy number-altered regions
# 2024-11-27
# Barbara Walkowiak bw18

# INPUT: 
# 1 PILEUP dt (.cnv.segment) for each sample 

# OUTPUT:
# 1 list of mutations to exclude based on presence in copy number germline variants 

# Script to identify regions of copy number alterations (deletions or duplications) in the germline
# Identify mutations in those regions where the VAF is 0.3 and so those can be explained as a germline copy number alteration 

###################################################################################################################################
# LIBRARIES
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

###################################################################################################################################
# INPUT DATA 

setwd('/Users/bw18/Desktop/1SB')

# Read the input data 
twins.data = Read10X(data.dir = "scRNAseq/NB13760628/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
twins = CreateSeuratObject(counts = data, project = "twins.sarcoma", min.cells = 3, min.features = 200)
twins

###################################################################################################################################
# BASIC QC AND PRE-PROCESSING

# Examine a few genes in this set (select 3 genes, look across the first 30 cells), this uses the direct output of Read10X
twins.data[c("CD14", "CD68", "PAX7"), 1:30] # active assay: RNA (3 features, 0 variable features)
twins.data[c("CD3D", "TCL1A", "MS4A1"), 1:30] # none of the features provided found in this assay 

# Okay let's try and do QC now 
head(twins@meta.data)
# orig.ident - this is the name of out project 
# nCount_RNA - this is the number of RNA counts 
# nFeature_RNA - this is the number of features that these counts map to (ie genes)

# add a column that shows content of mitochondria
# for this package, we add a column using '[[]]'
twins[["percent.mt"]] <- PercentageFeatureSet(twins, pattern = "^MT-")
twins@meta.data %>% head() # added mitochondrial content correctly 

# QC plots one can obtain 
VlnPlot(twins, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# standard QC plots to do 
# nFeature_RNA number of genes / transcripts identified 
# nCount_RNA number of mapped reads
# percent.mt proportion of mitochondrial reads 

# Based on this, cut-off of 15% should be reasonable 
twins = subset(twins, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

# Data normalization
# from GitHub tutorial: By default, we employ a global-scaling normalization method “LogNormalize” 
# that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result
# There is a different method implemented in Seurat: SCTransform (which one is better and why?)
twins = NormalizeData(twins, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature selection 
twins = FindVariableFeatures(twins, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(twins), 10)

# plot variable features with and without labels
plot1 = VariableFeaturePlot(twins) # quite a lot of the variable count - do we believe this?
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE)
# what it picks up: COL8A1, MT1G, ELN, HBB, PTGDS, POSTN, COL11A1, HBA2, HBA1 
# HBA1/2 are a bit awkward - this is literally encoding for haemoglobin?
# COL11A1 / COL8A1 - collagen genes 

# Scaling the data 
# apply a linear transformation (‘scaling’), standard pre-processing before dimensional reduction.
# we do this so that: (1) mean expression of each gene across all cells is 0
# (2) the variance in expression across cells is 1 
# in this way, the standardized data has mean expression 0 and variance 1
# therefore, highly expressed genes do not dominate all the data

all.genes = rownames(twins)
twins = ScaleData(twins, features = all.genes)

twins[["RNA"]]@counts[1:6, 1:6]
twins[["RNA"]]@data[1:6, 1:6]


