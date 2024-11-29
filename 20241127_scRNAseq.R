###################################################################################################################################
# SCRIPT 8

# Script to remove germline mutations in copy number-altered regions
# 2024-11-27
# Barbara Walkowiak bw18

# INPUT: 
# 1 scRNAseq CellRanger output 

# OUTPUT: 
# 1 basic scRNA-seq analysis (just teaching myself how to do this)

# Note: I am getting started by following this tutorial: 
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_2.html
# plus some adjustments (as I am using Seurat V5, cf this workshop is for Seurat V3)

###################################################################################################################################
# QUESTIONS AND THOUGHTS 

# There are three files as each was sequenced on a different lane
# How do I aggregate samples from each lanes? # should I somehow merge the counts? 
# Based on GitHub here: https://github.com/satijalab/seurat/discussions/4197, you can but you would need to go back to fastq - you cannot do this on the matrix data
# An alternative (which they seem to be against here?) is to use Seurat integrate 
# Seurat v5 enables integration via IntegrateLayers function - is this appropriate to use in this case?

# Do I check whether there are differences between runs (lanes)? > likely yes 
# I imagine this could be done by checking if samples from different runs cluster on a PCA etc.

# Should I exclude cells with too many counts (UMIs) or features? e.g., possible doublets?

###################################################################################################################################
# LIBRARIES
library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)

###################################################################################################################################
# INPUT DATA 

setwd('/Users/bw18/Desktop/1SB/scRNAseq')

# Samples 
# I have the following samples (these are tumour samples!)
# PD62341 fresh: CG_SB_NB13652544, CG_SB_NB13652545, CG_SB_NB13652546
# PD63383 fresh: CG_SB_NB13760628, CG_SB_NB13760629, CG_SB_NB13760630
# PD63383 frozen (nuclear): CG_SB_NB14599068, CG_SB_NB14599069, CG_SB_NB14599070

# Read the input data 
# PD62341 tumour 
tumourPD62341.1 = Read10X(data.dir = "data/CG_SB_NB13652544/filtered_feature_bc_matrix/") # lane 1
tumourPD62341.2 = Read10X(data.dir = "data/CG_SB_NB13652545/filtered_feature_bc_matrix/") # lane 2
tumourPD62341.3 = Read10X(data.dir = "data/CG_SB_NB13652546/filtered_feature_bc_matrix/") # lane 3

# PD63383 tumour 
tumourPD63383.1 = Read10X(data.dir = "data/CG_SB_NB13760628/filtered_feature_bc_matrix/") # lane 1
tumourPD63383.2 = Read10X(data.dir = "data/CG_SB_NB13760629/filtered_feature_bc_matrix/") # lane 2
tumourPD63383.3 = Read10X(data.dir = "data/CG_SB_NB13760630/filtered_feature_bc_matrix/") # lane 3

# PD62341 tumour (frozen, nuclear)
tumourPD63383.1n = Read10X(data.dir = "data/CG_SB_NB14599068/filtered_feature_bc_matrix/") # lane 1
tumourPD63383.2n = Read10X(data.dir = "data/CG_SB_NB14599069/filtered_feature_bc_matrix/") # lane 2
tumourPD63383.3n = Read10X(data.dir = "data/CG_SB_NB14599070/filtered_feature_bc_matrix/") # lane 3

###################################################################################################################################
# BASIC QC AND PRE-PROCESSING: trial for one sample

# run manual analysis on one sample to decide what kinds of thresholds are based to use for nr counts, nr features and mt content

# read the data into the Seurat object 
tPD62341.1 = CreateSeuratObject(counts = tumourPD62341.1, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)

# calculate the percentage mt content 
tPD62341.1[["percent.mt"]] = PercentageFeatureSet(tPD62341.1, pattern = "^MT-")

# plot the distribution of the mt content in each sample 
# 10 sounds like a reasonable threshold and should give you most of the cells of interest 
pdf('Results/20241127_hist_mt_content_PD62341tumour.pdf', height = 4.5, width = 4.5)
hist(tPD62341.1[["percent.mt"]] %>% unlist(), xlab = '% mt content', main = 'PD62341 tumour', breaks = 100)
abline(v = 10, col = 'red', lwd = 2.5)   
dev.off()

# add CellID to the metadata 
tPD62341.1@meta.data$CellID = rownames(tPD62341.1@meta.data) # add the cellID to the metadata
meta = tPD62341.1@meta.data # create a dt with all the metadata of interest, including cell ID 

# plot histograms of nCount_RNA_min and nFeature_RNA_min to decide the thresholds
# usually (from a few papers I looked at), cells with fewer than 300 genes or fewer than 1000 UMIs are removed 

# nCountRNA is the total number of molecules detected in each cell 
pdf('Results/20241127_hist_nCountRNA_PD62341tumour.pdf', height = 4, width = 8)
ncounts_rna = as.numeric(tPD62341.1[["nCount_RNA"]] %>% unlist())
par(mfrow=c(1,2))
hist(ncounts_rna[ncounts_rna < 2500], xlab = 'nCountRNA', main = 'PD62341 tumour', breaks = 100)
hist(ncounts_rna, xlab = 'nCountRNA', main = 'PD62341 tumour', breaks = 100)
dev.off() # 1000 sounds reasonable based on the hist 

# nFeatureRNA is the number of genes detected in each cell 
pdf('Results/20241127_hist_nFeatureRNA_PD62341tumour.pdf', height = 4, width = 8)
nfeatures_rna = as.numeric(tPD62341.1[["nFeature_RNA"]] %>% unlist())
par(mfrow=c(1,2))
hist(nfeatures_rna[ncounts_rna < 2000], xlab = 'nFeatureRNA', main = 'PD62341 tumour', breaks = 100)
hist(nfeatures_rna, xlab = 'nFeatureRNA', main = 'PD62341 tumour', breaks = 100)
dev.off() # 300 sounds reasonable based on the hist 

# we can also plot those 3 QC parameters and save them to results files
pdf('Results/20241127_vlnPlotQC_PD62341tumour.pdf', height = 4, width = 8)
par(mfrow=c(1,1))
VlnPlot(tPD62341.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off() # I don't see how those plots would help anyone choose the lower threshold for counts and features 

# which cells to keep?
keeprows = meta[meta$nCount_RNA > 1000,] # require number of UMIs to be higher than minimum (1000 UMIs)
keeprows = keeprows[keeprows$nFeature_RNA > 300,] # require number of features to be higher than minimum (300 genes)
keeprows = keeprows[keeprows$percent.mt <= 10,] # exclude cells with mitochondrial content above a threshold (10%)
cellids_keep = rownames(keeprows) # which CellIDs to keep 
tPD62341.1_filtered = subset(tPD62341.1, subset = CellID %in% cellids_keep) # select cells of interest

# which genes to keep?
keep_features=rownames(tPD62341.1_filtered) 
tPD62341.1_filtered <- subset(tPD62341.1_filtered, features=keep_features)
# I could exclude heat shock or ribosomal genes - variable takes in literature as to whether one should do this or not?
# I won't do this in this dataset but I will run separate analysis excluding those genes

# run standard Seurat pipeline 
tPD62341.1_filtered = NormalizeData(object = tPD62341.1_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
tPD62341.1_filtered = FindVariableFeatures(tPD62341.1_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
tPD62341.1_filtered = ScaleData(tPD62341.1_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
tPD62341.1_filtered = RunPCA(tPD62341.1_filtered, npcs = 70, features=VariableFeatures(object = tPD62341.1_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)

# Determine how many PCs to include in downstream analysis 
pdf('Results/20241127_elbowPlot_PD62341tumour.pdf', height = 4, width = 8)
ElbowPlot(tPD62341.1_filtered, ndims = 50) # use 50 PCs 
dev.off()

# plot variance explained (scree plot)
mat <- tPD62341.1_filtered[["RNA"]]$scale.data 
pca <- tPD62341.1_filtered[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2 
varExplained = eigValues / total_variance
varExplained %>% enframe(name = "PC", value = "varExplained" ) %>%
  ggplot(aes(x = PC, y = varExplained)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  ggtitle("scree plot, PD62341 tumour") 
ggsave(glue('Results/20241127_varianceExplained_tumourPD62341_1.pdf'), height = 4, width = 4)

# Clustering 
tPD62341.1_filtered = FindNeighbors(tPD62341.1_filtered, dims = 1:20, verbose = T)
tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = 0.6, verbose = T)
tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)
tPD62341.1_filtered = RunTSNE(tPD62341.1_filtered) # might as well run it but tSNE doesn't seem to be popular 

# Plot outcomes of clustering 
pdf('Results/20241127_umap_PD62341tumour.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "umap", label = TRUE)
dev.off() # that actually looks reasonable!
pdf('Results/20241127_tsne_PD62341tumour.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "tsne", label = TRUE)
dev.off()
pdf('Results/20241127_pca_PD62341tumour.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "pca", label = TRUE)
dev.off()

# Can I try clustering with different resolution parameter
for (res in seq(0.1, 1.2, by = 0.1)){
  tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = res, verbose = T)
  tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)
  pdf(glue('Results/20241127_umap_PD62341tumour_{res}_res.pdf'), height = 4, width = 4)
  plot = DimPlot(tPD62341.1_filtered, reduction = "umap", label = TRUE)
  print(plot)
  dev.off() 
} # After this I'm quite inclined to go with resolution 0.2-0.4 

# Identify features and look at markers 
tPD62341.1_filtered.markers <- FindAllMarkers(tPD62341.1_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# include genes that take part in the fusion, early embryonic, CD markers, synaptophysin, stuff Nathan looked at (to see if I am getting sth vaguely similar)
genes_of_interest = c('MYOD1', 'PAX3', 'PAX7', 'CYP11A1', 'CD14', 'CD68', 'SYP', 'MN1', 'ZNF341', 'NR5A1', 'KCNQ1', 'TP53') # MYOG is 0 for all so I got rid of this
pdf(glue('Results/20241127_featurePlot_PD62341tumour.pdf'), width=12, height=8)
FeaturePlot(tPD62341.1_filtered, features = genes_of_interest, order=TRUE)
dev.off()

# Heatmap of markers 
top10 = tPD62341.1_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(glue('Results/20241127_heatmapTop10_PD62341tumour.pdf'), width=12, height=8)
DoHeatmap(tPD62341.1_filtered, features = top10$gene) + NoLegend()
dev.off()

# Dotplot of genes of interest 
pdf(glue('Results/20241127_dotPlot_PD62341tumour_GOI.pdf'), width=10, height=6)
DotPlot(object = tPD62341.1_filtered, features = genes_of_interest)
dev.off()

# Dotplot of cluster markers

# Score by cell cycle expression (S and G2M genes accessed according to https://satijalab.org/seurat/articles/cell_cycle_vignette.html from October 2023)
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes
tPD62341.1_filtered = CellCycleScoring(tPD62341.1_filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# set.ident = TRUE sets the identity of the Seurat object to the cell-cycle phase

# Color UMAP by cell cycle phase (based on expression of cell cycle genes)
pdf('Results/20241127_umap_PD62341tumour_cellCycle.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "umap", label = TRUE, group.by = 'Phase')
dev.off() # don't see a separate clusters of proliferating cells :(
pdf('Results/20241127_pca_PD62341tumour_cellCycle.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "pca", label = TRUE, group.by = 'Phase')
dev.off() # kind of separate but not very clear? not sure if I should do anything about this 
# I guess that means that we don't need to regress out the variation related to cell cycle phase 

# Run analysis and clustering on dataset with removed ribosomal and heatshock genes
genes_to_exclude = read.csv('Data/excludeGenes.tsv', sep='\t', header = F) %>% unlist()
keep_features2=rownames(tPD62341.1_filtered)[which(!rownames(tPD62341.1_filtered) %in% genes_to_exclude)]

tPD62341.1_filtered2 <- subset(tPD62341.1_filtered, features=keep_features2)
tPD62341.1_filtered2 = NormalizeData(object = tPD62341.1_filtered2, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
tPD62341.1_filtered2 = FindVariableFeatures(tPD62341.1_filtered2, selection.method = "vst", nfeatures = 2000) # default 2k 
tPD62341.1_filtered2 = ScaleData(tPD62341.1_filtered2, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
tPD62341.1_filtered2 = RunPCA(tPD62341.1_filtered2, npcs = 70, features=VariableFeatures(object = tPD62341.1_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
tPD62341.1_filtered2 = FindNeighbors(tPD62341.1_filtered2, dims = 1:20, verbose = T)
tPD62341.1_filtered2 = FindClusters(tPD62341.1_filtered2, resolution = 0.6, verbose = T)
tPD62341.1_filtered2 = RunUMAP(tPD62341.1_filtered2, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)
tPD62341.1_filtered2 = RunTSNE(tPD62341.1_filtered2) # might as well run it but tSNE doesn't seem to be popular 
tPD62341.1_filtered2 = CellCycleScoring(tPD62341.1_filtered2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Plot outcomes of clustering (it doesn't look massively different to me!)
pdf('Results/20241127_umap_PD62341tumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "umap", label = TRUE)
dev.off() # that actually looks reasonable!
pdf('Results/20241127_tsne_PD62341tumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "tsne", label = TRUE)
dev.off()
pdf('Results/20241127_pca_PD62341tumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "pca", label = TRUE)
dev.off()
pdf('Results/20241127_umap_PD62341tumour_rmGenes_cellCycle.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "umap", label = TRUE, by='Phase')
dev.off() # that actually looks reasonable!

###################################################################################################################################

# A function that takes the Seurat Object and carries out basic pre-processing, clustering etc.
sc_preprocessing = function(data, sample_name, mt_threshold, feature_min, count_min, n_pcs, n_clusters, resolution, min_dist, nn){
  
  # Required arguments
  # data = obtained with Read10X function
  # sample_name = name of the sample to be used in plots etc.
  # mt_threshold = % of mt reads above which a cell is excluded
  
  # Initialize Seurat objects with the raw (non-normalized data)
  # Optional settings when reading in the file: 
  # min.cells = Include features detected in at least this many cells
  # min.features = Include cells where at least this many features are detected
  # I guess this can be a useful filtering step but we can leave that for now and do it w/ default (not setting more stringent thresholds)
  seurat = CreateSeuratObject(counts = data, project = "twins.sarcoma",names.field = 1, names.delim = "_", meta.data = NULL)
  
  # QC: calculate fraction of mitochondrial reads (high indicates poor-quality cells)
  # plot histogram to see if the threshold used makes sense
  seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = "^MT-")
  
  pdf(glue('Results/20241127_hist_mt_content_{sample_name}.pdf'), height = 4.5, width = 4.5)
  hist(tPD62341.1[["percent.mt"]] %>% unlist(), xlab = '% mt content', main = sample_name, breaks = 100)
  abline(v = mt_threshold, col = 'red', lwd = 2.5)   
  dev.off()
  
  seurat@meta.data$CellID = rownames(seurat@meta.data) # add the cellID to the metadata
  meta = seurat@meta.data # create a dt with all the metadata of interest, including cell ID 
  
  keeprows = meta[meta$nCount_RNA > nCountRNA_min,] # require number of UMIs to be higher than minimum 
  keeprows = keeprows[keeprows$nFeature_RNA > nFeature_RNA_min,] # require number of features to be higher than minimum 
  keeprows = keeprows[keeprows$percent.mt <= mt_threshold,] # exclude cells with mitochondrial content above a threshold
  
  cellids_keep = rownames(keeprows)
  seurat_filtered = subset(seurat, subset = CellID %in% cellids_keep) # select cells of interest 
  
  #Exclude Ribosomal, Mitochondrial and other unimportant Genes
  keep_features=rownames(data_seurat)[which(!rownames(data_seurat)%in%excludeGenes)]
  data_filtered <- subset(data_filtered, features=keep_features)
  data_filtered = NormalizeData(object = data_filtered, normalization.method = "LogNormalize", scale.factor = 1e4)
  data_filtered = FindVariableFeatures(data_filtered, selection.method = "vst", nfeatures = 2000)
  data_filtered = ScaleData(data_filtered, verbose = T) #Memory issues if we scale using all genes. 
  data_filtered = RunPCA(data_filtered, npcs = npcs1, features=VariableFeatures(object = data_filtered), verbose = F)
  
  data_filtered = FindNeighbors(data_filtered, dims = 1:npcs1, verbose = T)
  data_filtered = FindClusters(data_filtered, resolution = clustres, verbose = T)
  data_filtered = RunUMAP(data_filtered, dims = 1:npcs1, verbose = T, min.dist = finalMinDist, n.neighbors = finalNN)
  
  s.genes = cc.genes.updated.2019$s.genes
  g2m.genes = cc.genes.updated.2019$g2m.genes
  data_filtered = CellCycleScoring(data_filtered, s.features = s.genes, 
                                   g2m.features = g2m.genes, set.ident = TRUE)
  
  }

# run for all your samples
sc_preprocessing(tumourPD62341.1)

# at what point should I integrate samples from different lanes?

# how do I do cell type annotation?


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

# to access count data
# note: there has been a change from Seurat v4 to Seurat v5
# therefore, the previous syntax (used in Seurat v4) will not work anymore
# V4 syntax: twins[["RNA"]]@counts[1:6, 1:6]
# V4 syntax: twins[["RNA"]]@data[1:6, 1:6]

FetchData(object = twins, vars = c('HBA1', 'COL11A1'), layer = "counts")
# counts: raw count data, before normalization or any pre-processing 
FetchData(object = twins, vars = c('HBA1', 'COL11A1'), layer = "data")
# data: library size-normalized, log-transformed data 
FetchData(object = twins, vars = c('HBA1', 'COL11A1'), layer = "scale.data")
# scale.data: scaled data (we did the linear transformation here)

# scaling is only done on genes that will be used as input to PCA / clustering
# as default, Seurat only performs scaling on previously identified 2k variable features
# we want to make sure that mt-mapped genes are not included in the analysis 
twins = ScaleData(twins, vars.to.regress = "percent.mt")

###################################################################################################################################
# PCA AND FEATURE SELECTION 

# we can do PCA, which is classic dimensionality reduction based on linear transformation
twins = RunPCA(twins, features = VariableFeatures(object = twins), verbose = FALSE)
p1 = DimPlot(twins, reduction = "pca")
# looks like you may want to do more filtering on this going forward 

twins[['pca']]
# this provides you with an output of the general idea of how this was done

# what are the available methods for dimensionality reduction?
utils::methods(class = 'DimReduc')

# if we want to see the coordinates after doing the PCA projection
Embeddings(twins, reduction = "pca") %>% head()

# replicate the plot using ggplot2 function
p2 = Embeddings(twins, reduction = "pca") %>% 
  as.data.frame()%>% 
  ggplot(aes(x = PC_1, y = PC_2)) +
  geom_point(color = "red", size = 0.5) +
  theme_classic()
p2

# How many PCs should we include for downstream analysis?

# what is the idea behind this: Seurat clusters cells based on the PCA scores 
# each PC represents a 'meta-feature' that combines information across clusters
# top PCs represent a 'compression' of the dataset that extracts the most important info

# what is JackStraw?
# the idea is to permute a subset of data (1%) and rerun the PCA (we do this multiple times)
# using this test, we can identify significant PCs as those which have a strong enrichment of low p-values features
twins = JackStraw(twins, num.replicate = 100, dims = 50)
twins = ScoreJackStraw(twins, dims = 1:50)
JackStrawPlot(twins, dims = 1:30) # what are we expecting to see? 

# alternatively, we can use an Elbow plot
# the Elbow plot ranks PCs based on the percentage of variance explained by each one (like a scree plot!)
ElbowPlot(twins, ndims = 50) 
# looks like it falls off after 10-20
# the y axis here is STANDARD DEVIATION, not variance! 
# would SD or variance be more informative / appropriate in this case?

# code to reconstruct the scree plot:
mat = twins[["RNA"]]$scale.data # select scaled data ($ works for Seurat V5, see here: https://satijalab.org/seurat/articles/seurat5_essential_commands
pca = twins[["pca"]] # select the PCA coordinates 

# Get the total variance explained by each PC:
total_variance = sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

varExplained %>% enframe(name = "PC", value = "varExplained" ) %>%
  ggplot(aes(x = PC, y = varExplained)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  ggtitle("scree plot")

###################################################################################################################################
# CLUSTERING

# Great - the way this is doing clustering is basically KNN
# we embed cells in a graph structure (KNN)
# edges are drawn between cells of similar expression patterns
# we then try to partition this graph into 'communities' of connected cells
# the KNN is constructed based on the Euclidean distance in PCA space
# we then refine edge weights using Jaccard similarity 
twins = FindNeighbors(twins, dims = 1:20)

# next, cells are clustered using he Louvain algorithm or SLM
# these algorithms iteratively group cells together and the idea is that we optimize a specific function
# this is implemented in the 'FindClusters' function (we can adjust granularity with 'resolution' parameter)
twins = FindClusters(twins, resolution = 0.5) 
# resolution depends on the number of cells you have, you should check if the resolution you are using here is appropriate

# ID clusters for the first 5 cells
head(Idents(twins), 5) # you have 10 clusters, the first 5 cells are in the 1st and 3rd clusters

###################################################################################################################################
# NON-LINEAR DIMENSIONALITY REDUCTION (UMAP / tSNE)

# you can run both, tho please consider if you believe any of this?
twins = RunUMAP(twins, dims = 1:20)
twins = RunTSNE(twins, dims = 1:20)

# plot your UMAP and tSNE results
DimPlot(twins, reduction = "umap", label = TRUE)
DimPlot(twins, reduction = "tsne", label = TRUE)

# compare to the PCA output
DimPlot(twins, reduction = "pca") # clearly quite different to what UMAP is telling you

# we can plot the UMAP visualization results using ggplot2 as well
# extract cell embedding obtained with UMAP (the first two UMAP dimensions)
Embeddings(twins, reduction = "umap") %>% head() # you can see umap_1 and umap_2 coordinates here

# add coloring by cluster
umap_df = bind_cols(seurat_clusters = twins@meta.data$seurat_clusters, # access clusters from the meta.data
                     Embeddings(twins, reduction = "umap") %>% as.data.frame()) # access UMAP embeddings
ggplot(umap_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = seurat_clusters), size = 0.5) +
  theme_classic(base_size = 14) 

###################################################################################################################################
# Identify cluster biomarkers (differentially expressed genes)

# find markers of the 1st cluster
cluster1.markers <- FindMarkers(twins, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5) 

# find markers that differentiate cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(twins, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5) # ident.2 sets what you want to distinguish your cells from

# find markers for each cluster
twins.markers <- FindAllMarkers(twins, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# why do we set this kind of lfc threshold?

# we want to identify top markers for each cluster
twins.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# we can visualize genes that come up as top markers
VlnPlot(twins, features = c("EXOC3L2", "DIO3OS"))
# these are very high in cluster 0 and there are some reads in other clusters 

# data slots in the matrix
twins[["RNA"]]$counts[c("EXOC3L2", "DIO3OS"), 1:30]
twins[["RNA"]]$data[c("EXOC3L2", "DIO3OS"), 1:30]
# you can look at the normalized (not scaled) counts for each of these genes in the top 30 cells 

# create a feature plot 
FeaturePlot(twins, features = c('MYOD1', 'PAX3', 'PAX7', 'CD14', 'CD68', 'MYOG', 'ARID1A', 'PIK3CA'))
# does not find MYOD1 or MYOG - does it mean they are not there at all, i.e., were 0 in this tumour sample in all cells?

# Reorder data points, such that the data points where the gene is expressed are on top
p = FeaturePlot(twins, features = "CD14") # let's say we want to look at CD14 gene
p_after = p
p_after$data <- p_after$data[order(p_after$data$CD14),]
CombinePlots(plots = list(p, p_after))

###################################################################################################################################
# Expression heatmap for different clusters 
top10 <- twins.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(twins, features = top10$gene) + NoLegend()
# omg this looks like it did actually do something that is reasonably correct 

###################################################################################################################################
# Cell type annotation 

# Logic followed in the method that I am going to be using here:
# compare your own cell types to cell types in reference atlases, such at the HCA
# the method used in Seurat aims to identify 'anchors' between pairs of datasets
# 'anchors' represent pairwise correspondences b/n individual cells (one in each dataset)
# we hypothesize that each of those cells originate from the same biological state 
# anchors are used to harmonize datasets and transfer info from one dataset to another

# where do I obtain a suitable reference to do this for? need to figure out what to download from HCA I think 




