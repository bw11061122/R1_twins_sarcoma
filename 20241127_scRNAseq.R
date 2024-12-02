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

# Genes that can be excluded from clustering (e.g., ribosomal genes)
genes_to_exclude = read.csv('Data/excludeGenes.tsv', sep='\t', header = F) %>% unlist()

###################################################################################################################################
# BASIC QC AND PRE-PROCESSING: trial for one sample

# manual analysis on one sample to figure out thresholds and see how this works 

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
hist(nfeatures_rna[nfeatures_rna < 2000], xlab = 'nFeatureRNA', main = 'PD62341 tumour', breaks = 100)
hist(nfeatures_rna, xlab = 'nFeatureRNA', main = 'PD62341 tumour', breaks = 100)
dev.off() # 300 sounds reasonable based on the hist 

# we can also plot those 3 QC parameters and save them to results files
pdf('Results/20241127_vlnPlotQC_PD62341tumour.pdf', height = 4, width = 8)
par(mfrow=c(1,1))
VlnPlot(tPD62341.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off() # I don't see how those plots would help anyone choose the lower threshold for counts and features 

# plot QC features against one another 
pdf('Results/20241127_featureScatter_PD62341tumour.pdf', height = 4, width = 12)
plot1 = FeatureScatter(tPD62341.1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(tPD62341.1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# which cells to keep?
keeprows = meta[meta$nCount_RNA > 1000 & meta$nCount_RNA < 30000,] # require number of UMIs to be higher than minimum (1000 UMIs)
keeprows = keeprows[keeprows$nFeature_RNA > 300 & keeprows$nFeature_RNA < 8000,] # require number of features to be higher than minimum (300 genes)
keeprows = keeprows[keeprows$percent.mt <= 10,] # exclude cells with mitochondrial content above a threshold (10%)
cellids_keep = rownames(keeprows) # which CellIDs to keep 
tPD62341.1_filtered = subset(tPD62341.1, subset = CellID %in% cellids_keep) # select cells of interest

# which genes to keep?
keep_features=rownames(tPD62341.1_filtered) 
tPD62341.1_filtered = subset(tPD62341.1_filtered, features=keep_features)
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
mat = tPD62341.1_filtered[["RNA"]]$scale.data 
pca = tPD62341.1_filtered[["pca"]]
total_variance = sum(matrixStats::rowVars(mat))
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
tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = 0.4, verbose = T)
tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.3, n.neighbors = 30)
tPD62341.1_filtered = RunTSNE(tPD62341.1_filtered) # might as well run it but tSNE doesn't seem to be popular 

# Plot outcomes of clustering 
pdf('Results/20241127_umap_PD62341tumour.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "umap", label = TRUE)
dev.off() 
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
tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = 0.4, verbose = T)
tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.3, n.neighbors = 30) # with preferred settings
tPD62341.1_filtered.markers = FindAllMarkers(tPD62341.1_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# include genes that take part in the fusion, early embryonic, CD markers, synaptophysin, stuff Nathan looked at (to see if I am getting sth vaguely similar)
genes_of_interest = c('MYOD1', 'PAX3', 'PAX7', 'CYP11A1', 'CD14', 'CD68', 'SYP', 'MN1', 'ZNF341', 'NR5A1', 'KCNQ1', 'TP53') # MYOG is 0 for all so I got rid of this
pdf(glue('Results/20241127_featurePlot_PD62341tumour.pdf'), width=12, height=8)
FeaturePlot(tPD62341.1_filtered, features = genes_of_interest, order=TRUE)
dev.off()

# Heatmap of markers 
top10 = tPD62341.1_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(glue('Results/20241127_heatmapTop10_PD62341tumour.pdf'), width=12, height=10)
DoHeatmap(tPD62341.1_filtered, features = top10$gene) + NoLegend()
dev.off()

# Dotplot of cluster markers (4 for each cluster)
tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = 0.4, verbose = T)
tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.3, n.neighbors = 30) # with preferred settings
top4 = tPD62341.1_filtered.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
pdf(glue('Results/20241127_dotPlotTop4_PD62341tumour.pdf'), width=10, height=8)
DotPlot(tPD62341.1_filtered, features = top4$gene) + coord_flip()
dev.off()

# Dotplot of genes of interest 
pdf(glue('Results/20241127_dotPlot_PD62341tumour_GOI.pdf'), width=10, height=6)
DotPlot(object = tPD62341.1_filtered, features = genes_of_interest)
dev.off()

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
excludeGenes='^(EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+_',rownames(tPD62341.1_filtered),value=TRUE), # poorly characterized (what does this mean??)
                            grep('MALAT1',rownames(tPD62341.1_filtered),value=TRUE), # contamination - why are we excluding this? how do we know this is contamination
                            grep('^HB[BGMQDAZE][12]?_',rownames(tPD62341.1_filtered),value=TRUE), # contamination (apparently you should filter out blood markers)
                            grep('^MT-',rownames(tPD62341.1_filtered),value=TRUE), # mt genes  
                            grep(excludeGenes,rownames(tPD62341.1_filtered),value=TRUE) # housekeeping genes
))
genes_to_exclude = c(genes_to_exclude, coreExcludeGenes)
paste("Number of genes to exclude from analysis:", length(genes_to_exclude))

keep_features2=rownames(tPD62341.1_filtered)[which(!rownames(tPD62341.1_filtered) %in% genes_to_exclude)]
tPD62341.1_filtered2 = subset(tPD62341.1_filtered, features=keep_features2)
tPD62341.1_filtered2 = NormalizeData(object = tPD62341.1_filtered2, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
tPD62341.1_filtered2 = FindVariableFeatures(tPD62341.1_filtered2, selection.method = "vst", nfeatures = 2000) # default 2k 
tPD62341.1_filtered2 = ScaleData(tPD62341.1_filtered2, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
tPD62341.1_filtered2 = RunPCA(tPD62341.1_filtered2, npcs = 70, features=VariableFeatures(object = tPD62341.1_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
tPD62341.1_filtered2 = FindNeighbors(tPD62341.1_filtered2, dims = 1:20, verbose = T)
tPD62341.1_filtered2 = FindClusters(tPD62341.1_filtered2, resolution = 0.4, verbose = T)
tPD62341.1_filtered2 = RunTSNE(tPD62341.1_filtered2) # might as well run it but tSNE doesn't seem to be popular 

# Plot outcomes of clustering (it doesn't look massively different to me!)
pdf('Results/20241127_umap_PD62341tumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "umap", label = TRUE)
dev.off() 
pdf('Results/20241127_tsne_PD62341tumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "tsne", label = TRUE)
dev.off()
pdf('Results/20241127_pca_PD62341tumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "pca", label = TRUE)
dev.off()

tPD62341.1_filtered2 = CellCycleScoring(tPD62341.1_filtered2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
tPD62341.1_filtered2 = RunUMAP(tPD62341.1_filtered2, dims = 1:20, verbose = T, min.dist = 0.3, n.neighbors = 30)

pdf('Results/20241127_umap_PD62341tumour_rmGenes_cellCycle.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered2, reduction = "umap", label = TRUE, group.by = 'Phase')
dev.off() 

###################################################################################################################################
# BASIC QC, PREPROCESSING, CLUSTERING - FOR EACH SAMPLE SEPARATELY 

# Create a function that takes the Seurat Object and carries out basic pre-processing, clustering etc.
sc_preprocessing = function(data, sample_name, mt_threshold, feature_threshold, count_threshold, n_pcs, resolution, min_dist, nn){
  
  # Required arguments
  # data = obtained with Read10X function
  # sample_name = name of the sample to be used in plots etc.
  # mt_threshold = % of mt reads above which a cell is excluded
  
  # Initialize Seurat objects with the raw (non-normalized data)
  # Optional settings when reading in the file: 
  # min.cells = Include features detected in at least this many cells
  # min.features = Include cells where at least this many features are detected
  # read the data into the Seurat object 
  seurat = CreateSeuratObject(counts = data, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
  
  # calculate the percentage mt content 
  seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # plot the distribution of the mt content in each sample 
  # 10 sounds like a reasonable threshold and should give you most of the cells of interest 
  pdf(glue('Results/20241201_1_hist_mt_content_{sample_name}.pdf'), height = 4.5, width = 4.5)
  hist(seurat[["percent.mt"]] %>% unlist(), xlab = '% mt content', main = glue('{sample_name}'), breaks = 100)
  abline(v = 10, col = 'red', lwd = 2.5)   
  dev.off()
  
  # add CellID to the metadata 
  seurat@meta.data$CellID = rownames(seurat@meta.data) # add the cellID to the metadata
  meta = seurat@meta.data # create a dt with all the metadata of interest, including cell ID 
  
  # plot histograms of nCount_RNA_min and nFeature_RNA_min to decide the thresholds
  # usually (from a few papers I looked at), cells with fewer than 300 genes or fewer than 1000 UMIs are removed 
  
  # nCountRNA is the total number of molecules detected in each cell 
  pdf(glue('Results/20241201_1_hist_nCountRNA_{sample_name}.pdf'), height = 4, width = 8)
  ncounts_rna = as.numeric(seurat[["nCount_RNA"]] %>% unlist())
  par(mfrow=c(1,2))
  hist(ncounts_rna[ncounts_rna < 2500], xlab = 'nCountRNA', main = sample_name, breaks = 100)
  hist(ncounts_rna, xlab = 'nCountRNA', main = sample_name, breaks = 100)
  dev.off() # 1000 sounds reasonable based on the hist 
  
  # nFeatureRNA is the number of genes detected in each cell 
  pdf(glue('Results/20241201_1_hist_nFeatureRNA_{sample_name}.pdf'), height = 4, width = 8)
  nfeatures_rna = as.numeric(seurat[["nFeature_RNA"]] %>% unlist())
  par(mfrow=c(1,2))
  hist(nfeatures_rna[nfeatures_rna < 2000], xlab = 'nFeatureRNA', main = sample_name, breaks = 100)
  hist(nfeatures_rna, xlab = 'nFeatureRNA', main = sample_name, breaks = 100)
  dev.off() # 300 sounds reasonable based on the hist 

  # we can also plot those 3 QC parameters and save them to results files
  pdf(glue('Results/20241201_2_vlnPlotQC_{sample_name}.pdf'), height = 4, width = 8)
  print(VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off() # I don't see how those plots would help anyone choose the lower threshold for counts and features 
  
  # plot QC features against one another 
  pdf(glue('Results/20241201_3_featureScatter_{sample_name}.pdf'), height = 4, width = 8)
  plot1 = FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 = FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  dev.off()
  
  # which cells to keep?
  keeprows = meta[meta$nCount_RNA > count_threshold,] # require number of UMIs to be higher than minimum (1000 UMIs)
  keeprows = keeprows[keeprows$nFeature_RNA > feature_threshold,] # require number of features to be higher than minimum (300 genes)
  keeprows = keeprows[keeprows$percent.mt <= mt_threshold,] # exclude cells with mitochondrial content above a threshold (10%)
  cellids_keep = rownames(keeprows) # which CellIDs to keep 
  seurat_filtered = subset(seurat, subset = CellID %in% cellids_keep) # select cells of interest
  
  # which genes to keep?
  keep_features=rownames(seurat_filtered) 
  seurat_filtered = subset(seurat_filtered, features=keep_features)
  # I could exclude heat shock or ribosomal genes - variable takes in literature as to whether one should do this or not?
  # I won't do this in this dataset but I will run separate analysis excluding those genes
  
  # run standard Seurat pipeline 
  seurat_filtered = NormalizeData(object = seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
  seurat_filtered = FindVariableFeatures(seurat_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
  seurat_filtered = ScaleData(seurat_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
  seurat_filtered = RunPCA(seurat_filtered, npcs = n_pcs, features=VariableFeatures(object =seurat_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
  
  # Determine how many PCs to include in downstream analysis 
  pdf(glue('Results/20241201_4_elbowPlot_{sample_name}.pdf'), height = 4, width = 8)
  print(ElbowPlot(seurat_filtered, ndims = 50)) # use 50 PCs 
  dev.off()
  
  # plot variance explained (scree plot)
  mat = seurat_filtered[["RNA"]]$scale.data 
  pca = seurat_filtered[["pca"]]
  total_variance = sum(matrixStats::rowVars(mat))
  eigValues = (pca@stdev)^2 
  varExplained = eigValues / total_variance
  varExplained %>% enframe(name = "PC", value = "varExplained" ) %>%
    ggplot(aes(x = PC, y = varExplained)) + 
    geom_bar(stat = "identity") +
    theme_classic() +
    ggtitle(glue("scree plot {sample_name}"))
  ggsave(glue('Results/20241201_4_varianceExplained_{sample_name}.pdf'), height = 4, width = 4)
  
  # Clustering 
  seurat_filtered = FindNeighbors(seurat_filtered, dims = 1:20, verbose = T)
  seurat_filtered = FindClusters(seurat_filtered, resolution = resolution, verbose = T)
  seurat_filtered = RunUMAP(seurat_filtered, dims = 1:20, verbose = T, min.dist = min_dist, n.neighbors = nn)
  
  # Plot outcomes of clustering 
  pdf(glue('Results/20241201_5_umap_{sample_name}.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "umap", label = TRUE))
  dev.off() 
  
  pdf(glue('Results/20241201_5_pca_{sample_name}.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "pca", label = TRUE))
  dev.off()
  
  # Identify features and look at markers 
  seurat_filtered.markers = FindAllMarkers(seurat_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Show distirbution of cluster markers
  top4 = seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
  pdf(glue('Results/20241127_6_dotPlotTop4_{sample_name}.pdf'), width=10, height=6)
  print(DotPlot(tPD62341.1_filtered, features = top4$gene) + coord_flip())
  dev.off()
  
  # include genes that take part in the fusion, early embryonic, CD markers, synaptophysin, stuff Nathan looked at (to see if I am getting sth vaguely similar)
  genes_of_interest = c('MYOD1', 'PAX3', 'PAX7', 'CYP11A1', 'CD14', 'CD68', 'SYP', 'MN1', 'ZNF341', 'NR5A1', 'KCNQ1', 'TP53') # MYOG is 0 for all so I got rid of this
  pdf(glue('Results/20241201_6_featurePlot_{sample_name}.pdf'), width=12, height=8)
  print(FeaturePlot(seurat_filtered, features = genes_of_interest, order=TRUE))
  dev.off()
  
  # Heatmap of markers 
  top10 = seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  pdf(glue('Results/20241201_7_heatmapTop10_{sample_name}.pdf'), width=12, height=8)
  print(DoHeatmap(seurat_filtered, features = top10$gene) + NoLegend())
  dev.off()
  
  # Dotplot of genes of interest 
  pdf(glue('Results/20241201_8_dotPlot_{sample_name}_GOI.pdf'), width=10, height=6)
  print(DotPlot(object = seurat_filtered, features = genes_of_interest))
  dev.off()
  
  # Dotplot of cluster markers
  
  # Score by cell cycle expression (S and G2M genes accessed according to https://satijalab.org/seurat/articles/cell_cycle_vignette.html from October 2023)
  s.genes = cc.genes$s.genes
  g2m.genes = cc.genes$g2m.genes
  seurat_filtered = CellCycleScoring(seurat_filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  # set.ident = TRUE sets the identity of the Seurat object to the cell-cycle phase
  
  # Color UMAP by cell cycle phase (based on expression of cell cycle genes)
  pdf(glue('Results/20241201_9_umap_{sample_name}_cellCycle.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "umap", label = TRUE, group.by = 'Phase'))
  dev.off() # don't see a separate clusters of proliferating cells :(
  pdf(glue('Results/20241201_9_pca_{sample_name}_cellCycle.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "pca", label = TRUE, group.by = 'Phase'))
  dev.off() # kind of separate but not very clear? not sure if I should do anything about this 
  # I guess that means that we don't need to regress out the variation related to cell cycle phase 
  
  # Run analysis and clustering on dataset with removed ribosomal and heatshock genes
  keep_features2=rownames(seurat_filtered)[which(!rownames(seurat_filtered) %in% genes_to_exclude)]
  
  seurat_filtered2 = subset(seurat_filtered, features=keep_features2)
  seurat_filtered2 = NormalizeData(object = seurat_filtered2, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
  seurat_filtered2 = FindVariableFeatures(seurat_filtered2, selection.method = "vst", nfeatures = 2000) # default 2k 
  seurat_filtered2 = ScaleData(seurat_filtered2, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
  seurat_filtered2 = RunPCA(seurat_filtered2, npcs = 70, features=VariableFeatures(object = seurat_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
  seurat_filtered2 = FindNeighbors(seurat_filtered2, dims = 1:20, verbose = T)
  seurat_filtered2 = FindClusters(seurat_filtered2, resolution = resolution, verbose = T)
  seurat_filtered2 = RunUMAP(seurat_filtered2, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)
 
  # Plot outcomes of clustering   
  pdf(glue('Results/20241201_10_umap_{sample_name}_rmGenes.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered2, reduction = "umap", label = TRUE))
  dev.off()
  pdf(glue('Results/20241201_10_pca_{sample_name}_rmGenes.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered2, reduction = "pca", label = TRUE))
  dev.off()
  
  # Add cell cycle scoring
  seurat_filtered2 = CellCycleScoring(seurat_filtered2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  pdf(glue('Results/20241201_11_umap_{sample_name}_rmGenes_cellCycle.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered2, reduction = "umap", label = TRUE, group.by = 'Phase'))
  dev.off() 
  
}

# run for all your samples
sc_preprocessing(tumourPD62341.1, 'tPD62341.1', 10, 300, 1000, 50, 0.4, 0.3, 30) # recommended settings - I can tune the parameters but maybe get started with this 
sc_preprocessing(tumourPD62341.2, 'tPD62341.2', 10, 300, 1000, 50, 0.4, 0.3, 30)
sc_preprocessing(tumourPD62341.3, 'tPD62341.3', 10, 300, 1000, 50, 0.4, 0.3, 30)
sc_preprocessing(tumourPD63383.1, 'tPD63383.1', 10, 300, 1000, 50, 0.4, 0.3, 30)
sc_preprocessing(tumourPD63383.2, 'tPD63383.2', 10, 300, 1000, 50, 0.4, 0.3, 30)
sc_preprocessing(tumourPD63383.3, 'tPD63383.3', 10, 300, 1000, 50, 0.4, 0.3, 30)

###################################################################################################################################
# NUCLEAR RNA-seq 

# manual analysis on one nuclear sample to figure out thresholds and see how this works 

# read the data into the Seurat object 
tPD63383.nucl.1 = CreateSeuratObject(counts = tumourPD63383.1n, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)

# calculate the percentage mt content 
tPD63383.nucl.1[["percent.mt"]] = PercentageFeatureSet(tPD63383.nucl.1, pattern = "^MT-")

# plot the distribution of the mt content in each sample 
# there is a LOT of mitochondria in some sample actually ???? 
# not sure what the procedure is for dealing with those 
pdf('Results/20241202_hist_mt_content_PD63383tumour_nuclear.pdf', height = 4.5, width = 4.5)
hist(tPD63383.nucl.1[["percent.mt"]] %>% unlist(), xlab = '% mt content', main = 'PD63383n tumour', breaks = 100)
abline(v = 10, col = 'red', lwd = 2.5)   
dev.off()

# add CellID to the metadata 
tPD63383.nucl.1@meta.data$CellID = rownames(tPD63383.nucl.1@meta.data) # add the cellID to the metadata
meta = tPD63383.nucl.1@meta.data # create a dt with all the metadata of interest, including cell ID 

# plot histograms of nCount_RNA_min and nFeature_RNA_min to decide the thresholds
# usually (from a few papers I looked at), cells with fewer than 300 genes or fewer than 1000 UMIs are removed 

# nCountRNA is the total number of molecules detected in each cell 
pdf('Results/20241202_hist_nCountRNA_PD63383tumour_nuclear.pdf', height = 4, width = 8)
ncounts_rna = as.numeric(tPD63383.nucl.1[["nCount_RNA"]] %>% unlist())
par(mfrow=c(1,2))
hist(ncounts_rna[ncounts_rna < 2500], xlab = 'nCountRNA', main = 'PD63383n tumour', breaks = 100)
hist(ncounts_rna, xlab = 'nCountRNA', main = 'PD63383n tumour', breaks = 100)
dev.off() # 1000 sounds reasonable based on the hist 

# nFeatureRNA is the number of genes detected in each cell 
pdf('Results/20241202_hist_nFeatureRNA_PD63383tumour_nuclear.pdf', height = 4, width = 8)
nfeatures_rna = as.numeric(tPD63383.nucl.1[["nFeature_RNA"]] %>% unlist())
par(mfrow=c(1,2))
hist(nfeatures_rna[nfeatures_rna < 2000], xlab = 'nFeatureRNA', main = 'PD63383n tumour', breaks = 100)
hist(nfeatures_rna, xlab = 'nFeatureRNA', main = 'PD63383n tumour', breaks = 100)
dev.off() # > 100 and < 1000 sounds reasonable based on the hist 

# we can also plot those 3 QC parameters and save them to results files
pdf('Results/20241202_vlnPlotQC_PD63383tumour_nuclear.pdf', height = 4, width = 8)
par(mfrow=c(1,1))
VlnPlot(tPD63383.nucl.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off() # I don't see how those plots would help anyone choose the lower threshold for counts and features 
# very odd mt count - does this mean sth is wrong with this sample? should I keep it?

# which cells to keep?
keeprows = meta[meta$nCount_RNA > 50,] # require number of UMIs to be higher than minimum (1000 UMIs)
keeprows = keeprows[keeprows$nFeature_RNA > 10,] # require number of features to be higher than minimum (300 genes)
# keeprows = keeprows[keeprows$percent.mt <= 10,] # ignore mt filtering as this just looks odd
cellids_keep = rownames(keeprows) # which CellIDs to keep 
tPD63383.nucl.1_filtered = subset(tPD63383.nucl.1, subset = CellID %in% cellids_keep) # select cells of interest

paste('Number of cells left:', length(cellids_keep)) # 194 cells 
# should I not filter anything out then????

# which genes to keep?
keep_features=rownames(tPD63383.nucl.1_filtered) 
tPD63383.nucl.1_filtered = subset(tPD63383.nucl.1_filtered, features=keep_features)
# I could exclude heat shock or ribosomal genes - variable takes in literature as to whether one should do this or not?
# I won't do this in this dataset but I will run separate analysis excluding those genes

# run standard Seurat pipeline 
tPD63383.nucl.1_filtered = NormalizeData(object = tPD63383.nucl.1_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
tPD63383.nucl.1_filtered = FindVariableFeatures(tPD63383.nucl.1_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
tPD63383.nucl.1_filtered = ScaleData(tPD63383.nucl.1_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
tPD63383.nucl.1_filtered = RunPCA(tPD63383.nucl.1_filtered, npcs = 70, features=VariableFeatures(object = tPD63383.nucl.1_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)

# Determine how many PCs to include in downstream analysis 
pdf('Results/20241127_elbowPlot_PD63383ntumour.pdf', height = 4, width = 8)
ElbowPlot(tPD63383.nucl.1_filtered, ndims = 50) # use 50 PCs 
dev.off()

# plot variance explained (scree plot)
mat = tPD63383.nucl.1_filtered[["RNA"]]$scale.data 
pca = tPD63383.nucl.1_filtered[["pca"]]
total_variance = sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2 
varExplained = eigValues / total_variance
varExplained %>% enframe(name = "PC", value = "varExplained" ) %>%
  ggplot(aes(x = PC, y = varExplained)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  ggtitle("scree plot, PD63383n tumour") 
ggsave(glue('Results/20241127_varianceExplained_tumourPD63383n_1.pdf'), height = 4, width = 4)

# Clustering 
tPD63383.nucl.1_filtered = FindNeighbors(tPD63383.nucl.1_filtered, dims = 1:20, verbose = T)
tPD63383.nucl.1_filtered = FindClusters(tPD63383.nucl.1_filtered, resolution = 0.4, verbose = T)
tPD63383.nucl.1_filtered = RunUMAP(tPD63383.nucl.1_filtered, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)

# Plot outcomes of clustering 
pdf('Results/20241127_umap_PD63383ntumour.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered, reduction = "umap", label = TRUE)
dev.off() 
pdf('Results/20241127_pca_PD63383ntumour.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered, reduction = "pca", label = TRUE)
dev.off()

# Can I try clustering with different resolution parameter
for (res in seq(0.1, 1.2, by = 0.1)){
  tPD63383.nucl.1_filtered = FindClusters(tPD63383.nucl.1_filtered, resolution = res, verbose = T)
  tPD63383.nucl.1_filtered = RunUMAP(tPD63383.nucl.1_filtered, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)
  pdf(glue('Results/20241127_umap_PD63383ntumour_{res}_res.pdf'), height = 4, width = 4)
  plot = DimPlot(tPD63383.nucl.1_filtered, reduction = "umap", label = TRUE)
  print(plot)
  dev.off() 
} # After this I'm quite inclined to go with resolution 0.2-0.4 

# Identify features and look at markers 
tPD63383.nucl.1_filtered.markers = FindAllMarkers(tPD63383.nucl.1_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# include genes that take part in the fusion, early embryonic, CD markers, synaptophysin, stuff Nathan looked at (to see if I am getting sth vaguely similar)
genes_of_interest = c('MYOD1', 'PAX3', 'PAX7', 'CYP11A1', 'CD14', 'CD68', 'SYP', 'MN1', 'ZNF341', 'NR5A1', 'KCNQ1', 'TP53') # MYOG is 0 for all so I got rid of this
pdf(glue('Results/20241127_featurePlot_PD63383ntumour.pdf'), width=12, height=8)
FeaturePlot(tPD63383.nucl.1_filtered, features = genes_of_interest, order=TRUE)
dev.off()

# Heatmap of markers 
top10 = tPD63383.nucl.1_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(glue('Results/20241127_heatmapTop10_PD63383ntumour.pdf'), width=12, height=10)
DoHeatmap(tPD63383.nucl.1_filtered, features = top10$gene) + NoLegend()
dev.off()

# Dotplot of cluster markers (4 for each cluster)
top4 = tPD63383.nucl.1_filtered.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
pdf(glue('Results/20241127_dotPlotTop4_PD63383ntumour.pdf'), width=10, height=6)
DotPlot(tPD63383.nucl.1_filtered, features = top4$gene) + coord_flip()
dev.off()

# Dotplot of genes of interest 
pdf(glue('Results/20241127_dotPlot_PD63383ntumour_GOI.pdf'), width=10, height=6)
DotPlot(object = tPD63383.nucl.1_filtered, features = genes_of_interest)
dev.off()

# Score by cell cycle expression (S and G2M genes accessed according to https://satijalab.org/seurat/articles/cell_cycle_vignette.html from October 2023)
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes
tPD63383.nucl.1_filtered = CellCycleScoring(tPD63383.nucl.1_filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# set.ident = TRUE sets the identity of the Seurat object to the cell-cycle phase

# Color UMAP by cell cycle phase (based on expression of cell cycle genes)
pdf('Results/20241127_umap_PD63383ntumour_cellCycle.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered, reduction = "umap", label = TRUE, group.by = 'Phase')
dev.off() # don't see a separate clusters of proliferating cells :(
pdf('Results/20241127_pca_PD63383ntumour_cellCycle.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered, reduction = "pca", label = TRUE, group.by = 'Phase')
dev.off() # kind of separate but not very clear? not sure if I should do anything about this 
# I guess that means that we don't need to regress out the variation related to cell cycle phase 

# Run analysis and clustering on dataset with removed ribosomal and heatshock genes
genes_to_exclude = read.csv('Data/excludeGenes.tsv', sep='\t', header = F) %>% unlist()
excludeGenes='^(EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+_',rownames(rtoc),value=TRUE), # poorly characterized (what does this mean??)
                            grep('MALAT1',rownames(rtoc),value=TRUE), # contamination - why are we excluding this? how do we know this is contamination
                            grep('^HB[BGMQDAZE][12]?_',rownames(rtoc),value=TRUE), # contamination (apparently you should filter out blood markers)
                            grep('^MT-',rownames(rtoc),value=TRUE), # mt genes  
                            grep(excludeGenes,rownames(rtoc),value=TRUE) # housekeeping genes
))
genes_to_exclude = c(genes_to_exclud, coreExcludeGenes)
keep_features2=rownames(tPD63383.nucl.1_filtered)[which(!rownames(tPD63383.nucl.1_filtered) %in% genes_to_exclude)]

tPD63383.nucl.1_filtered2 = subset(tPD63383.nucl.1_filtered, features=keep_features2)
tPD63383.nucl.1_filtered2 = NormalizeData(object = tPD63383.nucl.1_filtered2, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
tPD63383.nucl.1_filtered2 = FindVariableFeatures(tPD63383.nucl.1_filtered2, selection.method = "vst", nfeatures = 2000) # default 2k 
tPD63383.nucl.1_filtered2 = ScaleData(tPD63383.nucl.1_filtered2, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
tPD63383.nucl.1_filtered2 = RunPCA(tPD63383.nucl.1_filtered2, npcs = 70, features=VariableFeatures(object = tPD63383.nucl.1_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
tPD63383.nucl.1_filtered2 = FindNeighbors(tPD63383.nucl.1_filtered2, dims = 1:20, verbose = T)
tPD63383.nucl.1_filtered2 = FindClusters(tPD63383.nucl.1_filtered2, resolution = 0.6, verbose = T)
tPD63383.nucl.1_filtered2 = RunUMAP(tPD63383.nucl.1_filtered2, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)
tPD63383.nucl.1_filtered2 = RunTSNE(tPD63383.nucl.1_filtered2) # might as well run it but tSNE doesn't seem to be popular 
tPD63383.nucl.1_filtered2 = CellCycleScoring(tPD63383.nucl.1_filtered2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Plot outcomes of clustering (it doesn't look massively different to me!)
pdf('Results/20241127_umap_PD63383ntumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered2, reduction = "umap", label = TRUE)
dev.off() 
pdf('Results/20241127_tsne_PD63383ntumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered2, reduction = "tsne", label = TRUE)
dev.off()
pdf('Results/20241127_pca_PD63383ntumour_rmGenes.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered2, reduction = "pca", label = TRUE)
dev.off()
pdf('Results/20241127_umap_PD63383ntumour_rmGenes_cellCycle.pdf', height = 4, width = 4)
DimPlot(tPD63383.nucl.1_filtered2, reduction = "umap", label = TRUE, group.by = 'Phase')
dev.off() 

# try running the same function on nuclear samples, but note that these may require different QC
# does filtering based on mt still make sense?
# should I lower the thresholds for feature and gene count (since data is more sparse)?
sc_preprocessing(tumourPD63383.1n, 'tPD63383.nuclear.1', 10, 200, 500, 50, 0.2, 0.5, 30)
sc_preprocessing(tumourPD63383.2n, 'tPD63383.nuclear.2', 10, 200, 500, 50, 0.2, 0.5, 30)
sc_preprocessing(tumourPD63383.3n, 'tPD63383.nuclear.3', 10, 200, 500, 50, 0.2, 0.5, 30)

###################################################################################################################################
# INTEGRATION OF SAMPLES FROM DIFFERENT LANES

# at what point should I integrate samples from different lanes?
# I think the most common way to do this is to integrate the .fastq reads at the stage of running CellRanger
# However, I also found this https://ucdavis-bioinformatics-training.github.io/2020-August-intro-scRNAseq/data_analysis/scRNA_Workshop-PART1_fixed#:~:text=Load%20the%20Cell%20Ranger%20Matrix,can%20aggregate%20them%20in%20R.
# which seems to suggest that you can integrate samples from multiple lanes when running Seurat (after CellRanger, working on CellRanger count matrix)

setwd('/Users/bw18/Desktop/1SB/scRNAseq')

# Samples 
# I have the following samples (these are tumour samples!)
# PD62341 fresh: CG_SB_NB13652544, CG_SB_NB13652545, CG_SB_NB13652546
# PD63383 fresh: CG_SB_NB13760628, CG_SB_NB13760629, CG_SB_NB13760630
# PD63383 frozen (nuclear): CG_SB_NB14599068, CG_SB_NB14599069, CG_SB_NB14599070

PD62341_ids = c("CG_SB_NB13652544", "CG_SB_NB13652545", "CG_SB_NB13652546")
PD63383_ids = c("CG_SB_NB13760628", "CG_SB_NB13760629", "CG_SB_NB13760630")

data.loc = "/Users/bw18/Desktop/1SB/scRNAseq"

PD62341.d = sapply(PD62341_ids, function(i){
  d10x = Read10X(file.path(data.loc, "Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
PD62341.data = do.call("cbind", PD62341.d)

PD63383.d = sapply(PD63383_ids, function(i){
  d10x = Read10X(file.path(data.loc, "Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
PD63383.data = do.call("cbind", PD63383.d)

# We can now run our function on the aggregated data and see what happens 
sc_preprocessing(PD62341.data, 'tPD62341.agg', 10, 300, 1000, 50, 0.4, 0.3, 30)
sc_preprocessing(PD63383.data, 'tPD63383.agg', 10, 300, 1000, 50, 0.4, 0.3, 30)

# What is the correct way of merging the datasets from both twins together?
# Can try just this but I am not convinced this is the way to do this 
tumour_sc_ids = c("CG_SB_NB13652544", "CG_SB_NB13652545", "CG_SB_NB13652546", "CG_SB_NB13760628", "CG_SB_NB13760629", "CG_SB_NB13760630")
tumour_sc.d = sapply(tumour_sc_ids, function(i){
  d10x = Read10X(file.path(data.loc, "Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

tumour_sc.data = do.call("cbind", tumour_sc.d)
sc_preprocessing(tumour_sc.data, 'tumour_sc.agg', 10, 300, 1000, 50, 0.4, 0.3, 30)

###################################################################################################################################
# Cell type annotation 

# One possible option is to look at expression of known markers

# I found https://www-nature-com.ezp.lib.cam.ac.uk/articles/s41388-024-03001-8 (March 2024, Oncogene)
# which reports the scRNA-seq of undifferentiated pleomorphic sarcoma (which is kind of what I have; pleomorphic means occurring in different forms)
# they identify: monocytes (macrophages), tumour cells, fibroblast, T cell, endothelial, mast, NK, pericyte cells
# perhaps do this for tumour from each twin separately and then can compare clusters 

# run pre-processing on aggregated tumour from one twin
PD62341.d = sapply(PD62341_ids, function(i){
  d10x = Read10X(file.path(data.loc, "Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
PD62341.data = do.call("cbind", PD62341.d) # that just treats this as having > 20k cells - is this the right way to do this?

PD62341t = CreateSeuratObject(counts = PD62341.data, project = "PD62341.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
PD62341t[["percent.mt"]] = PercentageFeatureSet(PD62341t, pattern = "^MT-")
PD62341t@meta.data$CellID = rownames(PD62341t@meta.data) # add the cellID to the metadata
meta = PD62341t@meta.data # create a dt with all the metadata of interest, including cell ID 
keeprows = meta[meta$nCount_RNA > 1000 & meta$nCount_RNA < 30000,] # require number of UMIs to be higher than minimum (1000 UMIs)
keeprows = keeprows[keeprows$nFeature_RNA > 300 & meta$nFeature_RNA < 1000,] # require number of features to be higher than minimum (300 genes)
keeprows = keeprows[keeprows$percent.mt <= 10,] # exclude cells with mitochondrial content above a threshold (10%)
cellids_keep = rownames(keeprows) # which CellIDs to keep 
PD62341t_filtered = subset(PD62341t, subset = CellID %in% cellids_keep) # select cells of interest
keepfeatures=rownames(PD62341t_filtered)[which(!rownames(PD62341t_filtered) %in% genes_to_exclude)]
PD62341t_filtered = subset(PD62341t_filtered, features=keepfeatures)
PD62341t_filtered = NormalizeData(object = PD62341t_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
PD62341t_filtered = FindVariableFeatures(PD62341t_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
PD62341t_filtered = ScaleData(PD62341t_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
PD62341t_filtered = RunPCA(PD62341t_filtered, npcs = 50, features=VariableFeatures(object = PD62341t_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
PD62341t_filtered = FindNeighbors(PD62341t_filtered, dims = 1:50, verbose = T)
PD62341t_filtered = FindClusters(PD62341t_filtered, resolution = 0.4, verbose = T)
PD62341t_filtered = RunUMAP(PD62341t_filtered, dims = 1:50, verbose = T, min.dist = 0.6, n.neighbors = 40)

pdf(glue('Results/20241201_umap_PD62341t_rmGenes_aggTumour.pdf'), height = 4, width = 4)
DimPlot(PD62341t_filtered, reduction = "umap", label = TRUE)
dev.off() 

# first, start by looking at expression of each of the possible marker genes 
marker_genes = c("PTPRC", "CD68", "CD14", "CD163", "AIF1", # monocyte 
                 "COL11A1", "COL3A1", # fibroblasts 
                 "CD8A", "CD4", "CD3D", # T cells
                 "GNLY", "NKG7", # NK cells
                 "PLVAP", "PECAM1", "VWF", # endothelial
                 "ACTA2", "NOTCH3", "RGS5", # pericyte 
                 "ZNF341", "MN1", # tumour cells 
                 "TMEM119", "C1Q", "PROS1", "IBA1")  

pdf(glue('Results/20241201_featurePlot_marker_genes_tPD62341_agg.pdf'), width=12, height=16)
FeaturePlot(PD62341t_filtered, features = marker_genes, order=TRUE)
dev.off()

###################################################################################################################################
# Cell type annotation (SEURAT TUTORIAL) - I TRIED BUT DOESN'T SEEM TO WORK WITH THE FETAL DATASET

# another approach is to maybe try and do data integration with Seurat 
# library(devtools)
# devtools::install_github('satijalab/seurat-data')
library(SeuratData)
AvailableData() # what kinds of datasets would be useful for this sample? some sarcoma? any fetal dataset?

# install a reference package and load the data from the reference 
InstallData("panc8")
panc8 = LoadData("panc8")
pancreas.ref = subset(panc8, tech %in% c("celseq2", "smartseq2")) # select scRNA-seq generated with Cel-seq and SMART-seq
pancreas.ref[["RNA"]] = split(pancreas.ref[["RNA"]], f = pancreas.ref$tech)

# data pre-processing (no integration)
pancreas.ref = NormalizeData(pancreas.ref)
pancreas.ref = FindVariableFeatures(pancreas.ref)
pancreas.ref = ScaleData(pancreas.ref)
pancreas.ref = RunPCA(pancreas.ref)
pancreas.ref = FindNeighbors(pancreas.ref, dims = 1:30)
pancreas.ref = FindClusters(pancreas.ref)
pancreas.ref = RunUMAP(pancreas.ref, dims = 1:30)
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))

# integrate datasets from different sequencing tech into a shared reference
pancreas.ref = IntegrateLayers(object = pancreas.ref, method = CCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.cca", verbose = FALSE)
pancreas.ref = FindNeighbors(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
pancreas.ref = FindClusters(pancreas.ref)
pancreas.ref = RunUMAP(pancreas.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(pancreas.ref, group.by = c("tech", "celltype"))

# create a query dataset (two other technologies)
pancreas.query = subset(panc8, tech %in% c("fluidigmc1", "celseq"))
pancreas.query = NormalizeData(pancreas.query)
pancreas.anchors = FindTransferAnchors(reference = pancreas.ref, query = pancreas.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions = TransferData(anchorset = pancreas.anchors, refdata = pancreas.ref$celltype, dims = 1:30)
pancreas.query = AddMetaData(pancreas.query, metadata = predictions)
pancreas.query$prediction.match = pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)

# use reference dataset to assign cell type to the query
pancreas.ref = RunUMAP(pancreas.ref, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
pancreas.query = MapQuery(anchorset = pancreas.anchors, reference = pancreas.ref, query = pancreas.query,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
p1 = DimPlot(pancreas.ref, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 = DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

# maybe we could try with fetal data as reference? could identify immune cells / fibroblasts / endothelial
remotes::install_github('satijalab/azimuth', ref = 'master')
InstallData("fetusref")
fetal = LoadData("fetusref", "azimuth") 
# I am getting an error where fetusref is installed but cannot be loaded 
# attempted solving by installing Azimuth package (also Satija lab) but doesn't seem like I can get all the packages reasonably easily so I think I'm going to drop this for now
fetal.ref = subset(fetal) # select scRNA-seq generated with Cel-seq and SMART-seq
fetal.ref[["RNA"]] = split(fetal.ref[["RNA"]], f = fetal.ref$tech)

# data pre-processing (no integration)
fetal.ref = NormalizeData(fetal.ref)
fetal.ref = FindVariableFeatures(fetal.ref)
fetal.ref = ScaleData(fetal.ref)
fetal.ref = RunPCA(fetal.ref)
fetal.ref = FindNeighbors(fetal.ref, dims = 1:30)
fetal.ref = FindClusters(fetal.ref)
fetal.ref = RunUMAP(fetal.ref, dims = 1:30)
DimPlot(fetal.ref, group.by = c("tech", "celltype"))

# integrate datasets from different sequencing tech into a shared reference
fetal.ref = IntegrateLayers(object = fetal.ref, method = CCAIntegration, orig.reduction = "pca",
                               new.reduction = "integrated.cca", verbose = FALSE)
fetal.ref = FindNeighbors(fetal.ref, reduction = "integrated.cca", dims = 1:30)
fetal.ref = FindClusters(fetal.ref)
fetal.ref = RunUMAP(fetal.ref, reduction = "integrated.cca", dims = 1:30)
DimPlot(fetal.ref, group.by = c("tech", "celltype"))


# Cell type annotation based on markers of each cluster 

# Let's agree that we are using clusters with the following parameters:
# N_PCs = 70, resolution = 0.2, min_dist = 0.3, n_neighbors = 30
PD62341t = CreateSeuratObject(counts = PD62341.data, project = "PD62341.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
PD62341t[["percent.mt"]] = PercentageFeatureSet(PD62341t, pattern = "^MT-")
PD62341t@meta.data$CellID = rownames(PD62341t@meta.data) # add the cellID to the metadata
meta = PD62341t@meta.data # create a dt with all the metadata of interest, including cell ID 
keeprows = meta[meta$nCount_RNA > 1000,] # require number of UMIs to be higher than minimum (1000 UMIs)
keeprows = keeprows[keeprows$nFeature_RNA > 300,] # require number of features to be higher than minimum (300 genes)
keeprows = keeprows[keeprows$percent.mt <= 10,] # exclude cells with mitochondrial content above a threshold (10%)
cellids_keep = rownames(keeprows) # which CellIDs to keep 
PD62341t_filtered = subset(PD62341t, subset = CellID %in% cellids_keep) # select cells of interest
keepfeatures=rownames(PD62341t_filtered)[which(!rownames(PD62341t_filtered) %in% genes_to_exclude)]
PD62341t_filtered = subset(PD62341t_filtered, features=keepfeatures)
PD62341t_filtered = NormalizeData(object = PD62341t_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
PD62341t_filtered = FindVariableFeatures(PD62341t_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
PD62341t_filtered = ScaleData(PD62341t_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
PD62341t_filtered = RunPCA(PD62341t_filtered, npcs = 70, features=VariableFeatures(object = PD62341t_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
PD62341t_filtered = FindNeighbors(PD62341t_filtered, dims = 1:70, verbose = T)
PD62341t_filtered = FindClusters(PD62341t_filtered, resolution = 0.2, verbose = T)
PD62341t_filtered = RunUMAP(PD62341t_filtered, dims = 1:70, verbose = T, min.dist = 0.3, n.neighbors = 30)
DimPlot(PD62341t_filtered, reduction = "umap", label = TRUE) # okay so we have 14 clusters for the moment 

pdf('Results/20241201_umap_PD62341t_14clusters.pdf')
DimPlot(PD62341t_filtered, reduction = "umap", label = TRUE) # okay so we have 14 clusters for the moment 
dev.off()

# identify top 5 genes for each cluster 
PD62341t_filtered.markers = FindAllMarkers(PD62341t_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top4 = PD62341t_filtered.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
DotPlot(PD62341t_filtered, features = top4$gene) + coord_flip()

top10 = PD62341t_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf('Results/20241201_dotplot_PD62341t_14clusters.pdf', height=16, width=10)
DotPlot(PD62341t_filtered, features = top4$gene) + coord_flip()
dev.off()

# also do a feature plot as this will help identify the easiest clusters (tumour, fibroblasts)
pdf(glue('Results/20241201_featurePlot_knownMarkers_tPD62341_agg_14clusters.pdf'), width=12, height=16)
FeaturePlot(PD62341t_filtered, features = marker_genes, order=TRUE)
dev.off()

# checking genes for each cluster 
# 0 - looks like something from the bone marrow. could be macrophages? weirdly also has some brain-related stuff. but also, very low PTPRC expression?
# 1 - tumour 
# 2 - endothelial 
# 3 - looks like some cells from the bone marrow 
# 4 - tumour 
# 5 - looks like meyloid cells (+ immune subset)
# 6 - fibroblasts?
# 7 - express PTPRC and some immune-enriched stuff, so likely another myeloid set
# 8 - fibroblasts?
# 9 - tumour 
# 10 - fibroblasts?
# 11 - quite unclear 
# 12 - quite unclear 
# 13 - looks endothelial

# I am not sure what clusters 0 and 3 are 
# Find all markers for those two clusters
cluster0.markers <- FindMarkers(PD62341t_filtered, ident.1 = 0)
# based on the stuff that's up (HBA1, HBA2, HGB1, ALAS2, TMCC2) - feels like this has to be erythrocytes
cluster3.markers <- FindMarkers(PD62341t_filtered, ident.1 = 3)
# also feels like this has to be erythroid cells 

# can we use markers for blood cells?
# PBMC markers from https://satijalab.org/seurat/articles/pbmc3k_tutorial.html 
pbmc_markers = c("CCR7", # naive T cell (CD4)
                 "CD14", "LYZ", # monocyte 
                 "S100A4", # memory T cell (CD4)
                 "CD8A", # CD8+ T cell 
                 "CD14", "CD68", # myeloid-like stuff I guess
                 "FCGR3A", "MS4A7", # monocyte 
                 "GNLY", "NKG7", # NK
                 "FCER1A", "CST3", # DC
                 "CD34", "RUNDC3A", "TMEM119",
                 "CD45",
                 "HBA1", "HBA2", "HBB", "ALAS2", "AHSP") # hemoglobin 

pdf(glue('Results/20241201_featurePlot_pbcmMarkers_tPD62341_agg.pdf'), width=12, height=12)
FeaturePlot(PD62341t_filtered, features = pbmc_markers, order=TRUE)
dev.off()
 
# Another way of doing this is what I found here: https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html

###################################################################################################################################
# Logistic regression on HCA 

# the first thing I can do is probably check the expression of some adrenal markers in my dataset
# markers from https://pmc.ncbi.nlm.nih.gov/articles/PMC10443814/#sd 
adrenal_markers = c("NR5A1", "CYP11A1", "MC2R",
                    "STAR", "CLRN1", "MIR202HG", "FAM166B", "CYP11B1", "MT3", 
                    "CCN3", "SULT2A1") # CCN3 is a DZ (definitive zone) marker, SULT2A1 is an FZ (fetal zone) marker 

pdf(glue('Results/20241201_featurePlot_adrenalMarkers_tPD62341_agg.pdf'), width=12, height=8)
FeaturePlot(PD62341t_filtered, features = adrenal_markers, order=TRUE)
dev.off()

# generally I feel like I am getting too many clusters
# tumour clearly expresses some adrenal markers (like CYP11A1 and SULT2A1)
# however, the 0-3 cluster also seems to express some of those 

# Logistic regression on HCA 
# the way to do this is described in Matthew Young's paper here:
# https://www.science.org/doi/10.1126/science.aat1699

# how this is done in the paper:
# To measure the similarity of tumor cells to either normal mature kidney cells or fetal cells,
# we trained a logistic regression model using elastic net regularization on the cellular identities
# defined by the clusters in normal mature and fetal epithelial and vascular compartments. In
# training this model we set alpha=0.99 to produce strong regularization but to prevent the
# exclusion of strongly co-linear genes (would be good for you to check what they mean by regularization)

# To obtain regression coefficients specific to each cluster in our training data we fit a series
# of N binomial logistic regression models, where N is the number of clusters in the training data
# (i.e., one-versus-rest binomial logistic regression). To prevent the observed frequencies of cells
# (which we do not expect to accurately reflect the true abundances in situ) from biasing the
# regression coefficients we use an offset for each model given by log(f / (1-f))
# where f is the fraction of cells in the cluster being trained.

# In each case, we performed 10-fold cross validation and selected the regularization co-
# efficient, lambda, to be as large as possible (i.e., as few non-zero coefficients as possible) such
# that the cross validated accuracy was within 1 standard deviation of the minimum.
# These models were then used to calculate a predicted similarity to each of the fetal and
# normal clusters for each cell in the tumor map. In calculating the predicted values, an offset of 0
# was used. Softmax normalization was not used to allow for the possibility that tumor cells do
# not resemble any of the fetal or normal cells in the training set. Predicted logits were then
# averaged within each tumor cluster and converted to probabilities for visualization.

# required libraries
library(glmnet)
library(ggplot2)
library(cowplot)
library(foreach)
library(doMC)
BiocManager::install(c("monocle"))
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
devtools::install_github("cole-trapnell-lab/garnett")

# Run logistic regression analysis 
# code from https://github.com/constantAmateur/scKidneyTumors/blob/master/postPipeline.R

# define required functions (as per https://github.com/constantAmateur/scKidneyTumors/blob/master/similarity.R)
# get offset (what does this actually mean?)
getPopulationOffset = function(y){
  if(!is.factor(y)) # convert the input to factor level 
    y=factor(y)
  if(length(levels(y))!=2)
    stop("y must be a two-level factor")
  off = sum(y==levels(y)[2])/length(y)
  off = log(off/(1-off))
  return(rep(off,length(y)))
}

# Do the OvR fit for every variable.  This just does a simple CV selection of regularisation amount.  Far from ideal, but should be good enough for the main conclusions.
multinomialFitCV = function(x,y,nParallel=1,...){
  fits = list()
  if(nParallel>1)
    registerDoMC(cores=nParallel)
  #Do them in order of size
  marks = names(sort(table(as.character(y))))
  for(mark in marks){
    message(sprintf("Fitting model for variable %s",mark))
    fac = factor(y==mark)
    #The two main modes of failure are too few positives and errors constructing lambda.  These should be handled semi-gracefully
    fits[[mark]] = tryCatch(
      cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,...),
      error = function(e) {
        tryCatch(
          cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,lambda=exp(seq(-10,-3,length.out=100)),...),
          error = function(e) {
            warning(sprintf("Could not fit model for variable %s",mark))
            return(NULL)
          })
      })
  }
  return(fits)
}

# Load training data
loadTrainingData = function(dirs){
  mDat = NULL
  toc = NULL
  for(dir in dirs){
    metadata = readRDS(file.path(dir,'metadata.RDS'))
    ttoc = readRDS(file.path(dir,'tableOfCounts.RDS'))
    rownames(metadata$meta.data) = paste0(metadata$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(metadata$meta.data)))
    colnames(ttoc) = rownames(metadata$meta.data)
    if(is.null(mDat)){
      mDat = metadata$meta.data
    }else{
      mDat = rbind(mDat,metadata$meta.data)
    }
    if(is.null(toc)){
      toc = ttoc
    }else{
      toc = cbind(toc,ttoc)
    }
  }
  return(list(toc=toc,mDat=mDat))
}

# where do I get the data from??
# I am considering possibly this https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
# (you can try and get the .Rds data files for the different organs that they have)

adrenal_rds = readRDS('Data/adrenal_classifier.Rds')


# read the data 
plotDir='/Results'
tfs = read.table('~/Projects/Common/Homo_sapiens_transcription_factors_gene_list.txt',sep='\t',header=TRUE)
rtoc = readRDS('~/scratch/KidneySC/globalRawTableOfCounts.RDS')

# exclude the following genes
excludeGenes='^(EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'
coreExcludeGenes = unique(c(grep('\\.[0-9]+_',rownames(rtoc),value=TRUE), # poorly characterized (what does this mean??)
                            grep('MALAT1',rownames(rtoc),value=TRUE), # contamination - why are we excluding this? how do we know this is contamination
                            grep('^HB[BGMQDAZE][12]?_',rownames(rtoc),value=TRUE), # contamination (apparently you should filter out blood markers)
                            grep('^MT-',rownames(rtoc),value=TRUE), # mt genes  
                            grep(excludeGenes,rownames(rtoc),value=TRUE) # housekeeping genes
))


# Training models 

###################
# Foetal Clusters
trainDat = loadTrainingData(file.path(plotDir,'barcodeLevelGroup_foetalEpitheliumAndVascularV3/batchCorrected/'))
trainDat$mDat$Trainer = trainDat$mDat$res.1

#Add in the MAST cells
extras = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_tumourImmune/batchCorrected/'))
#Get the cells and their class
cells = c(rownames(trainDat$mDat),rownames(extras$mDat[extras$mDat$res.1=='16',]))
classes = c(trainDat$mDat$res.1,rep('MAST',sum(extras$mDat$res.1=='16')))
#Any genes to exclude go here
excludeGenes=coreExcludeGenes
#Any genes that we wish to give extra weight should go here
includeGenes=c()

#Get and normalise the data
dat = Seurat::LogNormalize(rtoc[,cells])
dat = dat[(rownames(dat) %in% includeGenes) | (Matrix::rowSums(dat>0)>3 & !(rownames(dat)%in%excludeGenes)),]
dat = t(dat)
fitFoetalClusters = sim$multinomialFitCV(dat,classes,nParallel=10)
#Get the genes that are going to be informative
saveRDS(fitFoetalClusters,file.path(plotDir,'lrFoetalClustersV4.RDS'))

###################
# Normal Clusters
trainDat = loadTrainingData(file.path(plotDir,'barcodeLevelGroup_normalEpitheliumAndVascularWithoutProximalTubularV2/batchCorrected/'))
#Label and merge clusters

#annotMap = c(EN17='G',EN0='PT1',EN7='PT1',EN6='PT2',EN5='PT2',EN4='PT2',EN9='PT2',EN3='PT2',EN18='PT3',EN19='D',EN21='C1',EN22='C2',EN20='P',EN15='U1/2',EN25='U3/4',EN23='M',EN24='F',EN16='GE',EN12='DV',EN10='AV')
#Drop the ambiguous clusters and merge the ones that should really be the same
#trainDat$mDat = trainDat$mDat[paste0('EN',trainDat$mDat$res.1) %in% names(annotMap),]
#trainDat$mDat$Trainer = annotMap[paste0('EN',trainDat$mDat$res.1)]
trainDat$mDat$Trainer = trainDat$mDat$res.1

#Add in the MAST cells
extras = sim$loadTrainingData(file.path(plotDir,'barcodeLevelGroup_tumourImmune/batchCorrected/'))

#Get the cells and their class
cells = c(rownames(trainDat$mDat),rownames(extras$mDat[extras$mDat$res.1=='16',]))
classes = c(trainDat$mDat$res.1,rep('MAST',sum(extras$mDat$res.1=='16')))

#Any genes to exclude go here
excludeGenes=coreExcludeGenes

#Any genes that we wish to give extra weight should go here
includeGenes=c()

#Get and normalise the data
dat = Seurat::LogNormalize(rtoc[,cells])
dat = dat[(rownames(dat) %in% includeGenes) | (Matrix::rowSums(dat>0)>3 & !(rownames(dat)%in%excludeGenes)),]
dat = t(dat)
fitNormalClusters = multinomialFitCV(dat,classes,nParallel=10)

#Get the genes that are going to be informative
saveRDS(fitNormalClusters,file.path(plotDir,'lrNormalClustersWithControl.RDS'))

###################################################################################################################################
# Myeloid cells - what's up with those? 
# where are the myeloid cells and what can I say about the expression of CD14 and CD68?

# How do I get reads from those cells to do the genotyping?




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
twins[["percent.mt"]] = PercentageFeatureSet(twins, pattern = "^MT-")
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
cluster1.markers = FindMarkers(twins, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5) 

# find markers that differentiate cluster 5 from clusters 0 and 3
cluster5.markers = FindMarkers(twins, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5) # ident.2 sets what you want to distinguish your cells from

# find markers for each cluster
twins.markers = FindAllMarkers(twins, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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
p_after$data = p_after$data[order(p_after$data$CD14),]
CombinePlots(plots = list(p, p_after))

###################################################################################################################################
# Expression heatmap for different clusters 
top10 = twins.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(twins, features = top10$gene) + NoLegend()
# omg this looks like it did actually do something that is reasonably correct 





