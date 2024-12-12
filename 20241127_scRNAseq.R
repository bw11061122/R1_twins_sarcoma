###################################################################################################################################
# SCRIPT 6

# Script to remove germline mutations in copy number-altered regions
# November - December 2024
# Barbara Walkowiak bw18

# INPUT: 
# 1 scRNAseq CellRanger output
# 2 list of genes to exclude (from Matthew Young's code)

# OUTPUT: 
# 1 basic scRNA-seq analysis (QC + clustering) 

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
library(Matrix)
library(tibble)
library(glue)

###################################################################################################################################
# INPUT DATA 

setwd('/Users/bw18/Desktop/1SB')

# Samples 
# I have the following samples (these are tumour samples!)
# PD62341 fresh: CG_SB_NB13652544, CG_SB_NB13652545, CG_SB_NB13652546
# PD63383 fresh: CG_SB_NB13760628, CG_SB_NB13760629, CG_SB_NB13760630
# PD63383 frozen (nuclear): CG_SB_NB14599068, CG_SB_NB14599069, CG_SB_NB14599070

# Read the input data 
# PD62341 tumour 
tumourPD62341.1 = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB13652544/filtered_feature_bc_matrix/") # lane 1
tumourPD62341.2 = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB13652545/filtered_feature_bc_matrix/") # lane 2
tumourPD62341.3 = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB13652546/filtered_feature_bc_matrix/") # lane 3

# PD63383 tumour 
tumourPD63383.1 = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB13760628/filtered_feature_bc_matrix/") # lane 1
tumourPD63383.2 = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB13760629/filtered_feature_bc_matrix/") # lane 2
tumourPD63383.3 = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB13760630/filtered_feature_bc_matrix/") # lane 3

# PD62341 tumour (frozen, nuclear)
tumourPD63383.1n = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB14599068/filtered_feature_bc_matrix/") # lane 1
tumourPD63383.2n = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB14599069/filtered_feature_bc_matrix/") # lane 2
tumourPD63383.3n = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB14599070/filtered_feature_bc_matrix/") # lane 3

# Genes that can be excluded from clustering (e.g., ribosomal genes)
# List of genes from Matthew Young (github scRNAseq kidney paper https://github.com/constantAmateur/scKidneyTumors)
excludeGenes='^(EIF[0-9]|RPL[0-9]|RPS[0-9]|RPN1|POLR[0-9]|SNX[0-9]|HSP[AB][0-9]|H1FX|H2AF[VXYZ]|PRKA|NDUF[ABCSV]|ATP[0-9]|PSM[ABCDEFG][0-9]|UBA[0-9]|UBE[0-9]|USP[0-9]|TXN)'

###################################################################################################################################
# BASIC QC AND PRE-PROCESSING: processing for a single lane to establish appropriate thresholds 

# manual analysis on one sample to figure out thresholds and see how this works 

# read the data into the Seurat object 
tPD62341.1 = CreateSeuratObject(counts = tumourPD62341.1, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)

# calculate the percentage mt content 
tPD62341.1[["percent.mt"]] = PercentageFeatureSet(tPD62341.1, pattern = "^MT-")

# plot the distribution of the mt content in each sample 
# 10 sounds like a reasonable threshold and should give you most of the cells of interest 
pdf('Figures/F6/20241208_hist_mt_content_PD62341tumour_sample1.pdf', height = 4.5, width = 4.5)
hist(tPD62341.1[["percent.mt"]] %>% unlist(), xlab = '% mt content', main = 'PD62341 tumour', breaks = 100)
abline(v = 10, col = 'red', lwd = 2.5)   
dev.off()

# add CellID to the metadata 
tPD62341.1@meta.data$CellID = rownames(tPD62341.1@meta.data) # add the cellID to the metadata
meta = tPD62341.1@meta.data # create a dt with all the metadata of interest, including cell ID 

# plot histograms of nCount_RNA_min and nFeature_RNA_min to decide the thresholds
# usually (from a few papers I looked at), cells with fewer than 300 genes or fewer than 1000 UMIs are removed 

# nCountRNA is the total number of molecules detected in each cell 
pdf('Figures/F6/20241208_hist_nCountRNA_PD62341tumour_sample1.pdf', height = 4, width = 8)
ncounts_rna = as.numeric(tPD62341.1[["nCount_RNA"]] %>% unlist())
par(mfrow=c(1,2))
hist(ncounts_rna[ncounts_rna < 2500], xlab = 'nCountRNA', main = 'PD62341 tumour', breaks = 100)
hist(ncounts_rna, xlab = 'nCountRNA', main = 'PD62341 tumour', breaks = 100)
dev.off() # 1000 molecules sounds reasonable based on the hist 

# nFeatureRNA is the number of genes detected in each cell 
pdf('Figures/F6/20241208_hist_nFeatureRNA_PD62341tumour_sample1.pdf', height = 4, width = 8)
nfeatures_rna = as.numeric(tPD62341.1[["nFeature_RNA"]] %>% unlist())
par(mfrow=c(1,2))
hist(nfeatures_rna[nfeatures_rna < 2000], xlab = 'nFeatureRNA', main = 'PD62341 tumour', breaks = 100)
hist(nfeatures_rna, xlab = 'nFeatureRNA', main = 'PD62341 tumour', breaks = 100)
dev.off() # 300 features reasonable based on the hist 

# we can also plot those 3 QC parameters and save them to results files
pdf('Figures/F6/20241208_vlnPlotQC_PD62341tumour_sample1.pdf', height = 4, width = 8)
par(mfrow=c(1,1))
VlnPlot(tPD62341.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off() # I don't see how those plots would help anyone choose the lower threshold for counts and features 

# plot QC features against one another 
pdf('Figures/F6/20241208_featureScatter_PD62341tumour_sample1.pdf', height = 4, width = 12)
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
genes_to_exclude = unique(c(grep('\\.[0-9]+_',rownames(tPD62341.1_filtered),value=TRUE), # poorly characterized (what does this mean??)
                            grep('MALAT1',rownames(tPD62341.1_filtered),value=TRUE), # contamination - why are we excluding this? how do we know this is contamination
                            grep('^HB[BGMQDAZE][12]?_',rownames(tPD62341.1_filtered),value=TRUE), # contamination (apparently you should filter out blood markers)
                            grep('^MT-',rownames(tPD62341.1_filtered),value=TRUE), # mt genes  
                            grep(excludeGenes,rownames(tPD62341.1_filtered),value=TRUE) # housekeeping genes
))
keep_features=setdiff(rownames(tPD62341.1_filtered), genes_to_exclude)
paste('Number of features to keep:', length(keep_features)) # 33,057
tPD62341.1_filtered = subset(tPD62341.1_filtered, features=keep_features)
# I could exclude heat shock or ribosomal genes - variable takes in literature as to whether one should do this or not?
# I won't do this in this dataset but I will run separate analysis excluding those genes

# run standard Seurat pipeline 
tPD62341.1_filtered = NormalizeData(object = tPD62341.1_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
tPD62341.1_filtered = FindVariableFeatures(tPD62341.1_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
tPD62341.1_filtered = ScaleData(tPD62341.1_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
tPD62341.1_filtered = RunPCA(tPD62341.1_filtered, npcs = 70, features=VariableFeatures(object = tPD62341.1_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)

# Determine how many PCs to include in downstream analysis 
pdf('Figures/F6/20241208_elbowPlot_PD62341tumour_sample1.pdf', height = 4, width = 8)
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
ggsave(glue('Figures/F6/20241208_varianceExplained_tumourPD62341_sample1.pdf'), height = 4, width = 4)

# Clustering 
tPD62341.1_filtered = FindNeighbors(tPD62341.1_filtered, dims = 1:20, verbose = T)
tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = 0.4, verbose = T)
tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.3, n.neighbors = 30)
tPD62341.1_filtered = RunTSNE(tPD62341.1_filtered) # might as well run it but tSNE doesn't seem to be popular 

# Plot outcomes of clustering 
pdf('Figures/F6/20241208_umap_PD62341tumour_sample1.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "umap", label = TRUE)
dev.off() 
pdf('Figures/F6/20241208_tsne_PD62341tumour_sample1.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "tsne", label = TRUE)
dev.off()
pdf('Figures/F6/20241208_pca_PD62341tumour_sample1.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "pca", label = TRUE)
dev.off()

# Can I try clustering with different resolution parameter
for (res in seq(0.1, 1.2, by = 0.1)){
  tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = res, verbose = T)
  tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.5, n.neighbors = 30)
  pdf(glue('Figures/F6/20241208_umap_PD62341tumour_sample1_res{res}.pdf'), height = 4, width = 4)
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
pdf(glue('Figures/F6/20241208_featurePlot_PD62341tumour_sample1.pdf'), width=12, height=8)
FeaturePlot(tPD62341.1_filtered, features = genes_of_interest, order=TRUE)
dev.off()

# Heatmap of markers 
top10 = tPD62341.1_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf(glue('Figures/F6/20241208_heatmapTop10_PD62341tumour_sample1.pdf'), width=12, height=10)
DoHeatmap(tPD62341.1_filtered, features = top10$gene) + NoLegend()
dev.off()

# Dotplot of cluster markers (4 for each cluster)
tPD62341.1_filtered = FindClusters(tPD62341.1_filtered, resolution = 0.4, verbose = T)
tPD62341.1_filtered = RunUMAP(tPD62341.1_filtered, dims = 1:20, verbose = T, min.dist = 0.3, n.neighbors = 30) # with preferred settings
top4 = tPD62341.1_filtered.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC) %>% unique() # prevent duplicates
top4_markers = top4$gene %>% unique()
pdf(glue('Figures/F6/20241208_dotPlotTop4_PD62341tumour_sample1.pdf'), width=10, height=8)
DotPlot(tPD62341.1_filtered, features = top4_markers) + coord_flip()
dev.off()

# Dotplot of genes of interest 
pdf(glue('Figures/F6/20241208_dotPlot_PD62341tumour_sample1_GOI.pdf'), width=10, height=6)
DotPlot(object = tPD62341.1_filtered, features = genes_of_interest)
dev.off()

# Score by cell cycle expression (S and G2M genes accessed according to https://satijalab.org/seurat/articles/cell_cycle_vignette.html from October 2023)
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes
tPD62341.1_filtered = CellCycleScoring(tPD62341.1_filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# set.ident = TRUE sets the identity of the Seurat object to the cell-cycle phase

# Color UMAP by cell cycle phase (based on expression of cell cycle genes)
pdf('Figures/F6/20241208_umap_PD62341tumour_cellCycle_sample1.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "umap", label = TRUE, group.by = 'Phase')
dev.off() # don't see a separate clusters of proliferating cells :(
pdf('Figures/F6/20241208_pca_PD62341tumour_cellCycle_sample1.pdf', height = 4, width = 4)
DimPlot(tPD62341.1_filtered, reduction = "pca", label = TRUE, group.by = 'Phase')
dev.off() # kind of separate but not very clear? not sure if I should do anything about this 

###################################################################################################################################
# BASIC QC, PREPROCESSING, CLUSTERING - FOR EACH SAMPLE SEPARATELY 

# Create a function that takes the Seurat Object and carries out basic pre-processing, clustering etc.
sc_preprocessing = function(data, sample_name, mt_threshold, feature_threshold, count_threshold, n_pcs, resolution, min_dist, nn){
  
  # Required arguments
  # data = obtained with Read10X function
  # sample_name = name of the sample to be used in plots etc.
  # mt_threshold = % of mt reads above which a cell is excluded
  # feature_threshold = min nr of features (genes) identified in a cell
  # count_threshold = min nr of counts (molecules) identified in a cell
  
  # Initialize Seurat objects with the raw (non-normalized data)
  # Optional settings when reading in the file: 
  # min.cells = Include features detected in at least this many cells
  # min.features = Include cells where at least this many features are detected
  
  # read the data into the Seurat object 
  seurat = CreateSeuratObject(counts = data, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
  
  # calculate the percentage mt content 
  seurat[["percent.mt"]] = PercentageFeatureSet(seurat, pattern = "^MT-")
  
  # plot the distribution of the mt content in each sample 
  pdf(glue('Figures/F6/20241208_1_hist_mt_content_{sample_name}.pdf'), height = 4.5, width = 4.5)
  hist(seurat[["percent.mt"]] %>% unlist(), xlab = '% mt content', main = glue('{sample_name}'), breaks = 100)
  abline(v = mt_threshold, col = 'red', lwd = 2.5)   
  dev.off()
  
  # add CellID to the metadata 
  seurat@meta.data$CellID = rownames(seurat@meta.data) # add the cellID to the metadata
  meta = seurat@meta.data # create a dt with all the metadata of interest, including cell ID 
  
  # nCountRNA is the total number of molecules detected in each cell 
  pdf(glue('Figures/F6/20241208_1_hist_nCountRNA_{sample_name}.pdf'), height = 4, width = 8)
  ncounts_rna = as.numeric(seurat[["nCount_RNA"]] %>% unlist())
  par(mfrow=c(1,2))
  hist(ncounts_rna[ncounts_rna < 2500], xlab = 'nCountRNA', main = sample_name, breaks = 100)
  hist(ncounts_rna, xlab = 'nCountRNA', main = sample_name, breaks = 100)
  dev.off() # 1000 sounds reasonable based on the hist 
  
  # nFeatureRNA is the number of genes detected in each cell 
  pdf(glue('Figures/F6/20241208_1_hist_nFeatureRNA_{sample_name}.pdf'), height = 4, width = 8)
  nfeatures_rna = as.numeric(seurat[["nFeature_RNA"]] %>% unlist())
  par(mfrow=c(1,2))
  hist(nfeatures_rna[nfeatures_rna < 2000], xlab = 'nFeatureRNA', main = sample_name, breaks = 100)
  hist(nfeatures_rna, xlab = 'nFeatureRNA', main = sample_name, breaks = 100)
  dev.off() # 300 sounds reasonable based on the hist 

  # we can also plot those 3 QC parameters and save them to results files
  pdf(glue('Figures/F6/20241208_2_vlnPlotQC_{sample_name}.pdf'), height = 4, width = 8)
  print(VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off() # I don't see how those plots would help anyone choose the lower threshold for counts and features 
  
  # plot QC features against one another 
  pdf(glue('Figures/F6/20241208_3_featureScatter_{sample_name}.pdf'), height = 4, width = 8)
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
  genes_to_exclude = unique(c(grep('\\.[0-9]+_',rownames(seurat_filtered),value=TRUE), # poorly characterized (what does this mean??)
                              grep('MALAT1',rownames(seurat_filtered),value=TRUE), # contamination - why are we excluding this? how do we know this is contamination
                              grep('^HB[BGMQDAZE][12]?_',rownames(seurat_filtered),value=TRUE), # contamination (apparently you should filter out blood markers)
                              grep('^MT-',rownames(seurat_filtered),value=TRUE), # mt genes  
                              grep(excludeGenes,rownames(seurat_filtered),value=TRUE) # housekeeping genes
  ))
  keep_features=setdiff(rownames(seurat_filtered), genes_to_exclude)
  seurat_filtered = subset(seurat_filtered, features=keep_features)
  
  # run standard Seurat pipeline 
  seurat_filtered = NormalizeData(object = seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000) # default 1e4
  seurat_filtered = FindVariableFeatures(seurat_filtered, selection.method = "vst", nfeatures = 2000) # default 2k 
  seurat_filtered = ScaleData(seurat_filtered, vars.to.regress = "percent.mt") # Data scaling performed by default on the 2,000 variable features, some tutorials recommend regressing out mt content  
  seurat_filtered = RunPCA(seurat_filtered, npcs = n_pcs, features=VariableFeatures(object =seurat_filtered), verbose = F) # npcs = number of PCs to compute and store (50 by default)
  
  # Determine how many PCs to include in downstream analysis 
  pdf(glue('Figures/F6/20241208_4_elbowPlot_{sample_name}.pdf'), height = 4, width = 8)
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
  ggsave(glue('Figures/F6/20241208_4_varianceExplained_{sample_name}.pdf'), height = 4, width = 4)
  
  # Clustering 
  seurat_filtered = FindNeighbors(seurat_filtered, dims = 1:20, verbose = T)
  seurat_filtered = FindClusters(seurat_filtered, resolution = resolution, verbose = T)
  seurat_filtered = RunUMAP(seurat_filtered, dims = 1:20, verbose = T, min.dist = min_dist, n.neighbors = nn)
  
  # Plot outcomes of clustering 
  pdf(glue('Figures/F6/20241208_5_umap_{sample_name}.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "umap", label = TRUE))
  dev.off() 
  
  pdf(glue('Figures/F6/20241208_5_pca_{sample_name}.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "pca", label = TRUE))
  dev.off()
  
  # OUTPUT
  # Save the processed file with clusters
  saveRDS(seurat_filtered, file = glue("Out/F6/{sample_name}.rds"))
  
  # Identify features and look at markers 
  seurat_filtered.markers = FindAllMarkers(seurat_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # Show distribution of cluster markers
  top4 = seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
  top4_markers = top4$gene %>% unique()
  pdf(glue('Figures/F6/20241208_6_dotPlotTop4_{sample_name}.pdf'), width=10, height=6)
  print(DotPlot(tPD62341.1_filtered, features = top4_markers) + coord_flip())
  dev.off()
  
  # include genes that take part in the fusion, early embryonic, CD markers, synaptophysin, stuff Nathan looked at (to see if I am getting sth vaguely similar)
  genes_of_interest = c('MYOD1', 'PAX3', 'PAX7', 'CYP11A1', 'CD14', 'CD68', 'SYP', 'MN1', 'ZNF341', 'NR5A1', 'KCNQ1', 'TP53') # MYOG is 0 for all so I got rid of this
  pdf(glue('Figures/F6/20241208_6_featurePlot_{sample_name}.pdf'), width=12, height=8)
  print(FeaturePlot(seurat_filtered, features = genes_of_interest, order=TRUE))
  dev.off()
  
  # Heatmap of markers 
  top10 = seurat_filtered.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  pdf(glue('Figures/F6/20241208_7_heatmapTop10_{sample_name}.pdf'), width=12, height=8)
  print(DoHeatmap(seurat_filtered, features = top10$gene) + NoLegend())
  dev.off()
  
  # Dotplot of genes of interest 
  pdf(glue('Figures/F6/20241208_8_dotPlot_{sample_name}_GOI.pdf'), width=10, height=6)
  print(DotPlot(object = seurat_filtered, features = genes_of_interest))
  dev.off()
  
  # Score by cell cycle expression (S and G2M genes accessed according to https://satijalab.org/seurat/articles/cell_cycle_vignette.html from October 2023)
  s.genes = cc.genes$s.genes
  g2m.genes = cc.genes$g2m.genes
  seurat_filtered = CellCycleScoring(seurat_filtered, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  # set.ident = TRUE sets the identity of the Seurat object to the cell-cycle phase, I don't want that
  
  # Color UMAP by cell cycle phase (based on expression of cell cycle genes)
  pdf(glue('Figures/F6/20241208_9_umap_{sample_name}_cellCycle.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "umap", label = TRUE, group.by = 'Phase'))
  dev.off() # don't see a separate clusters of proliferating cells :(
  pdf(glue('Figures/F6/20241208_9_pca_{sample_name}_cellCycle.pdf'), height = 4, width = 4)
  print(DimPlot(seurat_filtered, reduction = "pca", label = TRUE, group.by = 'Phase'))
  dev.off() # kind of separate but not very clear? not sure if I should do anything about this 
  
}

###################################################################################################################################
# INTEGRATION OF SAMPLES FROM DIFFERENT LANES

# at what point should I integrate samples from different lanes?
# I think the most common way to do this is to integrate the .fastq reads at the stage of running CellRanger
# However, I also found this https://ucdavis-bioinformatics-training.github.io/2020-August-intro-scRNAseq/data_analysis/scRNA_Workshop-PART1_fixed#:~:text=Load%20the%20Cell%20Ranger%20Matrix,can%20aggregate%20them%20in%20R.
# which seems to suggest that you can integrate samples from multiple lanes when running Seurat (after CellRanger, working on CellRanger count matrix)

# Samples 
# I have the following samples (these are tumour samples!)
# PD62341 fresh: CG_SB_NB13652544, CG_SB_NB13652545, CG_SB_NB13652546
# PD63383 fresh: CG_SB_NB13760628, CG_SB_NB13760629, CG_SB_NB13760630
# PD63383 frozen (nuclear): CG_SB_NB14599068, CG_SB_NB14599069, CG_SB_NB14599070 - not used due to QC failure

PD62341_ids = c("CG_SB_NB13652544", "CG_SB_NB13652545", "CG_SB_NB13652546")
PD63383_ids = c("CG_SB_NB13760628", "CG_SB_NB13760629", "CG_SB_NB13760630")

data.loc = "/Users/bw18/Desktop/1SB"

PD62341.d = sapply(PD62341_ids, function(i){
  d10x = Read10X(file.path(data.loc, "scRNAseq/Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
PD62341.data = do.call("cbind", PD62341.d)

PD63383.d = sapply(PD63383_ids, function(i){
  d10x = Read10X(file.path(data.loc, "scRNAseq/Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})
PD63383.data = do.call("cbind", PD63383.d)

# Pre-processing of tumour data with specified parameters  
sc_preprocessing(PD62341.data, 'tPD62341.agg', 10, 300, 1000, 50, 0.2, 0.5, 50) 
sc_preprocessing(PD63383.data, 'tPD63383.agg', 10, 300, 1000, 50, 0.2, 0.5, 50) 

# Analyse datasets from both twins together 
tumour_sc_ids = c("CG_SB_NB13652544", "CG_SB_NB13652545", "CG_SB_NB13652546", "CG_SB_NB13760628", "CG_SB_NB13760629", "CG_SB_NB13760630")
tumour_sc.d = sapply(tumour_sc_ids, function(i){
  d10x = Read10X(file.path(data.loc, "scRNAseq/Data", i, "filtered_feature_bc_matrix"))
  colnames(d10x) = paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

tumour_sc.data = do.call("cbind", tumour_sc.d)
sc_preprocessing(tumour_sc.data, 'tumour_sc.agg', 10, 300, 1000, 50, 0.1, 0.5, 50)

###################################################################################################################################
# NUCLEAR RNA-seq 
# TL;DR conclusion from this analysis is that the samples failed (very poor quality, only 194 cells)
# therefore, I won't proceed with those samples in further analysis 

# read the data into the Seurat object 
tPD63383.nucl.1 = CreateSeuratObject(counts = tumourPD63383.1n, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
paste('Number of cells in the single nuclear dataset:', dim(tPD63383.nucl.1)[2]) # 194 so too few to be useful




