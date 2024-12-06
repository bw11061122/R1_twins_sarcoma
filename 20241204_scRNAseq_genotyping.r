###################################################################################################################################
# SCRIPT 9

# Script to analyse the scRNA-seq genotyping allele integrator output 
# 2024-12-04
# Barbara Walkowiak bw18

# INPUT: 
# 1 .tsv dataframes generated from alleleIntegrator
# allele integrator ran 04/12/2024; jobID 664854
# downloaded files from the farm via rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scRNA_genotyping/output/" /Users/bw18/Desktop/1SB/scRNAseq/Data/
# 2 data frame with final 255 mutations (ID + assignment to informative group on the phylogeny)
# 3 Seurat object with clusters (generated from tumour from either twin)

# OUTPUT:
# Plots with mutation identity 

###################################################################################################################################
# LIBRARIES 

# Load needed libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)
library(comprehenr) # list comprehension in R 
library(stringr)
library(beeswarm)
library(viridis)
library(grid)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_background = '#c8c8c8'

###################################################################################################################################
# INPUT FILES 

setwd('/Users/bw18/Desktop/1SB')

# load AlleleIntegrator output from all samples 
samples_scRNAseq = c('NB13652544', 'NB13652545', 'NB13652546',
                     'NB13760628', 'NB13760629', 'NB13760630',
                     'NB14599068', 'NB14599069', 'NB14599070')

ai_counts_dt = data.table()
for (sample in samples_scRNAseq){
  ai_counts = fread(paste0('scRNAseq/Data/', sample, '_somaticVar_scRNA_alleleCounts.tsv'), sep='\t', fill = TRUE, header=TRUE)
  ai_counts[, sample_ID := sample]
  ai_counts_dt = rbind(ai_counts_dt, ai_counts)
}

# load Seurat object with UMAP-based clustering 
tumour_PD62341_clusters = readRDS(file = "scRNAseq/Out/tPD62341.agg.rds")
tumour_PD63383_clusters = readRDS(file = "scRNAseq/Out/tPD63383.agg.rds")

# load list of final mutations used for the phylogeny
muts_assignment = fread('Data/20241205_muts_classes_255.csv', sep = ',', header = TRUE)

###################################################################################################################################
# Processing the AlleleIntegrator file 

# Number of cells in which a read spanning the mutation position was identified
paste('Number of cells with >= 1 read spanning mutation position:', length(ai_counts_dt[, barcode] %>% unlist() %>% unique()))
paste('Number of reads spanning mutation position (across all cells):', sum(ai_counts_dt[, Tot]))
paste('Number of cells with reads spanning > 1 mutant positions:', sum(duplicated(ai_counts_dt[, barcode])))

# Histogram of reads
mut_counts = data.table(table(ai_counts_dt[, Tot]))
pdf('scRNAseq/Results/20241205_genotyping_hist_Tot.pdf', height = 3, width = 4)
hist(ai_counts_dt[, Tot], xlab = 'Number of reads over a position per cell',
     main = 'Read number per cell') # one position in ~200 cells, interesting
dev.off()

table(ai_counts_dt[, Tot])

# Which samples do the cells come from?
table(ai_counts_dt[, sample_ID]) # identified cells from all samples 

# Add source of the sample (twin and sample type: fresh vs nuclear in PD63383)
ai_counts_dt[, twin := factor(fcase(
  sample_ID %in% c('NB13652544', 'NB13652545', 'NB13652546'), 'PD62341',
  sample_ID %in% c('NB13760628', 'NB13760629', 'NB13760630', 'NB14599068', 'NB14599069', 'NB14599070'), 'PD63383'
))]

ai_counts_dt[, sample := factor(fcase(
  sample_ID %in% c('NB13652544', 'NB13652545', 'NB13652546'), 'PD62341_fresh',
  sample_ID %in% c('NB13760628', 'NB13760629', 'NB13760630'), 'PD63383_fresh',
  sample_ID %in% c('NB14599068', 'NB14599069', 'NB14599070'), 'PD63383_nuclear'
))]

table(ai_counts_dt[, twin])
table(ai_counts_dt[, sample])

# Add mutation identity and features 
muts_assignment[, mut_ID0 := substr(mut_ID, 1, nchar(mut_ID)-4)]

# identify mutation ID 
ai_counts_dt[, Chrom := paste0('chr', chr)] # set the same chromosome name as mut_ID 
ai_counts_dt[, mut_ID0 := paste0(Chrom, '_', pos)]
ai_counts_dt = merge(ai_counts_dt, muts_assignment, by = 'mut_ID0')

# For each read, determine if it is mutant or wt
ai_counts_dt[, Ref := tstrsplit(mut_ID, '_', fixed=TRUE, keep=3)]
ai_counts_dt[, Alt := tstrsplit(mut_ID, '_', fixed=TRUE, keep=4)]
mut_cols = c('A', 'C', 'G', 'T')
ai_counts_dt[, Read := apply(ai_counts_dt[, c(mut_cols), with=FALSE], 1, function(row) {
  read = mut_cols[which(row > 0)]
  paste(read, collapse = ', ')
})]
ai_counts_dt[, status := factor(fcase(
  Read == Ref, 'wt',
  Read == Alt, 'mutant'
))]

# number of mutations which could possibly be detected (because there is a scRNA-seq spanning this position)
paste('Number of positions which have a somatic mutation with reads in scRNA-seq:', length(ai_counts_dt[, mut_ID] %>% unlist() %>% unique()))
# 93 (out of 255)

paste('Number of reads reporting wt variant:', sum(ai_counts_dt[status=='wt', Tot])) # 758 
paste('Number of reads reporting mutation:', sum(ai_counts_dt[status=='mutant', Tot])) # 25 
  
# is each mutation seen only once across all cells, or multiple times?
mut_counts = data.table(table(ai_counts_dt[, mut_ID]))
pdf('scRNAseq/Results/20241205_genotyping_hist_mutID.pdf', height = 5, width = 5)
hist(mut_counts[,N], xlab = 'Number of times mutation seen', main = 'Mutation occurrence', breaks = 40)
dev.off()

# for each mutation, calculate its VAF in the scRNA-seq data 
positions_identified = ai_counts_dt[, mut_ID] %>% unlist() %>% unique()
vafs_scrnaseq = data.table()
for (mut in positions_identified){
  dt = ai_counts_dt[mut_ID == mut]
  
  # get total nr reads, total nr of mutant reads across all scRNA-seq samples
  dep_all = dim(dt)[1]
  mtr_all = dim(dt[status=='mutant'])[1]
  
  # get total nr reads, total nr of mutant reads in PD62341 tumour
  dep_PD62341 = dim(dt[twin=='PD62341'])[1]
  mtr_PD62341 = dim(dt[status=='mutant'&twin=='PD62341'])[1]
  
  # get total nr reads, total nr of mutant reads in PD63383 tumour (this can have PD63383 normal cells)
  dep_PD63383 = dim(dt[twin=='PD63383'])[1]
  mtr_PD63383 = dim(dt[status=='mutant'&twin=='PD63383'])[1]
  
  vals = data.table(t(c(mut, mtr_all, dep_all, mtr_PD62341, dep_PD62341, mtr_PD63383, dep_PD63383)))
  vafs_scrnaseq = rbind(vafs_scrnaseq, vals)
}

setnames(vafs_scrnaseq, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7'), 
         c('mut_ID', 'mtr_all', 'dep_all', 'mtr_PD62341', 'dep_PD62341', 'mtr_PD63383', 'dep_PD63383'))
vafs_scrnaseq[, names(vafs_scrnaseq)[2:7] := lapply(.SD, as.numeric), .SDcols = !c(1)]
vafs_scrnaseq[, vaf_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all, mtr_all/dep_all)]
vafs_scrnaseq[, vaf_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all, mtr_all/dep_all)]
vafs_scrnaseq[, vaf_PD62341_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_PD62341, mtr_PD62341/dep_PD62341)]
vafs_scrnaseq[, vaf_PD62341_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_PD62341, mtr_PD62341/dep_PD62341)]
vafs_scrnaseq[, vaf_PD63383_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_PD63383, mtr_PD63383/dep_PD63383)]
vafs_scrnaseq[, vaf_PD63383_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_PD63383, mtr_PD63383/dep_PD63383)]

# add annotation - what kind of mutation is this
vafs_scrnaseq = merge(vafs_scrnaseq, muts_assignment, by = 'mut_ID')

# what mutations are those?
table(ai_counts_dt[, mut_class])
table(ai_counts_dt[status== 'wt', mut_class])
table(ai_counts_dt[status == 'mutant', mut_class]) # mostly tumour but you should get embryonic ones!

###################################################################################################################################
# Integrate analysis with clusters

# in which clusters are there cells with those reads? what about cells with mutant reads specifically?
meta_PD62341 = tumour_PD62341_clusters@meta.data
meta_PD63383 = tumour_PD63383_clusters@meta.data

# integrate with clusters
ai_counts_dt[, barcodeID :=tstrsplit(barcode, '-', fixed=TRUE, keep=1) ]
ai_counts_dt[, CellID := paste0(barcodeID, '-CG_SB_', sample_ID)]
ai_counts_dt$seurat_clusters_PD62341 = meta_PD62341$seurat_clusters[match(ai_counts_dt$cellID,meta_PD62341$CellID)]
ai_counts_dt$seurat_clusters_PD63383 = meta_PD63383$seurat_clusters[match(ai_counts_dt$cellID,meta_PD63383$CellID)]

# why are some cells apparently not captured? are they removed through QC?
# this affects 226/755 cells so weird if all removed via QC 

# EXAMPLE
# ATCCGAATCCAGTAGT-CG_SB_NB13652545 is the first barcode not assigned to either PD62341 / PD63383 clusters 
tumourPD62341.2 = Read10X(data.dir = "scRNAseq/Data/CG_SB_NB13652545/filtered_feature_bc_matrix/") # lane 2
tPD62341.2 = CreateSeuratObject(counts = tumourPD62341.2, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
tPD62341.2@meta.data$CellID = rownames(tPD62341.2@meta.data) 
"ATCCGAATCCAGTAGT-1" %in% tPD62341.2@meta.data$CellID # FALSE 
# so why does my scRNA-seq AlleleIntegrator output say otherwise?

# add information on variant presence in a given cluster 
# create a dataframe with cells assigned to PD62341 tumour clusters, report no read / wt read / mutant read
PD62341_reads = ai_counts_dt[sample=='PD62341_fresh' & !is.na(seurat_clusters_PD62341), c('CellID', 'status'), with=FALSE]
PD63383_reads = ai_counts_dt[sample=='PD63383_fresh' & !is.na(seurat_clusters_PD63383), c('CellID', 'status'), with=FALSE]

# remove duplicate cells because for the moment I don't know what to do with those
PD62341_reads = PD62341_reads[!duplicated(PD62341_reads[,CellID])]
PD63383_reads = PD63383_reads[!duplicated(PD63383_reads[,CellID])]

PD62341_barcodes_wt = PD62341_reads[status=='wt', CellID] %>% unlist()
PD62341_barcodes_mut = PD62341_reads[status=='mut', CellID] %>% unlist()
PD63383_barcodes_wt = PD63383_reads[status=='wt', CellID] %>% unlist()
PD63383_barcodes_mut = PD63383_reads[status=='mut', CellID] %>% unlist()

rownames(PD62341_reads) = PD62341_reads$CellID
tumour_PD62341_clusters = AddMetaData(tumour_PD62341_clusters, metadata = PD62341_reads)
rownames(PD63383_reads) = PD63383_reads$CellID
tumour_PD63383_clusters = AddMetaData(tumour_PD63383_clusters, metadata = PD63383_reads)

DimPlot(tumour_PD62341_clusters, reduction = "umap", group.by = "status", cols = c(col_tumour, col_normal, col_background))
DimPlot(tumour_PD63383_clusters, reduction = "umap", group.by = "status", cols = c(col_tumour, col_normal, col_background))

# OK, I'm afraid we will need some more umap customization 
umap_data = FetchData(tumour_PD62341_clusters, vars = c("umap_1", "umap_2", "status"))
umap_data = data.table(umap_data)
umap_data[, status2 := factor(fcase(
  status == 'mutant', 'mutant read',
  status == 'wt', 'wt read',
  is.na(status), 'no reads mapped'
))]
umap_data[, status2 := factor(status2, levels = c('no reads mapped', 'wt read', 'mutant read'))]
ggplot(setorder(umap_data, status2), aes(umap_1, umap_2, color = status2, order = status2))+
  geom_point()+
  scale_color_manual(values = c(col_background, col_normal, col_tumour))+
  theme_classic(base_size = 13)+
  coord_equal(ratio=1)+
  xlim(c(-20, 20))+
  ylim(c(-20, 20))+
  labs(title = 'Reads spanning mutant positions\nPD62341 tumour', x = 'UMAP 1', y = 'UMAP 2', col = 'Read')
ggsave('scRNAseq/Results/20241205_umap_all_mutations_PD62341tumour.pdf', height = 5, width = 5)
# it looks like there are way more than there actually are

# look at tumour mutations only 
# add information on variant presence in a given cluster 
# create a dataframe with cells assigned to PD62341 tumour clusters, report no read / wt read / mutant read
tumour_reads = ai_counts_dt[(!is.na(seurat_clusters_PD62341) | is.na(seurat_clusters_PD63383)) & mut_class %in% c('tumour', 'tumour_PD62341_only', 'tumour_PD63383_only'), c('CellID', 'status'), with=FALSE]
normal_reads = ai_counts_dt[(!is.na(seurat_clusters_PD62341) | is.na(seurat_clusters_PD63383)) & mut_class %in% c('PD62341_specific_embryonic', 'PD63383_specific_embryonic'), c('CellID', 'status'), with=FALSE]

# remove duplicate cells because for the moment I don't know what to do with those
tumour_reads = tumour_reads[!duplicated(tumour_reads[,CellID])]
normal_reads = normal_reads[!duplicated(normal_reads[,CellID])]

rownames(tumour_reads) = tumour_reads$CellID
tumour_PD62341_clusters = AddMetaData(tumour_PD62341_clusters, metadata = tumour_reads)

# OK, I'm afraid we will need some more umap customization 
umap_data = FetchData(tumour_PD62341_clusters, vars = c("umap_1", "umap_2", "status"))
umap_data = data.table(umap_data)
umap_data[, status2 := factor(fcase(
  status == 'mutant', 'mutant read',
  status == 'wt', 'wt read',
  is.na(status), 'no reads mapped'
))]
umap_data[, status2 := factor(status2, levels = c('no reads mapped', 'wt read', 'mutant read'))]
ggplot(setorder(umap_data, status2), aes(umap_1, umap_2, color = status2, order = status2))+
  geom_point()+
  scale_color_manual(values = c(col_background, col_normal, col_tumour))+
  theme_classic(base_size = 13)+
  coord_equal(ratio=1)+
  xlim(c(-20, 20))+
  ylim(c(-20, 20))+
  labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Tumour mutations positions\nPD62341 tumour', col = 'Read')
ggsave('scRNAseq/Results/20241205_umap_tumour_mutations.pdf', height = 5, width = 5)

rownames(normal_reads) = normal_reads$CellID
tumour_PD62341_clusters = AddMetaData(tumour_PD62341_clusters, metadata = normal_reads)
umap_data = FetchData(tumour_PD62341_clusters, vars = c("umap_1", "umap_2", "status"))
umap_data = data.table(umap_data)
umap_data[, status2 := factor(fcase(
  status == 'mutant', 'mutant read',
  status == 'wt', 'wt read',
  is.na(status), 'no reads mapped'
))]
umap_data[, status2 := factor(status2, levels = c('no reads mapped', 'wt read', 'mutant read'))]
ggplot(setorder(umap_data, status2), aes(umap_1, umap_2, color = status2, order = status2))+
  geom_point()+
  scale_color_manual(values = c(col_background, col_normal, col_tumour))+
  theme_classic(base_size = 13)+
  coord_equal(ratio=1)+
  xlim(c(-20, 20))+
  ylim(c(-20, 20))+
  labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Normal mutations positions\nPD62341 tumour', col = 'Read')
ggsave('scRNAseq/Results/20241205_umap_normal_mutations.pdf', height = 5, width = 5)
# okay can I check which clusters are those mutations in 



