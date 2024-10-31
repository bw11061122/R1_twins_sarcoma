# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-10-25
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# Cleaned up script to identify mutations of interest 
# Script to look at tumour-specific mutations and tumour evolution 

###################################################################################################################################
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

###################################################################################################################################
# INPUT FILES 
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241030_413.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 255
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Add column to indicate chromosomes lost in the tumour
twins_filtered_dt[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'
))]

###################################################################################################################################
# PLOT SETTINGS
# Specify colors for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#a45f03"
col_PD63383 = "#d4b17e"
col_tumour_PD62341 = "#980505"
col_tumour_PD63383 = "#eb6767"
col_normal_PD62341 = "#0785a5"
col_normal_PD63383 = "#70c6db"
col_bar = '#8c08d3'

######################################################################################################
# SAMPLES
# create lists of possible samples of interest
samples_names = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
samples_normal = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341h", "PD62341n")
samples_tumour = c("PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq", "PD62341b", "PD62341u")
samples_PD62341 = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
samples_PD63383 = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383ap", "PD63383aq", "PD63383bb")
samples_normal_PD62341 = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD62341h", "PD62341n")
samples_normal_PD63383 = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb")
samples_tumour_PD62341 = c("PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD62341b", "PD62341u")
samples_tumour_PD63383 = c( "PD63383ap", "PD63383aq")

samples_mtr = grep("MTR", names(twins_dt), value = TRUE)
samples_dep = grep("DEP", names(twins_dt), value = TRUE)
samples_vaf = grep("VAF", names(twins_dt), value = TRUE)

samples_normal_mtr = paste(samples_normal, 'MTR', sep='_')
samples_normal_dep = paste(samples_normal, 'DEP', sep='_')
samples_normal_vaf = paste(samples_normal, 'VAF', sep='_')

samples_tumour_mtr = paste(samples_tumour, 'MTR', sep='_')
samples_tumour_dep = paste(samples_tumour, 'DEP', sep='_')
samples_tumour_vaf = paste(samples_tumour, 'VAF', sep='_')

samples_normal_PD62341_mtr = paste(samples_normal_PD62341, 'MTR', sep='_')
samples_normal_PD62341_dep = paste(samples_normal_PD62341, 'DEP', sep='_')
samples_normal_PD62341_vaf = paste(samples_normal_PD62341, 'VAF', sep='_')

samples_tumour_PD62341_mtr = paste(samples_tumour_PD62341, 'MTR', sep='_')
samples_tumour_PD62341_dep = paste(samples_tumour_PD62341, 'DEP', sep='_')
samples_tumour_PD62341_vaf = paste(samples_tumour_PD62341, 'VAF', sep='_')

samples_normal_PD63383_mtr = paste(samples_normal_PD63383, 'MTR', sep='_')
samples_normal_PD63383_dep = paste(samples_normal_PD63383, 'DEP', sep='_')
samples_normal_PD63383_vaf = paste(samples_normal_PD63383, 'VAF', sep='_')

samples_tumour_PD63383_mtr = paste(samples_tumour_PD63383, 'MTR', sep='_')
samples_tumour_PD63383_dep = paste(samples_tumour_PD63383, 'DEP', sep='_')
samples_tumour_PD63383_vaf = paste(samples_tumour_PD63383, 'VAF', sep='_')

samples_PD62341_mtr = paste(samples_PD62341, 'MTR', sep='_')
samples_PD62341_dep = paste(samples_PD62341, 'DEP', sep='_')
samples_PD62341_vaf = paste(samples_PD62341, 'VAF', sep='_')

samples_PD63383_mtr = paste(samples_PD63383, 'MTR', sep='_')
samples_PD63383_dep = paste(samples_PD63383, 'DEP', sep='_')
samples_PD63383_vaf = paste(samples_PD63383, 'VAF', sep='_')

######################################################################################################
# MTR, DEP, VAF DATA 

# subset MTR, DEP and VAF data
twins_filtered_mtr = twins_filtered_dt[,c('mut_ID', samples_mtr), with=FALSE]
twins_filtered_dep = twins_filtered_dt[,c('mut_ID', samples_dep), with=FALSE]
twins_filtered_vaf = twins_filtered_dt[,c('mut_ID', samples_vaf), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf
twins_filtered_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_mtr]
twins_filtered_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_mtr]
twins_filtered_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_mtr]
twins_filtered_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = samples_mtr]

twins_filtered_dep[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_filtered_dep[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_filtered_dep[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_filtered_dep[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]

twins_filtered_vaf[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_vaf]
twins_filtered_vaf[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_vaf]
twins_filtered_vaf[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_vaf]
twins_filtered_vaf[,max_vaf := apply(.SD, 1, max), .SDcols = samples_vaf]

# determine this specifically in normal samples (in case of ploidy changes in the tumour)
twins_filtered_mtr[,mean_mtr_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_mtr]
twins_filtered_mtr[,median_mtr_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_mtr]
twins_filtered_mtr[,min_mtr_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_mtr]
twins_filtered_mtr[,max_mtr_normal := apply(.SD, 1, max), .SDcols = samples_normal_mtr]

twins_filtered_dep[,mean_dep_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_dep]
twins_filtered_dep[,median_dep_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_dep]
twins_filtered_dep[,min_dep_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_dep]
twins_filtered_dep[,max_dep_normal := apply(.SD, 1, max), .SDcols = samples_normal_dep]

twins_filtered_vaf[,mean_vaf_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_vaf]
twins_filtered_vaf[,median_vaf_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_vaf]
twins_filtered_vaf[,min_vaf_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_vaf]
twins_filtered_vaf[,max_vaf_normal := apply(.SD, 1, max), .SDcols = samples_normal_vaf]

######################################################################################################
# PRESENCE / ABSENCE 

# Add presence / absence based on MTR data
twins_filtered_mtr[, sum_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_mtr[, sum_normal := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_mtr[, sum_PD62341 := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_mtr[, sum_PD63383 := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_mtr[, sum_tumour_PD62341 := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_mtr[, sum_tumour_PD63383 := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_mtr[, sum_normal_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_mtr[, sum_normal_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# Add presence / absence based on VAF data
twins_filtered_vaf[, sum_tumour := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_filtered_vaf[, sum_normal := rowSums(.SD>=0.1), .SDcols = samples_normal_vaf]
twins_filtered_vaf[, sum_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_PD62341_vaf]
twins_filtered_vaf[, sum_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_PD63383_vaf]
twins_filtered_vaf[, sum_tumour_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_filtered_vaf[, sum_tumour_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]
twins_filtered_vaf[, sum_normal_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_filtered_vaf[, sum_normal_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

######################################################################################################
# Overview: heatmaps 

# create matrix that stores MTR values and specify column and row names
mat_mtr = as.matrix(twins_filtered_mtr[, c(samples_mtr), with=FALSE])
rownames(mat_mtr) = twins_filtered_mtr[,1] %>% unlist()  
colnames(mat_mtr) = tstrsplit(colnames(mat_mtr), '_MTR', fixed=TRUE, keep=1) %>% unlist()
col_order = c('PD63383aq', 'PD63383ap', 'PD62341ae', 'PD62341ag', 'PD62341aj', 
              'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
              'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 
              'PD62341v', 'PD63383t', 'PD63383u', 'PD63383w', 'PD63383ae', 
              'PD63383ak', 'PD63383bb')
mat_mtr = mat_mtr[, col_order]

# annotation data frame 
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 12)), Twin = c(rep('PD63383', 2), rep('PD62341', 14), rep('PD63383', 6)))
rownames(col_annotation) = colnames(mat_mtr)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

pdf('Results/20241030_p2_heatmap_status_413mut_mtr.pdf')
pheatmap(mat_mtr,
         cellwidth=10, cellheight=0.5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="413 mutations", 
         legend = T,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

pdf('Results/20241030_p2_heatmap_status_413mut_mtr_clustered.pdf')
pheatmap(mat_mtr,
         cellwidth=10, cellheight=0.5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="413 mutations", 
         legend = T,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         treeheight_row = 0,
         fontsize=11, cexCol=2) 
dev.off()

# now show by presence / absence (MTR >= 4)
twins_mtr_binary = data.table(twins_filtered_mtr[, c(samples_mtr), with=FALSE])

for (j in names(twins_mtr_binary)){
  rows_to_0 = which(twins_mtr_binary[[j]] < 4)
  set(twins_mtr_binary, rows_to_0, j, 0)
  rows_to_1 = which(twins_mtr_binary[[j]] >= 4)
  set(twins_mtr_binary, rows_to_1, j, 1)
}

mat_mtr_binary = as.matrix(twins_mtr_binary)
colnames(mat_mtr_binary) = tstrsplit(colnames(mat_mtr_binary), '_MTR', fixed=TRUE, keep=1) %>% unlist()

pdf('Results/20241030_p2_heatmap_status_413mut_mtr_clustered_binary.pdf')
pheatmap(mat_mtr_binary,
         cellwidth=10, cellheight=0.5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="413 mutations MTR", 
         legend = F,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         treeheight_row = 0,
         fontsize=11, cexCol=2) 
dev.off()

# the same for VAF
mat_vaf = as.matrix(twins_filtered_vaf[, c(samples_vaf), with=FALSE])
rownames(mat_vaf) = twins_filtered_vaf[,1] %>% unlist()  
colnames(mat_vaf) = tstrsplit(colnames(mat_vaf), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_order = c('PD63383aq', 'PD63383ap', 'PD62341ae', 'PD62341ag', 'PD62341aj', 
              'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
              'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 
              'PD62341v', 'PD63383t', 'PD63383u', 'PD63383w', 'PD63383ae', 
              'PD63383ak', 'PD63383bb')
mat_vaf = mat_vaf[, col_order]

pdf('Results/20241030_p2_heatmap_status_413mut_vaf.pdf')
pheatmap(mat_vaf,
         cellwidth=10, cellheight=0.5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="413 mutations", 
         legend = T,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

pdf('Results/20241030_p2_heatmap_status_413mut_vaf_clustered.pdf')
pheatmap(mat_vaf,
         cellwidth=10, cellheight=0.5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="413 mutations", 
         legend = T,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         treeheight_row = 0,
         fontsize=11, cexCol=2) 
dev.off()

pdf('Results/20241030_p2_heatmap_status_413mut_vaf_clustered_rownames.pdf', height = 85)
pheatmap(mat_vaf,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="413 mutations", 
         legend = T,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         treeheight_row = 0,
         fontsize=11, cexCol=2) 
dev.off()

# converted to binary (presence / absence by VAF 0.1)
twins_vaf_binary = data.table(twins_filtered_vaf[, c(samples_vaf), with=FALSE])

for (j in names(twins_vaf_binary)){
  rows_to_0 = which(twins_vaf_binary[[j]] < 0.1)
  set(twins_vaf_binary, rows_to_0, j, 0)
  rows_to_1 = which(twins_vaf_binary[[j]] >= 0.1)
  set(twins_vaf_binary, rows_to_1, j, 1)
}

mat_vaf_binary = as.matrix(twins_vaf_binary)
colnames(mat_vaf_binary) = tstrsplit(colnames(mat_vaf_binary), '_VAF', fixed=TRUE, keep=1) %>% unlist()

pdf('Results/20241030_p2_heatmap_status_413mut_vaf_clustered_binary.pdf')
pheatmap(mat_vaf_binary,
         cellwidth=10, cellheight=0.5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="413 mutations: VAF", 
         legend = F,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         treeheight_row = 0,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Mutations that define each twin 

# PD62341 
twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 1
# chr14_105458006_C_A favourite one 
twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 1] # 0
twins_filtered_mtr[sum_normal_PD62341 == 6 & sum_normal_PD63383 == 2] # 3
# chr13_19104709_A_G # poor mapping 
# chr15_49480646_T_A # looks okay
# chr17_33422229_C_A # looks okay 

twins_filtered_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0
twins_filtered_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 1] # 4
# always missing in PD62341v and tend to be present in PD63383bb (contamination from the tumour?)
# chr16_5479739_C_T # looks real
# chr1_148178303_A_G # relatively poor mapping
# chr2_95662131_G_A # looks okay
# chr3_50106043_C_T # looks okay 
twins_filtered_mtr[sum_normal_PD62341 == 5 & sum_normal_PD63383 == 2] # 1
# chr1_25421823_C_T # rather poorly mapped (double check but v not convinced)

twins_filtered_mtr[sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 3
# chr20_54999342_C_A # looks okay (but not sure it's specific enough?)
# chr3_62055057_C_G # looks good 
# chr3_62055077_G_C # looks good 
twins_filtered_mtr[sum_normal_PD62341 == 4 & sum_normal_PD63383 == 1] # 3
# chr16_88009874_G_A # looks okay
# chr20_44114996_C_T # looks okay
# chrX_106521135_C_A # tricky - double check 
twins_filtered_mtr[sum_normal_PD62341 == 4 & sum_normal_PD63383 == 2] # 1
# chr9_23627511_T_G # looks okay

# private to PD63383
twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 6] # 0
twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 6] # 3
# all present in PD62341v (spleen?), all at VAF ~0.1?
# chr11_34011887_C_T looks okay (vaf ~0.25 in PD63383)
# chr4_75704880_G_A looks okay (vaf ~0.25 in PD63383)
# chr6_165179306_G_A looks okay (vaf ~0.25 in PD63383)
twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 6] # 1
# chr3_165901319_C_A # poor mapping (but not outrageous so double check!)
# chr6_55120674_A_G # not very believable, double check

twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 5] # 1
# chr13_50815806_A_G # looks okay (3 reads in PD63383u but VAF 0.15)
twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 5] # 2
# chr4_74625500_G_T # oh that looks great! 
# chr7_73831920_C_T # looks okay
twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 5] # 1
# chr9_94544944_G_C not perfect but seems okay - double check

twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 4] # 1
# chr21_28681933_A_C # variable mapping (checked Blat)
twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 4] # 2
# chr21_40193588_G_A # looks good but doesn't seem to be real PD63383-specific 
# chr3_77633967_C_T # looks okay (and present in v which kind of makes sense)
# chrX_115066661_C_T # looks okay
twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 4] # 1
# chr3_195014637_A_G does not look like a PD63383-specific mutation
twins_filtered_vaf[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 4] # 1
# chr2_95747706_T_C looks fine but quite high PD62341

######################################################################################################
######################################################################################################
# Plots (including only high-quality mutations with support from MTR + VAF)
muts_PD62341 = c("chr14_105458006_C_A",  'chr15_49480646_T_A', "chr17_33422229_C_A", 
                 "chr16_5479739_C_T", "chr2_95662131_G_A", "chr3_50106043_C_T",
                 "chr3_62055057_C_G", "chr3_62055077_G_C", "chr20_44114996_C_T") 
muts_PD63383 = c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                 "chr13_50815806_A_G", "chr4_74625500_G_T", "chr7_73831920_C_T",
                 "chr3_77633967_C_T", "chrX_115066661_C_T")
mut_early = c(muts_PD62341, muts_PD63383)

mut_PD62341_dt = twins_vaf[mut_ID %in% muts_PD62341, 1:23]
mut_PD62341_melt = melt(mut_PD62341_dt, id.vars = 'mut_ID')
mut_PD62341_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD62341_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_PD62341_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD62341_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD62341_melt[, mut_ID := factor(mut_ID, levels = 
                                     c('chr16_5479739_C_T','chr15_49480646_T_A',
                                       'chr17_33422229_C_A','chr14_105458006_C_A',  
                                       'chr3_50106043_C_T', 'chr2_95662131_G_A',
                                       "chr3_62055057_C_G", "chr3_62055077_G_C", 
                                       "chr20_44114996_C_T"))]

# plot option 1: we can make VAF plots for different samples
for (sample_name in mut_PD62341_melt[,sample] %>% unlist() %>% unique()){
  ggplot(mut_PD62341_melt[sample==sample_name], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD62341)+
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = 'Mutation', y = 'VAF')+
    ggtitle(glue('VAF distribution in sample {sample_name}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_{sample_name}.pdf'), width=9, height=7)
}

# show distribution of VAFs for each mutation in each sample 
ggplot(mut_PD62341_melt, aes(x=value, y=mut_ID, color=sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 14)+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color='black', size = 0.7)+
  geom_vline(xintercept = 0.25, color='black', size = 0.7)
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_samples_xy.pdf'), width=9, height=7)

# identify spleen samples 
mut_PD62341_melt[, sample_type2 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(status, twin, sep = '_')
))]

mut_PD62341_melt[, sample_type3 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  sample == 'PD63383bb', 'normal_PD63383_skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(status, twin, sep = '_')
))]

ggplot(mut_PD62341_melt, aes(x=value, y=mut_ID, color=sample_type2))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_normal_PD62341, 'orange', col_normal_PD63383, 'darkgreen', col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 14)+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color='black', size = 0.7)+
  geom_vline(xintercept = 0.25, color='black', size = 0.7)
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_samples_xy_labelspleen.pdf'), width=9, height=7)

ggplot(mut_PD62341_melt[status=='normal'], aes(x=value, y=mut_ID, color=sample_type2))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_normal_PD62341, '#C42F0F', col_normal_PD63383, '#F99A49'))+
  theme_bw(base_size = 14)+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color='black', size = 0.7)+
  geom_vline(xintercept = 0.25, color='black', size = 0.7)
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_samples_xy_normal_labelspleen.pdf'), width=9, height=7)

ggplot(mut_PD62341_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_normal_PD62341, '#C42F0F', col_normal_PD63383, '#F99A49'))+
  theme_bw(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF distribution for PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))+
  geom_hline(yintercept = 0.5, color='black', size = 0.3)+
  geom_hline(yintercept = 0.25, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_samples_normal_labelspleen.pdf'), width=9, height=7)

ggplot(mut_PD62341_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type3))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_normal_PD62341, '#C42F0F', col_normal_PD63383, '#A149F9', '#F99A49'))+
  theme_bw(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF distribution for PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))+
  geom_hline(yintercept = 0.5, color='black', size = 0.3)+
  geom_hline(yintercept = 0.25, color='black', size = 0.3)
  #scale_y_continuous(breaks = NULL)
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_samples_normal_label_spleen_skin.pdf'), width=9, height=7)

# can we aggregate samples by sample type and plot the mean VAF (+ standard deviation)
# we should use a different color for each category 
mut_PD62341_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD62341_melt, FUN=mean))
mut_PD62341_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr16_5479739_C_T','chr15_49480646_T_A',
                                       'chr17_33422229_C_A','chr14_105458006_C_A',  
                                       'chr3_50106043_C_T', 'chr2_95662131_G_A',
                                       "chr3_62055057_C_G", "chr3_62055077_G_C", 
                                       "chr20_44114996_C_T"))]

ggplot(mut_PD62341_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for twin-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_aggregated.pdf'), width=9, height=7)

# what if I do this without contaminated samples 
samples_contaminated = c('PD62341v', 'PD62341aa', 'PD62341h', 'PD63383bb')
mut_PD62341_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD62341_melt[!sample %in% samples_contaminated], FUN=mean))
mut_PD62341_agg[, mut_ID := factor(mut_ID, levels = 
                                     c('chr16_5479739_C_T','chr15_49480646_T_A',
                                       'chr17_33422229_C_A','chr14_105458006_C_A',  
                                       'chr3_50106043_C_T', 'chr2_95662131_G_A',
                                       "chr3_62055057_C_G", "chr3_62055077_G_C", 
                                       "chr20_44114996_C_T"))]

ggplot(mut_PD62341_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type, ncol = 1)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for twin-specific mutations\n(excluded contaminated samples)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_aggregated_rmcontamination.pdf'), width=5, height=8)

ggplot(mut_PD62341_agg, aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=3)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD62341 mutations\n(excluded contaminated samples)'))
ggsave(glue('Results/20241030_p3_vaf_dist_PD62341_mut9_aggregated_rmcontamination_reorderxy.pdf'), width=6.5, height=8)

# PD63383 
muts_PD63383 = c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                 "chr13_50815806_A_G", "chr4_74625500_G_T", "chr7_73831920_C_T",
                 "chr3_77633967_C_T", "chrX_115066661_C_T")

mut_PD63383_dt = twins_vaf[mut_ID %in% muts_PD63383, 1:23]
mut_PD63383_melt = melt(mut_PD63383_dt, id.vars = 'mut_ID')
mut_PD63383_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD63383_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
mut_PD63383_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD63383_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD63383_melt[, mut_ID := factor(mut_ID, levels = 
                                      c("chr4_75704880_G_A", "chr11_34011887_C_T",  "chr6_165179306_G_A",
                                        "chr4_74625500_G_T", "chrX_115066661_C_T", "chr13_50815806_A_G", 
                                        "chr7_73831920_C_T", "chr3_77633967_C_T"))]

# plot option 1: we can make VAF plots for different samples
for (sample_name in mut_PD63383_melt[,sample] %>% unlist() %>% unique()){
  ggplot(mut_PD63383_melt[sample==sample_name], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD63383)+
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = 'Mutation', y = 'VAF')+
    ggtitle(glue('VAF distribution in sample {sample_name}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_{sample_name}.pdf'), width=9, height=7)
}

# can we aggregate samples by sample type and plot the mean VAF (+ standard deviation)
# we should use a different color for each category 
mut_PD63383_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD63383_melt, FUN=mean))
mut_PD63383_agg[, mut_ID := factor(mut_ID, levels = 
                                      c("chr4_75704880_G_A", "chr11_34011887_C_T",  "chr6_165179306_G_A",
                                        "chr4_74625500_G_T", "chrX_115066661_C_T", "chr13_50815806_A_G", 
                                        "chr7_73831920_C_T", "chr3_77633967_C_T"))]
mut_PD63383_agg[, sample_type := factor(sample_type, levels = 
                                          c('normal_PD63383', 'normal_PD62341',
                                            'tumour_PD63383', 'tumour_PD62341'))]

ggplot(mut_PD63383_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD63383, col_normal_PD62341, col_tumour_PD63383, col_tumour_PD62341))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for PD63383 specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_aggregated.pdf'), width=9, height=7)

# what if I do this without contaminating samples
samples_contaminated = c('PD62341v', 'PD62341aa', 'PD62341h', 'PD63383bb')
mut_PD63383_agg = data.table(aggregate(value~mut_ID+sample_type, mut_PD63383_melt[!sample %in% samples_contaminated], FUN=mean))
mut_PD63383_agg[, mut_ID := factor(mut_ID, levels = 
                                     c("chr4_75704880_G_A", "chr11_34011887_C_T",  "chr6_165179306_G_A",
                                       "chr4_74625500_G_T", "chrX_115066661_C_T", "chr13_50815806_A_G", 
                                       "chr7_73831920_C_T", "chr3_77633967_C_T"))]
mut_PD63383_agg[, sample_type := factor(sample_type, levels = 
                                          c('normal_PD63383', 'normal_PD62341',
                                            'tumour_PD63383', 'tumour_PD62341'))]

ggplot(mut_PD63383_agg, aes(x=mut_ID, y=value, col = sample_type))+
  geom_line(group = 1)+
  geom_point(size=1.5)+
  geom_area(aes(fill = mut_ID))+
  scale_color_manual(values = c(col_normal_PD63383, col_normal_PD62341, col_tumour_PD63383, col_tumour_PD62341))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  facet_wrap(~sample_type, ncol = 1)+
  theme(legend.position = "none")+
  ggtitle(glue('VAF distribution for PD63383 specific mutations\n(excluded contaminated samples)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_aggregated_rmcontamination.pdf'), width=5, height=8)

ggplot(mut_PD63383_agg, aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=3)+
  scale_color_manual(values = c(col_normal_PD63383, col_normal_PD62341, col_tumour_PD63383, col_tumour_PD62341))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD63383 mutations\n(excluded contaminated samples)'))+
  xlim(0, 0.4)+
  geom_vline(xintercept = 0.25, color = 'black', size = 0.7)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_aggregated_rmcontamination_reorderxy.pdf'), width=6.5, height=8)

# plot a distribution of PD63383 mutations in each sample 
mut_PD63383_melt[, mut_ID := factor(mut_ID, levels = 
                                     c("chr4_75704880_G_A", "chr11_34011887_C_T",  "chr6_165179306_G_A",
                                       "chr4_74625500_G_T", "chrX_115066661_C_T", "chr13_50815806_A_G", 
                                       "chr7_73831920_C_T", "chr3_77633967_C_T"))]
mut_PD63383_melt[, sample_type := factor(sample_type, levels = 
                                           c('normal_PD63383', 'normal_PD62341',
                                             'tumour_PD63383', 'tumour_PD62341'))]

ggplot(mut_PD63383_melt, aes(x=mut_ID, y=value, col = sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_normal_PD63383, col_normal_PD62341, col_tumour_PD63383, col_tumour_PD62341))+
  theme_bw(base_size = 14)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF distribution of PD63383-restricted mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.4))
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_dist_in_samples.pdf'), width=9, height=7)

ggplot(mut_PD63383_melt, aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=1.5)+
  scale_color_manual(values = c(col_normal_PD63383, col_normal_PD62341, col_tumour_PD63383, col_tumour_PD62341))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD63383 mutations\n(excluded contaminated samples)'))+
  xlim(0, 0.4)+
  geom_vline(xintercept = 0.25, color = 'black', size = 0.7)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_reorderxy.pdf'), width=6.5, height=8)

ggplot(mut_PD63383_melt[status=='normal'], aes(x=value, y=mut_ID, col = sample_type))+
  geom_point(size=1.5)+
  scale_color_manual(values = c(col_normal_PD63383, col_normal_PD62341, col_tumour_PD63383, col_tumour_PD62341))+
  theme_bw(base_size = 12)+
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black"))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD63383 mutations\n(excluded contaminated samples)'))+
  xlim(0, 0.4)+
  geom_vline(xintercept = 0.25, color = 'black', size = 0.7)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_reorderxy_normal.pdf'), width=6.5, height=8)

###################################################################################################################################
# Can we try to address the issue of possible twin-twin transfusion via the spleen?

# if this occurs, the spleen in PD62341 should be more similar to PD63383 than other PD62341 tissues
# at the same time, we expect spleen in PD63383 to be more similar to PD62341 than other PD63383 tissues 
# NB this is complicated because the tumour (which arose in PD62341) contaminated the skin sample (PD63383bb)

# check PD62341v and PD63383w in PD63383-specific mutations
mut_PD63383_melt[, sample_type := as.factor(paste(
  status, twin, sep = '_'))]

mut_PD63383_melt[, sample_type2 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(status, twin, sep = '_')
))]

mut_PD63383_melt[, sample_type3 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  sample == 'PD63383bb', 'normal_PD63383_skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(status, twin, sep = '_')
))]

ggplot(mut_PD63383_melt[status=='normal'], aes(x=value, y=mut_ID, color=sample_type2))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_normal_PD62341, '#C42F0F', col_normal_PD63383, '#F99A49'))+
  theme_bw(base_size = 14)+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF distribution for PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.4))+
  geom_vline(xintercept = 0.25, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_samples_xy_normal_labelspleen.pdf'), width=9, height=7)
# PD62341 spleen has more PD63383-specific mutations than other tissues
# PD63383 spleen does not seem to be 'special' among PD63383 tissues (not more similar to PD62341, ie not lower VAF)

ggplot(mut_PD63383_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_normal_PD62341, '#C42F0F', col_normal_PD63383, '#F99A49'))+
  theme_bw(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF distribution for PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))+
  geom_hline(yintercept = 0.5, color='black', size = 0.3)+
  geom_hline(yintercept = 0.25, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_mut8_samples_normal_labelspleen.pdf'), width=9, height=7)

ggplot(mut_PD63383_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type3))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_normal_PD62341, '#C42F0F', col_normal_PD63383, '#A149F9', '#F99A49'))+
  theme_bw(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF distribution for PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))+
  geom_hline(yintercept = 0.5, color='black', size = 0.3)+
  geom_hline(yintercept = 0.25, color='black', size = 0.3)
#scale_y_continuous(breaks = NULL)
ggsave(glue('Results/20241030_p3_vaf_dist_PD633831_mut8_samples_normal_label_spleen_skin.pdf'), width=9, height=7)

# ok this is as a thing to think about but I don't know how to proceed on this

###################################################################################################################################
# Rerun search for mutations excluding spleen samples due to the possibility of transfusion between twins 

# removed PD62341v and PD63383w, which are spleen samples from each twin
samples_normal_nospleen = c("PD62341q", "PD62341aa", "PD62341ad", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341h", "PD62341n")
samples_normal_ns_PD62341 = c("PD62341q", "PD62341aa", "PD62341ad", "PD62341h", "PD62341n")
samples_normal_ns_PD63383 = c("PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb")

samples_normal_ns_mtr = paste(samples_normal_nospleen, 'MTR', sep = '_')
samples_normal_ns_vaf = paste(samples_normal_nospleen, 'VAF', sep = '_')
samples_normal_ns_dep = paste(samples_normal_nospleen, 'DEP', sep = '_')

samples_normal_ns_PD62341_mtr = paste(samples_normal_ns_PD62341, 'MTR', sep = '_')
samples_normal_ns_PD62341_vaf = paste(samples_normal_ns_PD62341, 'VAF', sep = '_')
samples_normal_ns_PD62341_dep = paste(samples_normal_ns_PD62341, 'DEP', sep = '_')

samples_normal_ns_PD63383_mtr = paste(samples_normal_ns_PD63383, 'MTR', sep = '_')
samples_normal_ns_PD63383_vaf = paste(samples_normal_ns_PD63383, 'VAF', sep = '_')
samples_normal_ns_PD63383_dep = paste(samples_normal_ns_PD63383, 'DEP', sep = '_')

# subset MTR, DEP and VAF data
twins_ns_mtr = twins_filtered_dt[,c('mut_ID', samples_normal_ns_mtr), with=FALSE]
twins_ns_dep = twins_filtered_dt[,c('mut_ID', samples_normal_ns_dep), with=FALSE]
twins_ns_vaf = twins_filtered_dt[,c('mut_ID', samples_normal_ns_vaf), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf
twins_ns_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_ns_mtr]
twins_ns_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_ns_mtr]
twins_ns_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_ns_mtr]
twins_ns_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = samples_normal_ns_mtr]

twins_ns_dep[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_ns_dep]
twins_ns_dep[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_ns_dep]
twins_ns_dep[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_ns_dep]
twins_ns_dep[,max_dep := apply(.SD, 1, max), .SDcols = samples_normal_ns_dep]

twins_ns_vaf[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_ns_vaf]
twins_ns_vaf[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_ns_vaf]
twins_ns_vaf[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_ns_vaf]
twins_ns_vaf[,max_vaf := apply(.SD, 1, max), .SDcols = samples_normal_ns_vaf]

# Add presence / absence based on MTR data
twins_ns_mtr[, sum_normal := rowSums(.SD>=4), .SDcols = samples_normal_ns_mtr]
twins_ns_mtr[, sum_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_ns_PD62341_mtr]
twins_ns_mtr[, sum_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_ns_PD63383_mtr]

# Add presence / absence based on VAF data
twins_ns_vaf[, sum_tumour := rowSums(.SD>=0.1), .SDcols = samples_normal_ns_vaf]
twins_ns_vaf[, sum_PD62341 := rowSums(.SD>=0.1), .SDcols = samples_normal_ns_PD62341_vaf]
twins_ns_vaf[, sum_PD63383 := rowSums(.SD>=0.1), .SDcols = samples_normal_ns_PD63383_vaf]

###################################################################################################################################
# Identify mutations present only in one twin 

# specific to PD62341
twins_ns_mtr[sum_PD62341 == 5 & sum_PD63383 == 0] 
# chr14_105458006_C_A

twins_ns_mtr[sum_PD62341 == 4 & sum_PD63383 == 0] 
# chr20_54999342_C_A # looks good 
# chr3_62055057_C_G
# chr3_62055077_G_C

twins_ns_mtr[sum_PD62341 == 3 & sum_PD63383 == 0] 
# chr7_139659050_G_A # looks good 

twins_ns_mtr[sum_PD62341 == 2 & sum_PD63383 == 0] 
# "chr10_54938275_C_A" # insertions: double check
# "chr13_18565005_G_A" # poor mapping 
# "chr18_64267234_G_A" # variable 
# "chr18_71485446_G_A" # looks okay?
# "chr21_45729938_T_A" # looks okay
# "chr2_119170376_C_T" # looks okay
# "chr6_95827754_C_A" # variable mapping but looks okay

# specific to PD63383
twins_ns_mtr[sum_PD62341 == 0 & sum_PD63383 == 5] 
# chr11_34011887_C_T # looks okay
# chr4_75704880_G_A # looks okay
# chr6_165179306_G_A # looks okay

twins_ns_mtr[sum_PD62341 == 0 & sum_PD63383 == 4] 
# chr13_50815806_A_G # looks okay
# chr21_40193588_G_A # looks okay
# chr4_74625500_G_T # looks okay
# chr7_73831920_C_T # looks okay
# chrX_115066661_C_T # looks okay

twins_ns_mtr[sum_PD62341 == 0 & sum_PD63383 == 3] 
# chr12_18986629_G_A # not sure, double check
# chr21_28681933_A_C # poor mapping 
# chr3_77633967_C_T # looks okay
# chr7_149688370_C_T # looks okay
# chr7_63778593_T_G # poor mapping 

twins_ns_mtr[sum_PD62341 == 0 & sum_PD63383 == 2] 
# chr15_64724513_A_T # looks okay but feels like present in ~everything
# chr1_103587565_A_C # looks okay
# chr2_113525233_T_C # poor mapping 
# chr4_15905566_C_T # looks okay
# chr7_110349267_G_A # double check

###################################################################################################################################
# Presence / absence of mutations in either of the twins by VAF

# specific to PD62341
twins_ns_vaf[sum_PD62341 == 5 & sum_PD63383 == 0] 
# chr14_105458006_C_A

twins_ns_vaf[sum_PD62341 == 4 & sum_PD63383 == 0] 
# chr20_54999342_C_A # looks good 
# chr3_62055057_C_G
# chr3_62055077_G_C

twins_ns_vaf[sum_PD62341 == 3 & sum_PD63383 == 0] 
# chr13_18565005_G_A # poor mapping 
# chr8_46678612_T_C # fine but also present in the tumour!

twins_ns_vaf[sum_PD62341 == 2 & sum_PD63383 == 0] 
# chr10_100754461_C_T # looks okay
# chr20_44114996_C_T # looks okay
# chr21_45729938_T_A # looks okay
# chr2_57814739_A_C # looks okay
# chr6_95827754_C_A # looks okay
# chr7_139659050_G_A # looks okay 

# specific to PD63383
twins_ns_vaf[sum_PD62341 == 0 & sum_PD63383 == 5] 
# chr11_34011887_C_T # looks okay
# chr13_50815806_A_G # looks okay
# chr4_75704880_G_A # looks okay

twins_ns_vaf[sum_PD62341 == 0 & sum_PD63383 == 4] 
# chr21_40193588_G_A # looks okay
# chr4_74625500_G_T # looks okay
# chr6_165179306_G_A # looks okay
# chr7_73831920_C_T # looks okay
# chrX_115066661_C_T # looks okay

twins_ns_vaf[sum_PD62341 == 0 & sum_PD63383 == 3] 
# chr21_28681933_A_C # poor mapping 
# chr3_77633967_C_T # looks okay
# chr7_149688370_C_T # looks okay

twins_ns_vaf[sum_PD62341 == 0 & sum_PD63383 == 2] 
# chr11_18922568_G_T 
# chr19_7878006_A_G # why is this not thrown out?
# chr1_103587565_A_C # looks okay
# chr22_43290224_C_T # looks okay
# chr7_63778593_T_G # poor mapping 

###################################################################################################################################
###################################################################################################################################

# Why are these mutations not filtered out?
# chr7_63778593_T_G # poor mapping but I think not 'bad enough' to be picked up 
# chr19_7878006_A_G # it seems like the deletion is not there in the indel file? indel_all[Chrom=='chr19' & pos < 7898006 & pos > 7868006]
# chr21_28681933_A_C # poor mapping - decent but maps well to many places 
# chr7_150298369_A_T # poor mapping - actually that doesn't look too bad
# chr8_83543_A_G # poor mapping 
# chr13_18565005_G_A # poor mapping 
# chr2_113525233_T_C # poor mapping 
# chr13_19104709_A_G # poor mapping 
# chr15_49480646_T_A # poor mapping 


