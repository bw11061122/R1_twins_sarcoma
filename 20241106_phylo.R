###################################################################################################################################
# SCRIPT 4

# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-11-06
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# Cleaned up script to identify mutations of interest 
# Script to look at tumour-specific mutations and tumour evolution 

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

###################################################################################################################################
# INPUT FILES 

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241114_599.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 599
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Import dataframe with purity estimates
purity_dt = data.table(read.csv('Data/20241114_estimates_tumour_cont_27muts_median.csv'))

# Import list of driver genes (from Henry Lee-Six, 12/11/2024)
driver_genes_dt = data.table(read.csv('Data/HLS_fibromatoses_driver_list_with_fusions.csv', header=T))
driver_genes = driver_genes_dt[, gene] %>% unlist()

###################################################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#0ac368"
col_PD63383 = "#a249e8"
col_tumour_PD62341 = "#099272"
col_tumour_PD63383 = "#6F09D4"
col_normal_PD62341 = "#71D99B"
col_normal_PD63383 = "#C99DF6"
col_bar = '#e87811'


######################################################################################################
# SAMPLES

# Create lists of possible samples of interest
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
# CALCULATE VAF FOR AGGREGATED SAMPLES

# Aggregated VAF in tumour samples 
twins_filtered_dt[, dep_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_filtered_dt[, dep_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_filtered_dt[, dep_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_filtered_dt[, dep_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_filtered_dt[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]
twins_filtered_dt[, dep_all_tumour := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_filtered_dt[, dep_all := rowSums(.SD), .SDcols = samples_dep]

twins_filtered_dt[, mtr_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, mtr_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_filtered_dt[, mtr_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, mtr_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, mtr_all_normal := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_filtered_dt[, mtr_all_tumour := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, mtr_all := rowSums(.SD), .SDcols = samples_mtr]

twins_filtered_dt[, vaf_all_normal_PD62341 := mtr_all_normal_PD62341 / dep_all_normal_PD62341]
twins_filtered_dt[, vaf_all_normal_PD63383 := mtr_all_normal_PD63383 / dep_all_normal_PD63383]
twins_filtered_dt[, vaf_all_tumour_PD62341 := mtr_all_tumour_PD62341 / dep_all_tumour_PD62341]
twins_filtered_dt[, vaf_all_tumour_PD63383 := mtr_all_tumour_PD63383 / dep_all_tumour_PD63383]
twins_filtered_dt[, vaf_all_normal := mtr_all_normal / dep_all_normal]
twins_filtered_dt[, vaf_all_tumour := mtr_all_tumour / dep_all_tumour]
twins_filtered_dt[, vaf_all := mtr_all / dep_all]

# Add presence / absence based on MTR data
twins_filtered_dt[, sum_tumour_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, sum_normal_mtr := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_dt[, sum_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_dt[, sum_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_dt[, sum_tumour_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, sum_tumour_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, sum_normal_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, sum_normal_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# Add presence / absence based on VAF data
twins_filtered_dt[, sum_tumour_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_filtered_dt[, sum_normal_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_vaf]
twins_filtered_dt[, sum_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD62341_vaf]
twins_filtered_dt[, sum_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD63383_vaf]
twins_filtered_dt[, sum_tumour_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_filtered_dt[, sum_tumour_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]
twins_filtered_dt[, sum_normal_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_filtered_dt[, sum_normal_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

twins_filtered_dt[, sum_tumour_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_mtr', 'sum_tumour_vaf')]
twins_filtered_dt[, sum_normal_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_mtr', 'sum_normal_vaf')]
twins_filtered_dt[, sum_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_PD62341_mtr', 'sum_PD62341_vaf')]
twins_filtered_dt[, sum_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_PD63383_mtr', 'sum_PD63383_vaf')]
twins_filtered_dt[, sum_tumour_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_PD62341_mtr', 'sum_tumour_PD62341_vaf')]
twins_filtered_dt[, sum_tumour_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_PD63383_mtr', 'sum_tumour_PD63383_vaf')]
twins_filtered_dt[, sum_normal_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD62341_mtr', 'sum_normal_PD62341_vaf')]
twins_filtered_dt[, sum_normal_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD63383_mtr', 'sum_normal_PD63383_vaf')]

# Add a section for clean normal PD63383 samples (from script 2, evident that skin (PD63383bb) is contaminated)
samples_normal_PD63383_clean = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak")
samples_normal_PD63383_clean_dep = paste0(samples_normal_PD63383_clean, '_DEP')
samples_normal_PD63383_clean_mtr = paste0(samples_normal_PD63383_clean, '_MTR')
samples_normal_PD63383_clean_vaf = paste0(samples_normal_PD63383_clean, '_VAF')

twins_filtered_dt[, sum_normal_PD63383_clean_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_clean_mtr]
twins_filtered_dt[, sum_normal_PD63383_clean_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_clean_vaf]
twins_filtered_dt[, sum_normal_PD63383_clean_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD63383_clean_mtr', 'sum_normal_PD63383_clean_vaf')]

######################################################################################################
# Identify groups of mutations of interest 

# MUTATIONS OF INTEREST:
# Present in all samples
# Present in only a single sample
# Present in all normal samples 
# Tumour specific  
# Absent from PD63383 normal samples (PD62341 subclonal)

# Create lists of mutations 
muts = twins_filtered_dt[, mut_ID] %>% unlist()
muts_all_samples = twins_filtered_dt[sum_tumour_mtr_vaf==10 & sum_normal_mtr_vaf==12, mut_ID] %>% unlist()

muts_all_normal = twins_filtered_dt[sum_normal_mtr_vaf==12, mut_ID] %>% unlist()
muts_all_normal_nt = setdiff(muts_all_normal, muts_all_samples)

muts_all_tumour = twins_filtered_dt[sum_tumour_mtr_vaf==10, mut_ID] %>% unlist()
muts_all_tumour_nt = setdiff(muts_all_tumour, muts_all_samples)

muts_only_normal = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_mtr_vaf>1, mut_ID] %>% unlist() # exclude mutations present in a single sample
muts_only_normal_PD62341 = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_PD62341_mtr_vaf>1 & sum_normal_PD63383_mtr_vaf==0, mut_ID] %>% unlist() # exclude mutations present in a single sample
muts_only_normal_PD63383 = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_PD62341_mtr_vaf==0 & sum_normal_PD63383_mtr_vaf>1, mut_ID] %>% unlist() # exclude mutations present in a single sample

muts_only_tumour = twins_filtered_dt[sum_normal_mtr_vaf==0 & sum_tumour_mtr_vaf>1, mut_ID] %>% unlist()
muts_only_tumour_cont = twins_filtered_dt[sum_normal_mtr_vaf<=3 & sum_normal_mtr_vaf>=1 & sum_tumour_mtr_vaf>5, mut_ID] %>% unlist()
muts_tumour_PD63383 = twins_filtered_dt[sum_normal_PD63383_mtr_vaf<=1 & sum_normal_PD62341_mtr_vaf==0 & sum_tumour_PD63383_mtr_vaf==2, mut_ID] %>% unlist()
muts_only_tumour_PD63383 = setdiff(muts_tumour_PD63383, c(muts_only_tumour, muts_only_tumour_cont)) # 48 
muts_tumour_spec = c(muts_only_tumour, muts_only_tumour_cont, muts_only_tumour_PD63383) # 140

muts_one_sample_normal = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_mtr_vaf==1, mut_ID] %>% unlist()
muts_one_sample_tumour = twins_filtered_dt[sum_tumour_mtr_vaf==1 & sum_normal_mtr_vaf==0, mut_ID] %>% unlist()
muts_one_sample = c(muts_one_sample_normal, muts_one_sample_tumour)

muts_absent_PD63383_all = twins_filtered_dt[sum_normal_PD63383_clean_mtr_vaf==0, mut_ID] %>% unlist() # 282
muts_absent_PD63383 = setdiff(muts_absent_PD63383_all, c(muts_only_normal, muts_tumour_spec, muts_one_sample))

# Check number of mutations in each category 
paste('Mutations present in all samples:', length(muts_all_samples)) # 78

paste('Number of muts present in all normal samples:', length(muts_all_normal)) # 97 (632 - 535: makes sense)
paste('Number of muts present in all tumour samples:', length(muts_all_tumour)) # 144
paste('Number of muts present in all normal samples but not all tumour samples:', length(muts_all_normal_nt)) # 19 (97 - 78 = 19)
paste('Number of muts present in all tumour samples but not all normal samples:', length(muts_all_tumour_nt)) # 66 (144 - 78 = 66)

paste('Number of muts present in only normal samples:', length(muts_only_normal)) # 20
paste('Number of muts present in only PD62341 normal samples:', length(muts_only_normal_PD62341)) # 1
paste('Number of muts present in only PD63383 normal samples:', length(muts_only_normal_PD63383)) # 7

paste('Number of muts present in only tumour samples:', length(muts_only_tumour)) # 35
paste('Number of muts likely in only tumour samples:', length(muts_only_tumour_cont)) # 81
paste('Number of muts present in PD63383 tumour only:', length(muts_only_tumour_PD63383)) # 24 
paste('Number of muts likely to be tumour-specific:', length(c(muts_tumour_spec))) # 140

paste('Number of muts present in a single sample (normal):', length(muts_one_sample_normal)) # 38 
paste('Number of muts present in a single sample (tumour):', length(muts_one_sample_tumour)) # 99 
paste('Number of muts present in only a single sample:', length(c(muts_one_sample))) # 137

paste('Number of muts absent from PD63383 normal (not yet identified):', length(c(muts_absent_PD63383))) # 10

# Create list with all mutations that have been assigned to a group (groups)
muts_assigned = c(muts_all_samples, muts_all_normal_nt, muts_all_tumour_nt,
                  muts_only_normal, muts_tumour_spec, muts_one_sample, muts_absent_PD63383) %>% unique()
muts_unassigned = setdiff(muts, muts_assigned)

paste('Number of mutations assigned to a group:', length(muts_assigned)) # 449
paste('Number of mutations unassigned to a group:', length(muts_unassigned)) # 150 

######################################################################################################
# Heatmap to show all classes of mutations
mut_all_dt = twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE]
mut_all_mat = as.matrix(mut_all_dt[, c(samples_vaf), with=FALSE])
rownames(mut_all_mat) = mut_all_dt[,1] %>% unlist()  
colnames(mut_all_mat) = tstrsplit(colnames(mut_all_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)),
                            Purity = c(0, 0.1, 0.2, 0, 0.9, 0.3, 0.5, 0.6, 0.8, 0.8, 0, 0, 0, 0, 0, 0.9, 0.9, 0.3, 0.7, 0.4, 0, 0.5)) # fraction of tumour cells 
rownames(col_annotation) = colnames(mut_all_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), 
                         Twin = c(PD62341=col_PD62341, PD63383=col_PD63383),
                         Purity = colorRampPalette(c('#f2e7e7', 'darkred'))(11))

# heatmap
pdf('Results/20241114_p4_heatmap_all_mutations599.pdf')
pheatmap(mut_all_mat,
         cellwidth=10, cellheight=0.4,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations (599)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# save larger format to be able to see mutation names
pdf('Results/20241114_p4_heatmap_all_mutations599_large.pdf', height = 90)
pheatmap(mut_all_mat,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations (599)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# Heatmap: mutations assigned to a group
mut_assigned_dt = twins_filtered_dt[mut_ID %in% muts_assigned, c('mut_ID', samples_vaf), with=FALSE]
mut_assigned_mat = as.matrix(mut_assigned_dt[, c(samples_vaf), with=FALSE])
rownames(mut_assigned_mat) = mut_assigned_dt[,1] %>% unlist()  
colnames(mut_assigned_mat) = tstrsplit(colnames(mut_assigned_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

pdf('Results/20241114_p4_heatmap_assigned_mutations449.pdf')
pheatmap(mut_assigned_mat,
         cellwidth=10, cellheight=0.4,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Assigned mutations (449)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# Heatmap: mutations I didn't assign to any group
mut_unassigned_dt = twins_filtered_dt[!mut_ID %in% muts_assigned, c('mut_ID', samples_vaf), with=FALSE]
mut_unassigned_mat = as.matrix(mut_unassigned_dt[, c(samples_vaf), with=FALSE])
rownames(mut_unassigned_mat) = mut_unassigned_dt[,1] %>% unlist()  
colnames(mut_unassigned_mat) = tstrsplit(colnames(mut_unassigned_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

pdf('Results/20241114_p4_heatmap_unassigned_mutations150.pdf')
pheatmap(mut_unassigned_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Unassigned mutations (150)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

pdf('Results/20241114_p4_heatmap_unassigned_mutations150_large.pdf', height = 30)
pheatmap(mut_unassigned_mat,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Unassigned mutations (150)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Plot tri-nucleotide context of assigned and unassigned mutations

get_trinucs <- function(mybed, genome) {
  
  mybed$order <- 1:nrow(mybed) # order that rows are supplied in 
  gr <- GRanges(seqnames = mybed$Chrom, IRanges(start = mybed$Pos, width=1), ref=mybed$Ref, alt=mybed$Alt, order=mybed$order)
  # create a GRanges object with coordinates for the mutation
  
  if (all(substr(seqlevels(gr), 1, 3) != "chr")) {
    gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
  }
  
  seqinfo(gr) <- seqinfo(genome)[seqlevels(gr)]
  gr <- sort(gr)
  bases <- c("A", "C", "G", "T")
  trinuc_levels <- paste0(rep(bases, each = 16), rep(rep(bases, each = 4), 4), rep(bases, 16))
  get_trinuc <- function(seqname) {
    pos <- start(gr[seqnames(gr) == seqname])
    view <- Views(genome[[seqname]], start = pos - 1, end = pos + 1)
    ans <- factor(as.character(view), levels = trinuc_levels, labels = 1:64)
    return(as.numeric(ans))
  }
  trinuc <- sapply(seqlevels(gr), get_trinuc)
  gr$trinuc <- factor(unlist(trinuc, use.names = FALSE), levels = 1:64, labels = trinuc_levels)
  remove(trinuc)
  gr$REF <- gr$ref
  gr$ALT <- gr$alt
  gr$context <- gr$trinuc
  torc <- which(gr$ref %in% c("A", "G"))
  gr$REF[torc] <- as.character(reverseComplement(DNAStringSet(gr$REF[torc])))
  gr$ALT[torc] <- as.character(reverseComplement(DNAStringSet(gr$ALT[torc])))
  gr$context[torc] <- as.character(reverseComplement(DNAStringSet(gr$context[torc])))
  gr$class <- paste(gr$REF, gr$ALT, "in", gr$context, sep = ".")
  class_levels <- paste(rep(c("C", "T"), each = 48), rep(c("A", "G", "T", "A", "C", "G"), each = 16), "in", paste0(rep(rep(bases, each = 4), 6), rep(c("C", "T"), each = 48), rep(bases, 24)), sep = ".")
  gr$class <- factor(gr$class, levels = class_levels)
  
  # match them up with one another using the original order that I put in - I think that they may have been reshuffled.
  grdf <- as.data.frame(gr)
  grdf <- grdf[with(grdf, order(order)),]
  return(grdf$class)
} 

# Plot the trinucleotide context for unassigned mutations
mybed_fin = twins_dt[mut_ID %in% muts,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_fin = get_trinucs(mybed_fin, BSgenome.Hsapiens.UCSC.hg38)
dt_fin = twins_dt[mut_ID %in% muts]
dt_fin$trins_fin=trins_fin

# plot the distribution of different mutations across different contexts 
mut_sign_counts_fin = data.table(table(dt_fin[, trins_fin]))
setnames(mut_sign_counts_fin, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_fin[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_fin[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_fin[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_fin, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'All mutations (599)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_all_muts599.pdf', width = 7.5, height = 3.5)

# Plot the trinucleotide distribution for assigned mutations
mybed_asn = twins_dt[mut_ID %in% muts_assigned,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_asn = get_trinucs(mybed_asn, BSgenome.Hsapiens.UCSC.hg38)
dt_asn = twins_dt[mut_ID %in% muts_assigned]
dt_asn$trins_asn=trins_asn

# plot the distribution of different mutations across different contexts 
mut_sign_counts_asn = data.table(table(dt_asn[, trins_asn]))
setnames(mut_sign_counts_asn, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_asn[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_asn[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_asn[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_asn, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Assigned mutations (449)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_assigned_muts.pdf', width = 7.5, height = 3.5)

# Plot the trinucleotide context for unassigned mutations
mybed_usn = twins_dt[mut_ID %in% muts_unassigned,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_usn = get_trinucs(mybed_usn, BSgenome.Hsapiens.UCSC.hg38)
dt_usn = twins_dt[mut_ID %in% muts_unassigned]
dt_usn$trins_usn=trins_usn

# plot the distribution of different mutations across different contexts 
mut_sign_counts_usn = data.table(table(dt_usn[, trins_usn]))
setnames(mut_sign_counts_usn, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_usn[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_usn[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_usn[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_usn, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Unassigned mutations (150)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_unassigned_muts.pdf', width = 7.5, height = 3.5)

# Trinucleotide context for different groups
# Present in all normal samples
# Plot the trinucleotide context for unassigned mutations
mybed1 = twins_dt[mut_ID %in% muts_all_normal,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins1 = get_trinucs(mybed1, BSgenome.Hsapiens.UCSC.hg38)
dt1 = twins_dt[mut_ID %in% muts_all_normal]
dt1$trins1=trins1

# plot the distribution of different mutations across different contexts 
mut_sign_counts1 = data.table(table(dt1[, trins1]))
setnames(mut_sign_counts1, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts1[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts1[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts1[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts1, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations in all normal samples (97)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_allNormal_muts.pdf', width = 7.5, height = 3.5)

# Likely tumour specific 
mybed2 = twins_dt[mut_ID %in% muts_tumour_spec,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins2 = get_trinucs(mybed2, BSgenome.Hsapiens.UCSC.hg38)
dt2 = twins_dt[mut_ID %in% muts_tumour_spec]
dt2$trins2=trins2

# plot the distribution of different mutations across different contexts 
mut_sign_counts2 = data.table(table(dt2[, trins2]))
setnames(mut_sign_counts2, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts2[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts2[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts2[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts2, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations likely specific to the tumour (140)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_tumourSpec_muts.pdf', width = 7.5, height = 3.5)

# PD63383 tumour specific 
# Plot the trinucleotide context for unassigned mutations
mybed3 = twins_dt[mut_ID %in% muts_tumour_PD63383,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins3 = get_trinucs(mybed3, BSgenome.Hsapiens.UCSC.hg38)
dt3 = twins_dt[mut_ID %in% muts_tumour_PD63383]
dt3$trins3=trins3

# plot the distribution of different mutations across different contexts 
mut_sign_counts3 = data.table(table(dt3[, trins3]))
setnames(mut_sign_counts3, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts3[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts3[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts3[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts3, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations specific to PD63383 tumour (24)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_tumourPD63383_muts.pdf', width = 7.5, height = 3.5)

# PD63383 tumour specific 
# Plot the trinucleotide context for unassigned mutations
mybed4 = twins_dt[mut_ID %in% muts_one_sample,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins4 = get_trinucs(mybed4, BSgenome.Hsapiens.UCSC.hg38)
dt4 = twins_dt[mut_ID %in% muts_one_sample]
dt4$trins4=trins4

# plot the distribution of different mutations across different contexts 
mut_sign_counts4 = data.table(table(dt4[, trins4]))
setnames(mut_sign_counts4, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts4[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts4[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts4[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts4, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations present in a single sample (137)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_singleSample_muts.pdf', width = 7.5, height = 3.5)


######################################################################################################
# Inspect mutations identified in categories of interest 

muts_all_samples 

muts_all_normal_nt # 19
# "chr10_79512492_C_G" # poor quality 
# "chr13_38383628_T_C" # looks okay (some mis-mapping of minus strands and insertions) # 3 reads in PD63383aq
# "chr16_20575297_A_G" # poor quality 
# "chr16_88009831_A_G" # some insertions (VAF 0.09 in PD62341ap - has lots of insertion and only 3 reads with the mutation)
# "chr18_56506391_G_C" # insertions + split reads 
# "chr19_40853394_C_A" # poor quality mapping 
# "chr19_50081651_T_C" # poor mapping of (-) strand reads, PD62341 has 5 reads (3 on better-mapping reads)
# "chr1_103717717_G_A" # poor quality mapping 
# "chr1_13269066_T_C"  # poor quality mapping
# "chr1_21423048_C_A"  # poor quality mapping of (-) strand reads
# "chr1_38827952_C_A"  # loss of chr1 copy in PD62341am, PD63383ap, PD63383aq
# "chr1_47072883_A_T"  # poor quality mapping 
# "chr2_1440752_T_C"   # looks okay - 5 reads in PD62341ak and VAF 0.09, likely presence in all samples 
# "chr3_190409217_C_T" # looks okay - 0.064 in PD62341u, but high-quality reads, likely presence in all samples
# "chr6_94069453_A_T"  # all likely due to insertions
# "chr9_133062780_G_C" # poor mapping  
# "chr9_77762599_A_T"  # poor mapping 
# "chrX_140138202_C_T" # poor mapping of (-) strand reads (explains lower VAF in PD62341u)
# "chrX_66502763_G_A" # looks good (3 reads in am and ap, likely due to lower coverage in those samples)

muts_all_tumour_nt
# "chr10_100754461_C_T" 
# "chr11_114491131_C_T" 
# "chr11_134158143_C_T" 
# "chr11_4220977_T_C"   
# "chr11_4594537_A_G"  
# "chr13_42613199_A_G"  
# "chr13_57286266_C_T"  
# "chr13_60216504_A_C"  
# "chr14_104500332_C_T" 
# "chr14_105458006_C_A"
# "chr15_23462705_C_A"  
# "chr15_48353502_C_A"  
# "chr15_49480646_T_A"  
# "chr15_56691722_C_T"  
# "chr16_32621521_C_A" 
# "chr16_4003400_G_A"   
# "chr16_5479739_C_T"   
# "chr17_18402053_G_C"  
# "chr17_21307477_G_A"  
# "chr17_33422229_C_A" 
# "chr17_36165401_A_G"  
# "chr17_36209262_G_T"  
# "chr17_36388753_C_G"  
# "chr17_40061856_G_C"  
# "chr17_78256883_C_T" 
# "chr19_90678_G_A"     
# "chr1_149223310_G_A"  
# "chr1_153050724_T_C" 
# "chr1_157363251_C_T"  
# "chr1_213216928_T_C" 
# "chr1_69984580_G_A"   
# "chr21_10551234_A_C"  
# "chr22_17774668_G_A"  
# "chr22_30553442_G_T"  
# "chr2_133311125_G_T" 
# "chr2_199565129_G_A"  
# "chr2_95662131_G_A"   
# "chr3_50106043_C_T"   
# "chr3_62055057_C_G"   
# "chr3_62055077_G_C"  
# "chr4_147908994_T_C"  
# "chr4_178731458_A_G"  
# "chr5_131368207_T_C"  
# "chr5_136349883_C_T"  
# "chr5_176771316_G_A" 
# "chr5_176771327_T_A"  
# "chr5_179514426_A_G"  
# "chr5_37790015_G_A"   
# "chr5_44907911_G_A"   
# "chr6_101192680_A_G" 
# "chr6_32384367_G_A"   
# "chr6_95827754_C_A"   
# "chr7_102544208_T_C"  
# "chr7_24047089_G_C"   
# "chr7_55402890_G_T"  
# "chr7_64978551_T_C"   
# "chr7_76580290_G_T"   
# "chr7_76798258_T_C"   
# "chr7_77510284_A_G"   
# "chr8_142914761_A_G" 
# "chr8_7418832_G_C"   
# "chr9_11364991_T_A"   
# "chr9_42329754_C_T"   
# "chrX_66719643_C_G"   
# "chrX_68803487_C_T"  
# "chrX_798031_G_A" 

muts_only_normal # 20
# "chr11_34011887_C_T" # looks okay
# "chr13_24392029_C_G" # poor quality 
# "chr13_50815806_A_G" # looks okay
# "chr17_20476029_C_T" # poor mapping 
# "chr18_56262973_G_A" # poor quality (mapping reverse strands + indels)
# "chr1_103587565_A_C" # looks real (3 PD63383 samples)
# "chr21_28681933_A_C" # poor quality 
# "chr21_40193588_G_A" # looks okay
# "chr3_165901319_C_A" # to consider - mapping is quite okay but maps to several places on blat
# "chr3_77633967_C_T" # looks okay
# "chr4_15905566_C_T" # looks okay (present at low levels, but mapping is good quality, I do not note any insertions, no double-mapping on blat) 
# "chr4_74625500_G_T" # looks okay 
# "chr4_75704880_G_A" # looks okay
# "chr6_159851462_T_C" # poor quality (deletions) 
# "chr6_165179306_G_A" # looks okay
# "chr7_110349267_G_A" # plausible (+ strand reads map with poor quality)
# "chr7_143846161_G_A" # poor mapping 
# "chr7_149688370_C_T" # looks okay
# "chr7_73831920_C_T"  # looks okay
# "chrX_115066661_C_T" # looks okay 

######################################################################################################
# Accounting for contamination (use purity estimates)

# VAF_observed = fraction_normal * VAF_normal + fraction_tumour * VAF_tumour 
# VAF_tumour = (VAF_observed - fraction_normal * VAF_normal) / fraction_tumour
# VAF_normal = (VAF_observed - fraction_tumour * VAF_tumour) / fraction_normal

# fraction_normal = fraction of normal cells, fraction_tumour = fraction of tumour cells 
# VAF_normal = VAF in normal cells, VAF_tumour = VAF in tumour cells 

# first, I need to obtain VAF values for aggregated normal and tumour samples
twins_agg_vaf = merge(twins_filtered_vaf[, c('mut_ID', samples_vaf), with=FALSE], merge(twins_filtered_mtr[, c('mut_ID', samples_mtr), with=FALSE], twins_filtered_dep[, c('mut_ID', samples_dep), with=FALSE], by = 'mut_ID'), by = 'mut_ID')

twins_agg_vaf[, dep_normal_all := rowSums(.SD), .SDcols = samples_normal_dep]
twins_agg_vaf[, dep_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_agg_vaf[, dep_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_agg_vaf[, dep_tumour_all := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_agg_vaf[, dep_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_agg_vaf[, dep_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]

twins_agg_vaf[, mtr_normal_all := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_agg_vaf[, mtr_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_agg_vaf[, mtr_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_agg_vaf[, mtr_tumour_all := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_agg_vaf[, mtr_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_agg_vaf[, mtr_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]

twins_agg_vaf[, vaf_normal_all := mtr_normal_all / dep_normal_all]
twins_agg_vaf[, vaf_normal_PD62341 := mtr_normal_PD62341 / dep_normal_PD62341]
twins_agg_vaf[, vaf_normal_PD63383 := mtr_normal_PD63383 / dep_normal_PD63383]
twins_agg_vaf[, vaf_tumour_all := mtr_tumour_all / dep_tumour_all]
twins_agg_vaf[, vaf_tumour_PD62341 := mtr_tumour_PD62341 / dep_tumour_PD62341]
twins_agg_vaf[, vaf_tumour_PD63383 := mtr_tumour_PD63383 / dep_tumour_PD63383]

twins_filtered_vaf_adj = data.table(twins_filtered_vaf)
twins_filtered_vaf_adj = merge(twins_filtered_vaf_adj, twins_agg_vaf[, c('mut_ID', 'vaf_tumour_all', 'vaf_normal_PD62341'), with=FALSE], by = 'mut_ID')

for (sample in samples_normal){
  sample_vaf = paste0(sample, '_VAF')
  sample_diff = paste0(sample, '_adj')
  sample_purity = as.numeric(purity_dt[sample==sample_vaf, purity_est])
  twins_filtered_vaf_adj[, (sample_diff) := (get(sample_vaf) - vaf_tumour_all * (1-sample_purity)) / sample_purity]
}

for (sample in samples_tumour){
  sample_vaf = paste0(sample, '_VAF')
  sample_diff = paste0(sample, '_adj')
  sample_purity = as.numeric(purity_dt[sample==sample_vaf, purity_est])
  twins_filtered_vaf_adj[, (sample_diff) := (get(sample_vaf) - vaf_normal_PD62341 * (1-sample_purity)) / sample_purity]
}

samples_PD62341_adj = paste0(samples_PD62341, '_adj')
samples_PD63383_adj = paste0(samples_PD63383, '_adj')
samples_normal_adj = paste0(samples_normal, '_adj')
samples_normal_PD62341_adj = paste0(samples_normal_PD62341, '_adj')
samples_normal_PD63383_adj = paste0(samples_normal_PD63383, '_adj')
samples_tumour_adj = paste0(samples_tumour, '_adj')
samples_tumour_PD62341_adj = paste0(samples_tumour_PD62341, '_adj')
samples_tumour_PD63383_adj = paste0(samples_tumour_PD63383, '_adj')

twins_filtered_vaf_adj[, sum_normal_adj := rowSums(.SD > 0.1), .SDcols = samples_normal_adj]
twins_filtered_vaf_adj[, sum_normal_PD62341_adj := rowSums(.SD > 0.1), .SDcols = samples_normal_PD62341_adj]
twins_filtered_vaf_adj[, sum_normal_PD63383_adj := rowSums(.SD > 0.1), .SDcols = samples_normal_PD63383_adj]
twins_filtered_vaf_adj[, sum_tumour_adj := rowSums(.SD > 0.1), .SDcols = samples_tumour_adj]
twins_filtered_vaf_adj[, sum_tumour_PD62341_adj := rowSums(.SD > 0.1), .SDcols = samples_tumour_PD62341_adj]
twins_filtered_vaf_adj[, sum_tumour_PD63383_adj := rowSums(.SD > 0.1), .SDcols = samples_tumour_PD63383_adj]

twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 6 & sum_normal_PD63383_adj == 0] # 2
# chr14_105458006_C_A # looks okay
# chr17_33422229_C_A # looks okay
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 5 & sum_normal_PD63383_adj == 0] # 3
# chr16_5479739_C_T # looks okay
# chr2_95662131_G_A # looks okay
# chr3_50106043_C_T # looks okay
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 4 & sum_normal_PD63383_adj == 0] # 2
# chr3_62055057_C_G # looks okay
# chr3_62055077_G_C # looks okay
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 3 & sum_normal_PD63383_adj == 0] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 2 & sum_normal_PD63383_adj == 0] # 0 
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 1 & sum_normal_PD63383_adj == 0] # 11
# everything in h / aa at lower VAF - how do I know this is contamination or real? (NB some in h are at adjusted VAF 0.3, surely this has to be real?)
# chr10_100754461_C_T # looks okay
# chr11_134158143_C_T # looks okay
# chr13_60216504_A_C # looks okay
# chr16_32621521_C_A # looks okay
# chr17_78256883_C_T # looks okay
# chr1_69984580_G_A # looks okay
# chr22_30553442_G_T # looks okay
# chr2_133311125_G_T # looks okay
# chr4_147908994_T_C # looks okay
# chrX_66719643_C_G # looks okay
# chrX_68803487_C_T # looks okay

# okay that's what we expected that looks fine (still some contamination but low)
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 6] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 5] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 4] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 1] # 0
# chr15_56691722_C_T # looks okay (it doesn't come up before bc now adjusting in PD62341 worked)
# chr9_11364991_T_A # looks okay

# still some contamination but low  
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 6 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 5 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 4 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 3 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 2 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 1 & sum_normal_PD63383_adj == 1] # 0
# h and bb
# chr15_23462705_C_A # has mapping issues so this may screw VAF 
# chr5_136349883_C_T # looks okay

######################################################################################################
# Basic checks: number of mutations once adjusted for contamination 

# Create lists of mutations 
muts_all_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj==12 & sum_tumour_adj==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist())) 
muts_normal_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist())) 
muts_tumour_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_tumour_adj>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist()))
muts_normal_all_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()))
muts_tumour_all_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_tumour_adj==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()))
muts_tumour_only_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==0, mut_ID] %>% unlist()))
muts_tumour_all_only_adj = Reduce(intersect, list(muts_tumour_all_adj, muts_tumour_only_adj)) %>% unlist()
muts_normal_only_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_tumour_adj==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==0, mut_ID] %>% unlist()))
muts_normal_all_only_adj = Reduce(intersect, list(muts_normal_all_adj, muts_normal_only_adj)) %>% unlist()
muts_normal_all_nt_adj = setdiff(muts_normal_all_adj, muts_all_adj)
muts_tumour_all_nt_adj = setdiff(muts_tumour_all_adj, muts_all_adj)

paste('Mutations present in all samples:', length(muts_all_adj)) # 344
paste('Number of muts present in min 1 normal sample:', length(muts_normal_adj)) # 945
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour_adj)) # 984
paste('Number of muts present in all normal samples:', length(muts_normal_all_adj)) # 594
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_adj)) # 413
paste('Number of muts present in only tumour samples:', length(muts_tumour_only_adj)) # 33
paste('Number of muts present in all tumour samples and only tumour samples:', length(muts_tumour_all_only_adj)) # 1
paste('Number of muts present in only normal samples:', length(muts_normal_only_adj)) # 12
paste('Number of muts in all tumour but not all normal:', length(muts_tumour_all_nt_adj)) # 69
paste('Number of muts in all normal but not all tumour:', length(muts_normal_all_nt_adj)) # 250
paste('Number of muts present in all normal samples and max 1 tumour:', length(muts_normal_all_only)) # 0 # makes sense since tumour is from normal

######################################################################################################
# what is the correlation between adjusted and un-adjusted VAF?
for (sample in samples_names){
  sample_vaf = paste0(sample, '_VAF')
  sample_adj = paste0(sample, '_adj')
  vaf = twins_filtered_vaf_adj[, ..sample_vaf] %>% unlist()
  adj = twins_filtered_vaf_adj[, ..sample_adj] %>% unlist()
  dt = data.table(cbind(vaf, adj))
  ggplot(dt, aes(x = vaf, y = adj))+
    geom_point(size=1.5)+
    theme_classic(base_size = 12)+
    labs(x = 'VAF', y = 'Adjusted VAF')+
    ggtitle(glue('{sample}'))
  ggsave(glue('Results/20241109_p5_vaf_vs_adj_{sample}.pdf'), height = 3, width = 3)
}

######################################################################################################
# CLASSES OF MUTATIONS 1: MUTATIONS PRESENT IN ALL NORMAL SAMPLES 
# Inspect mutations that are present in all samples
# Why were these not called as germline?

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_all] %>% unlist(), 
     xlab = 'VAF (agg normal samples)', main = '632 muts: all normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_all]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_PD62341.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD62341] %>% unlist(), 
     xlab = 'VAF (agg normal PD62341)', main = '632 muts: PD62341 normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD62341]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_PD63383.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD63383] %>% unlist(), 
     xlab = 'VAF (agg normal PD63383)', main = '632 muts: PD63383 normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD63383]), col = 'purple', lwd = 2.5)
dev.off()

# conclusion: these are clearly not germline based on their VAF 
# are these artefacts or real?

# plot VAF in PD62341 vs PD63383 
# if real / informative, likely enriched just in one twin
ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_normal_PD63383))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'VAF (agg PD62341 normal)', y = 'VAF (agg PD63383 normal)')+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  ggtitle(glue('632 mutations (all normal samples)'))+
  coord_equal(ratio=1)
ggsave(glue('Results/20241109_vaf_normal_PD62341_vs_PD63383.pdf'), heigh = 4, width = 4)

# identify enriched mutations 

# log fold change 
twins_agg_vaf[, vaf_ratio := vaf_normal_PD62341 / vaf_normal_PD63383]
twins_agg_vaf[, vaf_ratio_log := log2(vaf_ratio)]

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_ratio_PD62341toPD63383.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio] %>% unlist(), 
     xlab = 'VAF ratio PD62341 / PD63383', main = '632 muts: VAF ratio', breaks=30)
dev.off()

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_ratio_PD62341toPD63383_log.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_log] %>% unlist(), 
     xlab = 'log2(VAF ratio PD62341 / PD63383)', main = '632 muts: log2(VAF ratio)', breaks=30)
abline(v = -0.5, col = 'red', lwd = 2.5)
abline(v = 0.5, col = 'red', lwd = 2.5)
dev.off()

# get a p-value for comparison (KS test)
# Note: I think it doesn't hurt to do this test, BUT importantly with 6 samples per group you have very little power 
ks_test <- function(a, b) {ks.test(a, b, alternative = c("two.sided"))$p.value}
twins_agg_vaf[, ks_vaf := apply(.SD, 1, function(row) {
  a = as.numeric(row[samples_normal_PD62341_vaf])
  b = as.numeric(row[samples_normal_PD63383_vaf])
  ks_test(a, b)}), .SDcols = c(samples_normal_PD62341_vaf, samples_normal_PD63383_vaf)]
twins_agg_vaf[, ks_vaf_log := -1 * log10(ks_vaf)] # -log10(p value) for better display

twins_agg_vaf[, ks_vaf_adj := p.adjust(ks_vaf, method = 'BH')] # using Benjamini-Hochberg to adjust for multiple testing 
twins_agg_vaf[, ks_vaf_adj_log := -1 * log10(ks_vaf_adj)] # -log10(p value) for better display

paste('Minimum p-value (2-sided KS test) for normal PD62341 and PD63383:', min(twins_agg_vaf[mut_ID %in% muts_normal_all, ks_vaf])) # 0.0021645023
twins_agg_vaf[mut_ID %in% muts_normal_all & ks_vaf < 0.05, mut_ID] # 20
# "chr10_47676090_A_G" 
# "chr12_31851985_G_A" 
# "chr12_8792660_G_A"  
# "chr15_20220521_G_A"
# "chr15_20262122_C_T" 
# "chr15_22581980_C_A" 
# "chr15_81721403_G_A" 
# "chr15_81761276_G_A"
# "chr15_82483067_G_A" 
# "chr17_22069771_A_G" 
# "chr17_63871774_G_A" 
# "chr18_14868156_C_T"
# "chr19_53041299_G_C" 
# "chr1_16864604_G_A"  
# "chr1_16912588_C_T"  
# "chr1_38827952_C_A" 
# "chr3_198145129_C_A" 
# "chr3_75677255_A_G"  
# "chr7_143850635_G_T" 
# "chr8_123201027_T_C"
              
twins_agg_vaf[, mut_cat := factor(fcase(
  vaf_normal_PD62341 > 0.45 & vaf_normal_PD63383 > 0.45, 'putative germline', 
  vaf_ratio_log < -0.5 & ks_vaf < 0.05, 'PD63383 enriched',
  vaf_ratio_log > 0.5 & ks_vaf < 0.05, 'PD62341 enriched',
  ks_vaf > 0.05 | (vaf_ratio_log > -0.5 & vaf_ratio_log < 0.5), 'no difference'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_cat), aes(x = vaf_ratio_log, y = ks_vaf_log, col = mut_cat))+
  geom_point(size=1.5)+
  theme_classic(base_size = 12)+
  labs(x = 'log2fc (VAF PD62341 / VAF PD63383)', y = '-log10(p value), 2 sided KS test', col = 'Mutation category')+
  scale_color_manual(values = c('grey', col_PD62341, col_PD63383, 'darkred'))+
  ggtitle(glue('632 mutations (all normal samples)'))
ggsave(glue('Results/20241109_volcano_normal_PD62341_PD63383_col_cat.pdf'), width = 6, height = 4)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_normal_PD63383, col = mut_cat))+
  geom_point(size=1.5, alpha = 1)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF agg PD62341 normal', y = 'VAF agg PD63383 normal', col = 'Mutation category')+
  scale_color_manual(values = c('grey', col_PD62341, col_PD63383, 'darkred'))+
  ggtitle(glue('632 mutations (all normal samples)'))
ggsave(glue('Results/20241109_vaf_normal_PD62341_PD63383_col_cat.pdf'), width = 6, height = 4)

# inspect those mutations in detail 
twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat == 'putative germline'] # 5
# "chr12_131428827_C_A" # mapping is funny (checked on blat - not great)
# "chr15_85186362_A_G"  # again, mapping is funny
# "chr17_21783609_C_A" # mapping issues 
# "chr5_127735154_A_G" # mapping issues 
# "chr7_4912517_G_A" # mapping issues 
twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat == 'PD62341 enriched'] # 3
# "chr15_82483067_G_A" # kind of okay? double check
# "chr1_16864604_G_A"  # poor mapping (reads map somewhere else)
# "chr8_123201027_T_C" # does indeed look real
twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat == 'PD63383 enriched'] # 1
# "chr1_38827952_C_A" # looks real 

# plot the distribution of VAFs across samples for the good looking mutations
twins_agg_vaf_sub = twins_agg_vaf[, c('mut_ID', samples_vaf), with=FALSE]
twins_agg_vaf_melt = melt(twins_agg_vaf_sub, id.vars = 'mut_ID')
twins_agg_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_agg_vaf_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v',
  'PD63383ae', 'PD63383ak', 'PD63383bb', 'PD63383t', 'PD63383u', 'PD63383w',
  'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u',
  'PD63383ap', 'PD63383aq'))]
twins_agg_vaf_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
twins_agg_vaf_melt[, status := factor(fcase(sample %in% samples_normal, 'normal', sample %in% samples_tumour, 'tumour'))]
twins_agg_vaf_melt[, sample_type := factor(fcase(status == 'tumour', 'tumour', 
                                                 status == 'normal' & twin == 'PD62341', 'PD62341 normal',
                                                 status == 'normal' & twin == 'PD63383', 'PD63383 normal'))]

ggplot(twins_agg_vaf_melt[mut_ID=='chr8_123201027_T_C'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ggtitle(glue('Mutation: chr8:123201027, T>C'))+
  ylim(c(0, 0.7))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr8_PD62341_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr15_82483067_G_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ggtitle(glue('Mutation: chr15:82483067, G>A'))+
  ylim(c(0, 0.7))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr15_PD62341_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr1_38827952_C_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr1_PD63383_enriched.pdf'), height = 3, width = 6.5)

# try a relaxed threshold (l2fc 0.5 is a lot at the end of the day)
twins_agg_vaf[, mut_cat_2 := factor(fcase(
  vaf_ratio_log < -0.3, 'PD63383 enriched',
  vaf_ratio_log > 0.3, 'PD62341 enriched',
  vaf_ratio_log > -0.3 & vaf_ratio_log < 0.3, 'no difference'))]

twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat_2 == 'PD62341 enriched' & ks_vaf < 0.05] # 6
# chr15_82483067_G_A # kind of okay? double check
# chr1_16864604_G_A # poor mapping (reads map somewhere else)
# chr8_123201027_T_C # does indeed look real

# chr10_47676090_A_G # not great 
# chr15_22581980_C_A # poor mapping 
# chr17_63871774_G_A # maybe? to check with Henry 
twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat_2 == 'PD63383 enriched' & ks_vaf < 0.05] # 1
# chr1_38827952_C_A 

# chr12_8792660_G_A # looks okay, but not sure the difference is huge
# chr15_20262122_C_T # looks okay 
# chr15_81721403_G_A # looks okay
# chr18_14868156_C_T # poor mapping (maps to several places in the genome)
# chr19_53041299_G_C # looks okay

ggplot(twins_agg_vaf_melt[mut_ID=='chr17_63871774_G_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutation: chr17:63871774, G>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr17_PD62341_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr12_8792660_G_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutation: chr12:8792660, G>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr12_PD63383_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr15_20262122_C_T'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutation: chr15:20262122, C>T'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr15_1_PD63383_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr15_81721403_G_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutation: chr15:81721403, G>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr15_2_PD63383_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr19_53041299_G_C'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutation: chr19:53041299, G>C'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr19_PD63383_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr4_189788936_T_C'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutation: chr4:189788936, T>C '))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr4_tumour_enriched.pdf'), height = 3, width = 6.5)

# Are other non-germline, non-enriched mutations artefacts?
# Found mutations on chr12 that all look great - what is going on?
mut_chr12_cluster = c("chr12_31851749_C_T",  "chr12_31851985_G_A",  "chr12_31895265_G_A",
                      "chr12_31852349_C_G",  "chr12_31852949_G_A",  "chr12_31856122_C_T", 
                      "chr12_31867266_G_A",  "chr12_31869715_A_G",  "chr12_31872355_C_T", 
                      "chr12_31879431_C_T",  "chr12_31880061_G_A",  "chr12_31883693_C_T", 
                      "chr12_31885515_G_A",  "chr12_31885689_A_G",  "chr12_31887208_C_T" )

ggplot(twins_agg_vaf_melt[mut_ID %in% mut_chr12_cluster], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutation on chr12 31,850,000-31,900,000'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.6))
ggsave(glue('Results/20241109_dist_samples_mut_chr12.pdf'), width = 4.5, height = 3)

pdf('Results/20241109_p5_hist_chr12_PD62341_PD63383.pdf')
par(mfrow = c(2, 1))
hist(twins_agg[mut_ID %in% mut_chr12_cluster, PD63383_normal_vaf], xlab = 'VAF', xlim = c(0, 1), main = 'PD63383 chr12', breaks = 6)
hist(twins_agg[mut_ID %in% mut_chr12_cluster, PD62341_normal_vaf], xlab = 'VAF', xlim = c(0, 1), main = 'PD62341 chr12', breaks = 6)
dev.off()

# plot the coverage vs VAF for these mutations
twins_agg_vaf[, mut_chr12 := factor(fcase(
  mut_ID %in% mut_chr12_cluster, 'segment chr12',
  !mut_ID %in% mut_chr12_cluster, 'other'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_chr12), aes(x = dep_normal_PD62341, y = vaf_normal_PD62341, col = mut_chr12))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred'))+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD62341 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD62341.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_chr12), aes(x = dep_normal_PD63383, y = vaf_normal_PD63383, col = mut_chr12))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred'))+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD63383 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD63383.pdf'), height=5, width=7)

# similar thing on chr 14?
mut_chr14_cluster = c("chr14_19725920_G_T",  "chr14_19735467_A_C",  "chr14_19770060_C_T",  "chr14_19778921_G_T", 
                      "chr14_19803896_T_A",  "chr14_19814466_C_T",  "chr14_19823326_A_T", "chr14_19880026_C_G", 
                      "chr14_19881983_C_T",  "chr14_19887911_C_T",  "chr14_19895402_G_A",  "chr14_19899009_A_G", 
                      "chr14_19899033_A_G",  "chr14_19903456_A_G")


twins_agg_vaf[, mut_cluster := factor(fcase(
  mut_ID %in% mut_chr12_cluster, 'segment chr12',
  mut_ID %in% mut_chr14_cluster, 'segment chr14',
  !mut_ID %in% c(mut_chr12_cluster, mut_chr14_cluster), 'other'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_cluster), aes(x = dep_normal_PD62341, y = vaf_normal_PD62341, col = mut_cluster))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred', '#0db793'))+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD62341 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD62341_clusters.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_cluster), aes(x = dep_normal_PD63383, y = vaf_normal_PD63383, col = mut_cluster))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred', '#0db793'))+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD63383 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD63383_clusters.pdf'), height=5, width=7)

# can I plot the distribution of these mutations for each chromosome
twins_agg_vaf[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
twins_agg_vaf[, pos := tstrsplit(mut_ID, '_', fixed=TRUE, keep=2)]
twins_agg_vaf[, pos := as.numeric(pos)]

Chrom = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
          "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
          "chr20", "chr21", "chr22", "chrX")
Chrom_length = as.numeric(c(248956422, 242193529, 198295559, 190214555,
                            181528259, 170805979, 159345973, 145138636,
                            139394717, 133797422, 135086622, 133275309,
                            114364328, 107043718, 101991189, 90338345,
                            83257441, 80373285, 58617616, 64444167,
                            46709983, 50818468, 150040895))
chr_lengths_dt = data.table(cbind(Chrom, Chrom_length))

twins_agg_vaf = merge(twins_agg_vaf, chr_lengths_dt, by = 'Chrom')
twins_agg_vaf[, Chrom := factor(Chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                              "chr20", "chr21", "chr22", "chrX"))]

for (chr in Chrom){
  dt = twins_agg_vaf[Chrom == chr & mut_ID %in% muts_normal_all]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = vaf_normal_PD62341))+
    geom_point(size=2.5)+
    theme_classic(base_size = 12)+
    scale_color_manual(values = '#7b13d4')+
    labs(x = 'Genomic position', y = 'VAF (PD62341 normal)')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241109_dist_across_genome_632_{chr}_PD62341.pdf'), height=3, width=4.5)
}

for (chr in Chrom){
  dt = twins_agg_vaf[Chrom == chr & mut_ID %in% muts_normal_all]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = vaf_normal_PD63383))+
    geom_point(size=2.5)+
    theme_classic(base_size = 12)+
    scale_color_manual(values = '#7b13d4')+
    labs(x = 'Genomic position', y = 'VAF (PD63383 normal)')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241109_dist_across_genome_632_{chr}_PD63383.pdf'), height=3, width=4.5)
}

# compare to coverage across all mutations in the set
twins_agg_vaf[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  mut_ID %in% muts_normal_only, 'only normal samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_all, muts_normal_only, muts_tumour_all_only, muts_tumour_only), 'other'))]

ggplot(twins_agg_vaf, aes(x = dep_normal_PD62341, y = vaf_normal_PD62341, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Mutation class')+
  ggtitle(glue('PD62341 normal samples (1002 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002_PD62341.pdf'), height=5, width=7)

ggplot(twins_agg_vaf, aes(x = dep_normal_PD63383, y = vaf_normal_PD63383, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Mutation class')+
  ggtitle(glue('PD63383 normal samples (1002 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002_PD63383.pdf'), height=5, width=7)

# Coverage histograms across different mutation classes
twins_dt[, sum_normal_PD63383_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_dt[, sum_normal_PD63383_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_dt[, sum_normal_PD62341_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_dt[, sum_normal_PD62341_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_dt[, sum_normal_PD63383_vaf := sum_normal_PD63383_mtr / sum_normal_PD63383_dep]
twins_dt[, sum_normal_PD62341_vaf := sum_normal_PD62341_mtr / sum_normal_PD62341_dep]
twins_dt[, sum_normal_mtr := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_dt[, sum_normal_dep := rowSums(.SD), .SDcols = samples_normal_dep]
twins_dt[, sum_normal_vaf := sum_normal_mtr / sum_normal_dep]

par(mfrow = c(1,1))

pdf('Results/20241109_p6_hist_PD62341_allmuts.pdf')
hist(twins_dt[sum_normal_dep < 600, sum_normal_PD62341_dep], 
     xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, all mutations (362,814)', breaks = 50)
abline(v = median(twins_dt[sum_normal_dep < 600, sum_normal_PD62341_dep]), col = 'red', lwd = 2)
dev.off()

pdf('Results/20241109_p6_hist_PD63383_allmuts.pdf')
hist(twins_dt[sum_normal_dep < 600, sum_normal_PD63383_dep], 
     xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, all mutations (362,814)', breaks = 50)
abline(v = median(twins_dt[sum_normal_dep < 600, sum_normal_PD63383_dep]), col = 'red', lwd = 2)
dev.off()

pdf('Results/20241109_p6_hist_PD62341_632muts.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all,dep_normal_PD62341] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, 632 normal')
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, dep_normal_PD62341]), col = 'blue', lwd = 2) # 264.5
dev.off()

pdf('Results/20241109_p6_hist_PD63383_632muts.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all,dep_normal_PD63383] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, 632 normal')
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, dep_normal_PD63383]), col = 'blue', lwd = 2) # 264.5
dev.off()

pdf('Results/20241109_p6_hist_PD62341_1002muts.pdf')
hist(twins_agg_vaf[,dep_normal_PD62341] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, 1002 normal')
abline(v = median(twins_agg_vaf[,dep_normal_PD62341]), col = 'purple', lwd = 2) # 264.5
dev.off()

pdf('Results/20241109_p6_hist_PD63383_1002muts.pdf')
hist(twins_agg_vaf[,dep_normal_PD63383] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, 1002 normal')
abline(v = median(twins_agg_vaf[,dep_normal_PD63383]), col = 'purple', lwd = 2) # 264.5
dev.off()

# Better visualization: position on chr vs log distance to the next mutation
# okay take each mutation and calculate the distance to the next mutation 
twins_mut_clusters = twins_agg_vaf[mut_ID %in% muts_normal_all, c('mut_ID', 'Chrom'), with=FALSE]

# plot all mutations that you identified (unless removed due to quality issues)
filters_qual = c('f1_mappedY', 'f2_FailedIndelNearby30', 'f3_lowDepthNormal', 'f3_highDepthNormal', 'f3_noReadsMapped',
                                         'f4_mtr4_presentInOne', 'f6_mtrAndVaf', 'f8_strandBiasMutOnly', 'f9_lowQualRatio')
twins_dt_filters[, sum_req_filters_qual := rowSums(.SD), .SDcols = filters_qual] 
twins_agg_allmuts = twins_dt_filters[sum_req_filters_qual==0, c('mut_ID', 'Chrom', samples_dep, samples_mtr), with=FALSE]

twins_agg_allmuts[, dep_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_agg_allmuts[, dep_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_agg_allmuts[, dep_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_agg_allmuts[, dep_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_agg_allmuts[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]
twins_agg_allmuts[, dep_all_tumour := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_agg_allmuts[, dep_all := rowSums(.SD), .SDcols = samples_dep]

twins_agg_allmuts = merge(twins_agg_allmuts, chr_lengths_dt, by = 'Chrom')
twins_agg_allmuts[, Chrom := factor(Chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                                  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                                  "chr20", "chr21", "chr22", "chrX"))]
# label mutations present in clusters 
twins_filtered_dt[, pos := tstrsplit(mut_ID, '_', fixed = TRUE, keep = 2)]
twins_filtered_dt[, pos := as.numeric(pos)]

twins_dt[, pos := tstrsplit(mut_ID, '_', fixed=TRUE, keep=2)]
twins_dt[, pos := as.numeric(pos)]

# define a function to find clusters of mutations 
find_clusters = function(numbers, range=100000, min_count = 10){
  numbers = sort(numbers)
  clusters = list()
  current_cluster = c()
  
  for (i in seq_along(numbers)){
    if (length(current_cluster)==0 || numbers[i] - tail(current_cluster, 1) <= range){
      current_cluster = c(current_cluster, numbers[i])
    } else {
      if (length(current_cluster) >= min_count){
        clusters = append(clusters, list(current_cluster))
      }
      current_cluster = c(numbers[i])
    }
  }
  if (length(current_cluster) >= min_count){
    clusters = append(clusters, list(current_cluster))
  }
  
  return(clusters %>% unlist())
}  

muts_clusters = c()

for (chr in twins_filtered_dt[, Chrom] %>% unique){
  pos_chr = sort(twins_filtered_dt[Chrom == chr, pos] %>% unlist())
  clusters_chr = find_clusters(pos_chr)
  muts_clusters_chr = twins_filtered_dt[Chrom == chr & pos > min(clusters_chr) & pos < max(clusters_chr), mut_ID] %>% unlist()
  muts_clusters = c(muts_clusters, muts_clusters_chr)
}

muts_clusters = muts_clusters %>% unlist()

twins_agg_allmuts[, pos := tstrsplit(mut_ID, '_', fixed = TRUE, keep = 2)]
twins_agg_allmuts[, pos := as.numeric(pos)]

twins_agg_allmuts[, mut_cat := factor(fcase(
  mut_ID %in% muts_clusters, 'clusters',
  !mut_ID %in% muts_clusters, 'other'))]

twins_agg_allmuts[, mut_cat2 := factor(fcase(
  mut_ID %in% muts_clusters & mut_ID %in% muts_normal_all, 'clusters, normal samples',
  mut_ID %in% muts_clusters & !mut_ID %in% muts_normal_all, 'clusters, other samples',
  !mut_ID %in% muts_clusters & mut_ID %in% muts_normal_all, 'no clusters, normal samples',
  !mut_ID %in% muts_clusters & !mut_ID %in% muts_normal_all, 'other'))]

my_colors = c('purple', 'grey')
names(my_colors) = levels(twins_agg_allmuts[, mut_cat])

twins_agg_allmuts[, mut_cat3 := factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  !mut_ID %in% muts_normal_all, 'other'))]

my_colors3 = c('#cc7f09', 'grey')
names(my_colors3) = levels(twins_agg_allmuts[, mut_cat3])

for (chr in Chrom){
  dt = twins_agg_allmuts[Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt %>% arrange(mut_cat3), aes(x = pos, y = dep_all, col = mut_cat3))+
    geom_point(size=1.5, alpha = 0.8)+
    theme_classic(base_size = 12)+
    scale_color_manual(values = c(my_colors3))+
    labs(x = 'Genomic position', y = 'Total coverage (agg all samples)', col = 'Mutation')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))
  ggsave(glue('Results/20241109_dep_across_genome_{chr}_col3.pdf'), height=3, width=6.5)
}

paste('Number of mutations in clusters:', length(muts_clusters)) # 353
paste('Number of mutations in clusters:', sum(muts_normal_all %in% muts_clusters)) # 353

######################################################################################################
# Effects of those copy number changes on protein-coding genes

gene_effects_632 = twins_dt[mut_ID %in% muts_normal_all & Gene != '-', c('mut_ID', 'Gene', 'Effect')]
genes_632 = data.table(table(gene_effects_632[, 'Gene']))
genes_632_clusters = genes_632[N > 1]
# there are some genes which would have been amplified if my hypothesis is true
# I don't see anything that screams childhood cancer though 

######################################################################################################
# How do mutations present in all normal samples look like in tumour samples?

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_tumour_all))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD62341 normal)', y = 'VAF (tumour)')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD62341_normal_vs_tumour_632.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD63383, y = vaf_tumour_all))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD63383 normal)', y = 'VAF (tumour)')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD63383_normal_vs_tumour_632.pdf'), height=5, width=7)

twins_agg_vaf[, vaf_ratio_PD62341 := vaf_normal_PD62341 / vaf_tumour_all]
twins_agg_vaf[, vaf_ratio_PD62341_log := log2(vaf_ratio_PD62341)]
twins_agg_vaf[, vaf_ratio_PD63383 := vaf_normal_PD63383 / vaf_tumour_all]
twins_agg_vaf[, vaf_ratio_PD63383_log := log2(vaf_ratio_PD63383)]

pdf('Results/20241107_p1_hist_vaf_ratio_PD62341_normal_vs_tumour_632.pdf')
hist(twins_agg[mut_ID %in% muts_normal_all, vaf_ratio_PD62341] %>% unlist(), xlab = 'Ratio VAF PD62341 normal / tumour', main = '632 mutations (all normal)')
dev.off()

pdf('Results/20241107_p1_hist_vaf_ratio_PD62341_normal_vs_tumour_632_log.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_PD62341_log] %>% unlist(), xlab = 'log2(Ratio VAF PD62341 normal / tumour)', main = '632 mutations (all normal)')
dev.off()

pdf('Results/20241107_p1_hist_vaf_ratio_PD63383_normal_vs_tumour_632.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_PD63383] %>% unlist(), xlab = 'Ratio VAF PD63383 normal / tumour', main = '632 mutations (all normal)')
dev.off()

pdf('Results/20241107_p1_hist_vaf_ratio_PD63383_normal_vs_tumour_632_log.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_PD63383_log] %>% unlist(), xlab = 'log2(Ratio VAF PD63383 normal / tumour)', main = '632 mutations (all normal)')
dev.off()

twins_agg_vaf[, mut_class_n2t_ratio_PD62341 := factor(fcase(
  vaf_ratio_PD62341_log < -0.5, 'enriched in tumour',
  vaf_ratio_PD62341_log > 0.5, 'enriched in normal (PD62341)',
  vaf_ratio_PD62341_log >= -0.5 & vaf_ratio_PD62341_log <= 0.5, 'no difference'))]

twins_agg_vaf[, mut_class_n2t_ratio_PD63383 := factor(fcase(
  vaf_ratio_PD63383_log < -0.5, 'enriched in tumour',
  vaf_ratio_PD63383_log > 0.5, 'enriched in normal (PD63383)',
  vaf_ratio_PD63383_log >= -0.5 & vaf_ratio_PD63383_log <= 0.5, 'no difference'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_tumour_all, col = mut_class_n2t_ratio_PD62341))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD62341 normal)', y = 'VAF (tumour)', col = 'Mutation category')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  scale_color_manual(values = c(col_PD62341, col_tumour, 'grey'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD62341_normal_vs_tumour_632_col.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD63383, y = vaf_tumour_all, col = mut_class_n2t_ratio_PD63383))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD63383 normal)', y = 'VAF (tumour)', col = 'Mutation category')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  scale_color_manual(values = c(col_PD63383, col_tumour, 'grey'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD63383_normal_vs_tumour_632_col.pdf'), height=5, width=7)

twins_agg_vaf[mut_ID %in% muts_normal_all & mut_class_n2t_ratio_PD62341 == 'enriched in tumour']
# chr1_16545599_G_C
# chr15_27127735_T_G
# chr4_189788936_T_C 

twins_agg_vaf[mut_ID %in% muts_normal_all & mut_class_n2t_ratio_PD63383 == 'enriched in tumour']
# chr1_16545599_G_C
# chr1_16864604_G_A
# chr1_86000_A_C

twins_agg_vaf[mut_ID %in% muts_normal_all & mut_class_n2t_ratio_PD63383 == 'enriched in normal (PD63383)' & mut_class_n2t_ratio_PD62341 == 'enriched in normal (PD62341)']
# "chr1_103717717_G_A" 
# "chr1_13225673_C_T"  
# "chr1_13225701_G_A"  
# "chr1_1646526_A_C"  
# "chr1_16652002_G_A"  
# "chr1_16772756_C_T"  
# "chr1_16889637_G_A"  
# "chr1_16892573_C_T" 
# "chr1_16905610_C_G"  
# "chr1_1697008_G_T"   
# "chr1_25409654_T_A"  
# "chr1_47072883_A_T" 
# "chr1_77402385_G_A"  
# "chr13_69946674_A_G" 
# "chr13_69946688_A_G" 
# "chr13_69946752_A_G"
# "chr18_44721988_C_T" 
# "chr18_55923631_A_G" 
# "chr18_56506391_G_C" 
# "chr9_39553748_C_T"

######################################################################################################
# Plot distribution of mutations present in all tumour samples across the genome 

for (chr in Chrom){
  dt = twins_agg[Chrom == chr & mut_ID %in% muts_tumour_all]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = tumour_vaf))+
    geom_point(size=2.5)+
    theme_classic(base_size = 12)+
    scale_color_manual(values = '#7b13d4')+
    labs(x = 'Genomic position', y = 'VAF (agg tumour)')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241109_dist_across_genome_641_tumour_{chr}.pdf'), height=3, width=4.5)
}

twins_agg[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all_nt, 'all normal, but not all tumour',
  mut_ID %in% muts_tumour_all_nt, 'all tumour, but not all normal samples',
  mut_ID %in% muts_all, 'all samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  mut_ID %in% muts_normal_only, 'only normal samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_all, muts_normal_only, muts_tumour_all_only, muts_tumour_only), 'other'))]

# make sure colors are always the same
myColors <- c('darkred', 'purple', '#07a94b', '#f3782b',  '#2fb1f3', 'grey')
names(myColors) <- levels(twins_agg[, mut_class])

for (chr in Chrom){
  dt = twins_agg[Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = tumour_vaf, col = mut_class))+
    geom_point(size=2.5)+
    theme_classic(base_size = 10)+
    scale_color_manual(values = myColors)+
    labs(x = 'Genomic position', y = 'VAF (agg tumour)', col = 'Mutation category')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241109_dist_across_genome_1002_tumour_{chr}_cat.pdf'), height=3, width=5)
}

######################################################################################################
# CLASSES OF MUTATIONS 2: present in all normal but not all tumour samples

paste('Number of mutations present in all normal but not all tumour samples:', length(muts_normal_all_nt))
twins_agg_vaf[mut_ID %in% muts_normal_all_nt]

# Are the sets for un-adjusted and adjusted VAF very different?
setdiff(muts_normal_all_nt, muts_normal_all_nt_adj) # 11
setdiff(muts_normal_all_nt_adj, muts_normal_all_nt) # 21
Reduce(intersect, list(muts_normal_all_nt, muts_normal_all_nt_adj)) # 46

# plot a heatmap just so you know what's going on
mut_normal_nt = twins_filtered_vaf_adj[mut_ID %in% muts_normal_all_nt, c('mut_ID', samples_vaf), with=FALSE]
mut_normal_nt_mat = as.matrix(mut_normal_nt[, c(samples_vaf), with=FALSE])
rownames(mut_normal_nt_mat) = mut_normal_nt[,1] %>% unlist()  
colnames(mut_normal_nt_mat) = tstrsplit(colnames(mut_normal_nt_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_normal_nt_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_all_normal_some_tumour.pdf')
pheatmap(mut_normal_nt_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="57 muts: all normal, but not all tumour", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# heatmap on VAF-adjusted values
mut_normal_nt_adj = twins_filtered_vaf_adj[mut_ID %in% muts_normal_all_nt, c('mut_ID', samples_normal_adj, samples_tumour_adj), with=FALSE]
mut_normal_nt_mat_adj = as.matrix(mut_normal_nt_adj[, c(samples_normal_adj, samples_tumour_adj), with=FALSE])
rownames(mut_normal_nt_mat_adj) = mut_normal_nt_adj[,1] %>% unlist()  
colnames(mut_normal_nt_mat_adj) = tstrsplit(colnames(mut_normal_nt_mat_adj), '_adj', fixed=TRUE, keep=1) %>% unlist()

pheatmap(mut_normal_nt_mat_adj,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="57 muts: all normal, but not all tumour", 
         legend = T,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 

# how many tumour samples is the mutation usually absent from
pdf('Results/20241109_p7_hist_muts_all_normal_some_tumour_vaf.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_normal_all_nt, sum_tumour] %>% unlist(),
     xlab = 'Number of tumour samples with mutation', main = 'Muts in all normal but not all tumour (57), VAF')
dev.off()

pdf('Results/20241109_p7_hist_muts_all_normal_some_tumour_mtr.pdf')
hist(twins_filtered_mtr[mut_ID %in% muts_normal_all_nt, sum_tumour] %>% unlist(),
     xlab = 'Number of tumour samples with mutation', main = 'Muts in all normal but not all tumour (57), MTR')
dev.off()

# what chromosomes are these mutations on? could these be lost in some tumour clones?
chr_counts_normal_nt = data.table(table(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, Chrom]))
chr_counts_normal_nt[, V1 := factor(V1, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5',
                                                        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                                        'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                                        'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                                        'chr21', 'chr22', 'chrX'))]
ggplot(chr_counts_normal_nt, aes(x = V1, y = N))+
  geom_bar(stat = 'identity')+
  theme_classic(base_size = 14)+
  labs(x = 'Chromosome', y = 'Frequency')+
  ggtitle(glue('Muts in normal but not all tumour samples (57)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('Results/20241109_p4_muts_normal_nt_dist_chr.pdf', height = 5, width = 7.5)

muts_normal_all_nt
# "chr10_79512492_C_G" # slightly odd mapping   
# "chr12_129106160_C_T" # some mapping issues 
# "chr12_55594747_G_T" # variable mapping 
# "chr16_20575297_A_G" # variable mapping   
# "chr16_88009831_A_G" # looks okay
# "chr17_18843733_C_T" # some mapping issues 
# "chr22_24001246_G_A" # ay reads in 1 direction map terribly wrong
# "chr2_1440752_T_C" # looks great and such a fun region!!!   
# "chr2_20583359_T_C" # lots of insertions so cannot trust it  
# "chr3_190409217_C_T" # very interesting mutation pattern 
# "chr5_784531_G_A" # issues with read pair mates mapping     
# "chr6_43580159_A_G" # some mapping issues   
# "chr6_94069453_A_T" # insertions so am thinking not great? 
# "chr7_126168639_G_A" # mapping is a bit funny   
# "chr7_143742625_A_G" # mapping is a bit funny but kind of okay? 
# "chr9_133062780_G_C" # mapping issues 
# "chr9_64775433_A_G" # mapping issues yet again  
# "chr9_77762599_A_T" # mapping issues   
# "chrX_140138202_C_T" # extreme strand bias in wt reads  
# "chrX_50060796_A_G" # also poor mapping   
# "chrX_65021869_T_A" # issues with mate mapping 
# "chrX_66502763_G_A" # looks great

# chromosome 18
# "chr18_44721988_C_T" - looks great  
# "chr18_55923631_A_G" - can we please have a look at this one? 
# "chr18_56506391_G_C" - can we please have a look at this one? 

# chromosome 13
# "chr13_38383628_T_C" - looks fine 
# "chr13_69946674_A_G" - segments not correctly aligned  
# "chr13_69946688_A_G" - segments not correctly aligned  
# "chr13_69946752_A_G" - segments not correctly aligned  

# these mutations likely represent cases where a segment of the chr with the mutation was lost from some tumour samples
# but then, those mutations are not germline, so are they informative for the normal phylogeny?
pdf('Results/20241109_p6_hist_57muts_vaf.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_all], xlim = c(0, 0.6),
     xlab = 'VAF (agg normal)', main = '57 mutations (all normal, some tumour)', breaks = 20)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_all]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p6_hist_57muts_vaf_PD62341.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD62341], xlim = c(0, 0.6),
     xlab = 'VAF (agg normal)', main = '57 mutations (all normal, some tumour), PD62341', breaks = 20)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD62341]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p6_hist_57muts_vaf_PD63383.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD63383], xlim = c(0, 0.6),
     xlab = 'VAF (agg normal)', main = '57 mutations (all normal, some tumour), PD63383', breaks = 20)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD63383]), col = 'purple', lwd = 2.5)
dev.off()

# compare mutations between PD62341 and PD63383
twins_agg_vaf[, ratio_PD62341_PD63383_normal := vaf_normal_PD62341 / vaf_normal_PD63383]
twins_agg_vaf[, ratio_PD62341_PD63383_normal_log := log2(ratio_PD62341_PD63383_normal)]
twins_agg_vaf[, mut_class_ratio_n2n := factor(fcase(
  vaf_normal_PD62341 > 0.35 & vaf_normal_PD63383 > 0.35, 'putative germline',
  ratio_PD62341_PD63383_normal_log < -0.5, 'enriched in PD63383',
  ratio_PD62341_PD63383_normal_log > 0.5, 'enriched in PD62341',
  ratio_PD62341_PD63383_normal_log > -0.5 & ratio_PD62341_PD63383_normal_log < 0.5, 'no difference'))]

myColors <- c(col_PD62341, col_PD63383, 'grey', 'darkred')
names(myColors) <- levels(twins_agg_vaf[, mut_class_ratio_n2n])

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all_nt], aes(x = vaf_normal_PD62341, y = vaf_normal_PD63383, col = mut_class_ratio_n2n))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (agg PD62341 normal)', y = 'VAF (agg PD63383 normal)', col = 'Mutation category')+
  ggtitle(glue('57 mutations (all normal, not all tumour)'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.6))+
  ylim(c(0, 0.6))+
  scale_color_manual(values = myColors)
ggsave('Results/20241109_vaf_PD62341_vs_PD63383_normal_57muts.pdf', height = 5, width = 7.5)

######################################################################################################
# CLASSES OF MUTATIONS 3: MUTATIONS PRESENT IN ONLY NORMAL SAMPLES 

# examine the smaller sets of mutations
muts_normal_only_adj # 12
# "chr11_34011887_C_T" # looks real
# "chr13_50815806_A_G" # looks real
# "chr1_103587565_A_C" # looks good
# "chr21_40193588_G_A" # looks real
# "chr3_165901319_C_A" # good enough
# "chr3_77633967_C_T" # looks okay
# "chr4_15905566_C_T" # looks real
# "chr4_74625500_G_T" # looks real 
# "chr4_75704880_G_A" # looks good
# "chr6_159851462_T_C" # absolutely not real
# "chr7_149688370_C_T" # looks okay
# "chrX_115066661_C_T" # looks real

muts_normal_only
# "chr11_34011887_C_T" # looks real
# "chr13_50815806_A_G" # looks real
# "chr1_103587565_A_C" # looks good
# "chr21_40193588_G_A" # looks real
# "chr4_15905566_C_T" # looks real
# "chr4_74625500_G_T" # looks real 
# "chr4_75704880_G_A" # looks good
# "chr6_159851462_T_C" # absolutely not real
# "chr7_149688370_C_T" # looks okay
# "chrX_115066661_C_T" # looks real

setdiff(muts_normal_only, muts_normal_only_adj)
# found only w/o adjusting for purity
# "chr17_20476029_C_T" # maps to several places
# "chr21_28681933_A_C" # not excellent
# "chr7_110349267_G_A" # bit questionable (maps to several places)
# "chr7_73831920_C_T" # looks okay

# Can I look at the VAF of each of those mutations
twins_vaf_adj_sub = twins_filtered_vaf_adj[, c('mut_ID', samples_normal_adj, samples_tumour_vaf), with=FALSE]
twins_vaf_adj_melt = melt(twins_vaf_adj_sub, id.vars = 'mut_ID')
twins_vaf_adj_melt[, sample := tstrsplit(variable, '_', keep = 1, fixed = TRUE)]
twins_vaf_adj_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_adj_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_adj_melt[, sample_type := paste(status, twin, sep = ', ')]
twins_vaf_adj_melt[, sample_type1 := factor(fcase(
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

ggplot(twins_vaf_adj_melt[mut_ID %in% muts_normal_only_adj], aes(x = value, y = mut_ID, col = sample_type1))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only normal samples, adj VAF): 12'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples_adj.pdf', height = 5, width = 7.5)

ggplot(twins_vaf_adj_melt[mut_ID %in% muts_normal_only], aes(x = value, y = mut_ID, col = sample_type1))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only normal samples): 16'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples.pdf', height = 5, width = 7.5)

# Use non-adjusted VAF values
twins_vaf_sub = twins_filtered_vaf[, c('mut_ID', samples_vaf), with=FALSE]
twins_vaf_melt = melt(twins_vaf_sub, id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_', keep = 1, fixed = TRUE)]
twins_vaf_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_melt[, sample_type := paste(status, twin, sep = ', ')]
twins_vaf_melt[, sample_type1 := factor(fcase(
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

ggplot(twins_vaf_melt[mut_ID %in% muts_normal_only], aes(x = value, y = mut_ID, col = sample_type1))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('16 mutations present only in normal samples'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples_nonadj_16.pdf', height = 5, width = 7.5)

muts_normal_only_val = c("chr11_34011887_C_T", "chr13_50815806_A_G", "chr1_103587565_A_C", "chr21_40193588_G_A",
                         "chr3_165901319_C_A", "chr3_77633967_C_T", "chr4_15905566_C_T", "chr4_74625500_G_T", 
                         "chr4_75704880_G_A", "chr7_149688370_C_T", "chrX_115066661_C_T", "chr7_73831920_C_T")
ggplot(twins_vaf_melt[mut_ID %in% muts_normal_only_val], aes(x = value, y = mut_ID, col = sample_type1))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('12 mutations only in normal samples (verified)'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples_nonadj_12.pdf', height = 5, width = 7.5)

twins_vaf_melt[, sample_type2 := factor(fcase(
  sample == 'PD62341v', 'PD62341, normal: spleen',
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]
ggplot(twins_vaf_melt[mut_ID %in% muts_normal_only_val], aes(x = value, y = mut_ID, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('12 mutations only in normal samples (verified)'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, 'darkgoldenrod', col_PD63383, col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples_nonadj_12_spleen.pdf', height = 5, width = 7.5)

twins_vaf_melt[, sample_type3 := factor(fcase(
  sample == 'PD62341v', 'PD62341, normal: spleen',
  sample == 'PD63383bb', 'PD63383, normal: skin',
  sample == 'PD63383w', 'PD63383, normal: spleen',
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

ggplot(twins_vaf_melt[mut_ID %in% muts_normal_only_val] %>% arrange(sample_type3), aes(x = value, y = mut_ID, col = sample_type3))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('12 mutations only in normal samples (verified)'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, 'darkgoldenrod', col_PD63383, '#0d7dc1', '#efb669', col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples_nonadj_12_spleen_skin.pdf', height = 5, width = 7.5)

######################################################################################################
# CLASSES OF MUTATIONS 4: MUTATIONS PRESENT IN ALL TUMOUR BUT ONLY SOME NORMAL SAMPLES

# lists of relevant mutations 
muts_tumour_all_nt
muts_tumour_all_nt_adj

# how many normal samples are these usually present in?
pdf('Results/20241109_p7_hist_nr_normal_66muts_tumour_nt.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt, sum_normal],
     xlab = 'Number of normal samples', main = '66 mutations (all tumour, not all normal)')
dev.off()

# heatmap to have a look at how this looks like generally
mut_tumour_nt = twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt, c('mut_ID', samples_vaf), with=FALSE]
mut_tumour_nt_mat = as.matrix(mut_tumour_nt[, c(samples_vaf), with=FALSE])
rownames(mut_tumour_nt_mat) = mut_tumour_nt[,1] %>% unlist()  
colnames(mut_tumour_nt_mat) = tstrsplit(colnames(mut_tumour_nt_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_tumour_nt_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_all_tumour_some_normal.pdf')
pheatmap(mut_tumour_nt_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="66 muts: all tumour, but not all normal", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

pdf('Results/20241109_p6_heatmap_all_tumour_some_normal_rownames.pdf', height = 80)
pheatmap(mut_tumour_nt_mat,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="66 muts: all tumour, but not all normal", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# heatmap for adjusted vaf values
mut_tumour_nt_adj = twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj, c('mut_ID', samples_normal_adj, samples_tumour_adj), with=FALSE]
mut_tumour_nt_mat_adj = as.matrix(mut_tumour_nt_adj[, c(samples_normal_adj, samples_tumour_adj), with=FALSE])
rownames(mut_tumour_nt_mat_adj) = mut_tumour_nt_adj[,1] %>% unlist()  
colnames(mut_tumour_nt_mat_adj) = tstrsplit(colnames(mut_tumour_nt_mat_adj), '_adj', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_tumour_nt_mat_adj)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_all_tumour_some_normal_vaf_adj.pdf')
pheatmap(mut_tumour_nt_mat_adj,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="66 muts: all tumour, but not all normal, adj VAF", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# are there any mutations present in PD63383 normal samples?
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 0] # 11 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 1] # 21
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 2] # 1

# are there any mutations that are only shared with normal PD63383?
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD62341 == 0] # 0

# get out mutations absent from PD63383 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1] # 33 

# which normal PD62341 samples are these shared with?
pdf('Results/20241109_p7_hist_nr_normal_66muts_tumour_nt_PD62341.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1, sum_normal_PD62341],
     xlab = 'Number of normal PD62341 samples with mutation', main = '66 mutations (all tumour, not all normal)') # 33 
dev.off()

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 2] # 8 
# likely all contamination
# chr10_100754461_C_T aa, h
# chr14_104500332_C_T aa, h
# chr15_23462705_C_A aa, h
# chr17_40061856_G_C aa, h
# chr4_147908994_T_C aa, h
# chr6_95827754_C_A aa, h
# chrX_66719643_C_G aa, h
# chrX_68803487_C_T aa, h
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 3] # 1 
# chr5_44907911_G_A q, aa, h (adjusted VAF in q = 0.086, so maybe not in the end)
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 4] # 2 
# chr3_62055057_C_G q, aa, h, n
# chr3_62055077_G_C q, aa, h, n 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 5] # 3 
# chr16_5479739_C_T q, aa, ad, h, n
# chr2_95662131_G_A q, aa, ad, h, n
# chr3_50106043_C_T q, aa, ad, h, n
# absent from v - but it is likely that MTR is lower due to twin-twin contamination 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 6] # 2
# chr14_105458006_C_A # all PD62341
# chr17_33422229_C_A # all PD62341

twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 1] # 13
# only sharing with aa or h
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 4] # 2
# chr3_62055057_C_G q, aa, h, n 
# chr3_62055077_G_C q, aa, h, n 
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 5] # 3
# chr16_5479739_C_T q, aa, ad, h, n
# chr2_95662131_G_A q, aa, ad, h, n
# chr3_50106043_C_T q, aa, ad, h, n
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 6] # 2
# chr14_105458006_C_A # all PD62341
# chr17_33422229_C_A # all PD62341

twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 1] # 13
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 4] # 2
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 5] # 3
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 6] # 2

twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 4] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 5] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 6] # 1
# chr15_49480646_T_A (bb and t - why t?)

# plot VAF 

######################################################################################################
# MUTATION CLASSES 5: MUTATIONS IN ONLY TUMOUR SAMPLES

ggplot(twins_vaf_adj_melt[mut_ID %in% muts_tumour_only], aes(x = value, y = mut_ID, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only tumour samples): 31'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour_PD62341, col_tumour_PD63383))
ggsave('Results/20241109_p4_muts_only_tumour_samples_vaf_adj.pdf', height = 6.5, width = 7.5)

twins_vaf_melt = melt(twins_filtered_vaf[, c('mut_ID', samples_vaf), with=FALSE], id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_vaf_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_melt[, sample_type := paste(status, twin, sep = ', ')]
twins_vaf_melt[, sample_type1 := factor(fcase(
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

ggplot(twins_vaf_melt[mut_ID %in% muts_tumour_only], aes(x = value, y = mut_ID, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only tumour samples): 31'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour_PD62341, col_tumour_PD63383))
ggsave('Results/20241109_p4_muts_only_tumour_samples_vaf.pdf', height = 6.5, width = 7.5)

######################################################################################################
# MUTATION CLASSES 4B: MUTATIONS ONLY IN PD63383 TUMOUR

# what about mutations present in PD63383 only?
twins_filtered_vaf[sum_tumour_PD63383 > 0 & sum_tumour_PD62341 == 0] 
twins_filtered_vaf_adj[sum_tumour_PD63383 > 0 & sum_tumour_PD62341 == 0] 
# everything looks great unless stated otherwise 
# the same set of mutations comes up in adj and non-adj VAF
# "chr10_100308403_C_T" 
# "chr10_31333522_G_A"  
# "chr11_80115755_C_T"  
# "chr12_43132316_G_A" 
# "chr12_80762031_G_A"  
# "chr19_43348626_C_A" # not as good but probably fine  
# "chr1_160088990_G_A"  
# "chr1_51714096_G_C"  
# "chr22_25281865_C_T" 
# "chr22_43290224_C_T" # double check
# "chr2_72939205_C_T"   
# "chr2_82147172_T_C"   
# "chr3_137508691_C_A" 
# "chr4_179587218_G_A"  
# "chr5_157248612_A_G" # present in only 1 sample (VAF = 0.09, but 4 reads)
# "chr5_28014472_C_T"   
# "chr5_54017866_G_T" 
# "chr6_86549064_C_T"  
# "chr6_92179408_A_G"   
# "chr7_13677323_A_G"   
# "chr7_148446413_G_A"  
# "chr7_49732059_G_A"  
# "chr9_98275090_C_T"   
# "chrX_117418608_G_A"  
# "chrX_119915655_G_A"  
# "chrX_36779694_C_T"  
# "chrX_48453603_G_A"   
# "chrX_70352624_G_A"   
# "chrX_72671786_A_T" # looks plausible but really maps to EVERYWHERE on Blat    
# "chrX_77681162_G_A"

muts_PD63383_tumour_only = twins_filtered_vaf_adj[sum_tumour_PD63383 > 0 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist() 

# in PD63383 only but not tumour only 
# these are due to persisting contamination in the skin (even after adjusting for contamination)
setdiff(muts_PD63383_tumour_only, muts_tumour_only) # 18

# in tumour only but not only in PD63383
# are these mutations only present in PD62341 tumour?
setdiff(muts_tumour_only, muts_PD63383_tumour_only) # 19
twins_filtered_vaf[mut_ID %in% setdiff(muts_tumour_only, muts_PD63383_tumour_only), 
         c('mut_ID', 'sum_tumour', 'sum_tumour_PD62341')]
# chr8_21444725_G_A # present in both PD63383 tumours and 2 PD62341 tumour samples
# chr7_152562130_G_A # present in both PD63383 tumours and 5 PD62341 tumour samples

twins_filtered_vaf[mut_ID %in% setdiff(muts_tumour_only, muts_PD63383_tumour_only), c('mut_ID', samples_tumour_PD62341_vaf), with=FALSE]
# note that this is with un-adjusted values and some samples we estimate to be of very low purity
# "chr10_72139813_C_T" ak, ap // ak, ap
# "chr11_43029745_A_T" ae, aj, b // ae, aj, b, u
# "chr12_19555916_C_T" ae, aj // ae, aj, u
# "chr13_28953051_G_T" ae, ag, aj // ae, ag, aj, b, u
# "chr13_62110424_G_A" am, ap // am, ap
# "chr15_68707110_G_A" ak, ap // ak, ap
# "chr16_17805832_C_T" ae, ag, aj, b, u // ae, ag, aj, b, u
# "chr1_14644486_C_A" ae, aj // ae, aj, b
# "chr1_38311094_G_T" ak, ap // ak, ap
# "chr2_231591536_G_T" ae, am // ae, am, u
# "chr3_189240074_C_A" am // ak, am, ap
# "chr3_95470911_T_A" ak // ak, ap
# "chr5_58182780_T_A" ak, ap // ak, ap
# "chr7_152562130_G_A" ae, aj, am, ap, b // ae, aj, am, ak, ap, b
# "chr7_153784080_C_G" ae, ag, am, b // ae, ag, ak, am, ap, b
# "chr8_21444725_G_A" am, b // am, b, u
# "chr8_91252357_A_C" ak // ak, ap
# "chrX_64561555_C_T" ak, ap // ag, ak, ap
# "chrX_9804641_G_A" ak, ap // ak, ap

# can I get a heatmap of this?
PD62341_subclonal = twins_filtered_vaf_adj[mut_ID %in% setdiff(muts_tumour_only, muts_PD63383_tumour_only), c('mut_ID', samples_tumour_vaf), with=FALSE]
PD62341_subclonal_mat = as.matrix(PD62341_subclonal[, c(samples_tumour_vaf), with=FALSE])
rownames(PD62341_subclonal_mat) = PD62341_subclonal[,1] %>% unlist()  
colnames(PD62341_subclonal_mat) = tstrsplit(colnames(PD62341_subclonal_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Twin = c(rep('PD62341', 6), rep('PD63383', 2), rep('PD62341', 2)))
rownames(col_annotation) = colnames(PD62341_subclonal_mat)
annotation_colors = list(Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# annotation data frame 
pdf('Results/20241109_p6_heatmap_PD62341_subclonal.pdf')
pheatmap(PD62341_subclonal_mat,
         cellwidth=10, cellheight=5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="19 PD62341 subclonal mutations", 
         legend = T,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Analysis of tumour evolution: coverage vs VAF plot

twins_agg_vaf[, mut_class := factor(fcase(
  mut_ID %in% muts_all, 'all samples',
  mut_ID %in% muts_normal_all_nt, 'all normal samples, not all tumour',
  mut_ID %in% muts_tumour_all_nt, 'all tumour samples, not all normal',
  mut_ID %in% muts_normal_only, 'normal only',
  mut_ID %in% muts_tumour_only, 'tumour only',
  !mut_ID %in% c(muts_tumour_only, muts_normal_only, muts_all, 
                 muts_normal_all_nt, muts_tumour_all_nt), 'other')) ]
twins_agg_vaf[, mut_class := factor(mut_class, levels = c(
  'all samples', 'all normal samples, not all tumour', 'all tumour samples, not all normal',
  'normal only', 'tumour only', 'other'))]

myColors <- c('lightblue', 'purple', '#07a94b', '#f3782b',  'darkred', 'grey')
names(myColors) <- levels(twins_agg_vaf[, mut_class])

# Plot coverage vs VAF 
ggplot(twins_agg_vaf %>% arrange(mut_class), aes(x = dep_tumour_all, y = vaf_tumour_all, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  scale_color_manual(values = myColors)+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002.pdf'), height = 5, width = 6)

ggplot(twins_agg_vaf %>% arrange(mut_class), aes(x = dep_normal_all, y = vaf_normal_all, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all normal samples)', y = 'VAF (all normal samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across normal samples'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  scale_color_manual(values = myColors)
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002_normal.pdf'), height = 5, width = 6)

# color by type of mutation
twins_agg_vaf[, c('Ref', 'Alt') := tstrsplit(mut_ID, '_', fixed=T, keep=3:4)]
twins_agg_vaf[, mut_type := paste0(Ref, '>', Alt)]
twins_agg_vaf[, mut_type := as.factor(fcase(
  mut_type == 'A>T', 'T>A',
  mut_type == 'A>C', 'T>G',
  mut_type == 'G>C', 'C>G',
  mut_type == 'G>T', 'C>A',
  mut_type == 'A>G', 'T>C',
  mut_type == 'G>A', 'C>T',
  mut_type == 'C>A', 'C>A',
  mut_type == 'C>G', 'C>G',
  mut_type == 'C>T', 'C>T',
  mut_type == 'T>A', 'T>A',
  mut_type == 'T>C', 'T>C',
  mut_type == 'T>G', 'T>G'))]

ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = mut_type))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_mut_type.pdf'))

# color by chromosome 
twins_agg_vaf[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = Chrom))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_chrom.pdf'))

twins_agg_vaf[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'))]

ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = loss))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_loss.pdf'))

# add number of tumour samples the mutation is present in 
twins_agg_vaf = merge(twins_agg_vaf, twins_filtered_vaf[, c('mut_ID', 'sum_tumour'), with=FALSE], by = 'mut_ID')
twins_agg_vaf[, sum_tumour := as.factor(sum_tumour)]
ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = sum_tumour))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = '# tumour samples\nmut detected in')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_sum_tumour.pdf'))

######################################################################################################
# Plot VAF and coverage across the genome 
twins_agg_vaf[, pos := tstrsplit(mut_ID, '_', fixed=T, keep=2)]
twins_agg_vaf[, pos := as.numeric(pos)]
twins_agg_vaf[, Chrom_nr := tstrsplit(Chrom, 'chr', fixed=T, keep=2)]
twins_agg_vaf[, Chrom_nr := factor(Chrom_nr, levels = 
                                        c('1', '2', '3', '4', '5',
                                          '6', '7', '8', '9', '10',
                                          '11', '12', '13', '14', '15',
                                          '16', '17', '18', '19', '20',
                                          '21', '22', 'X'))]

ggplot(twins_agg_vaf, aes(x = pos, y = vaf_tumour_all, col=mut_class))+
  geom_point(size=1.5)+
  labs(x = 'Position', y = 'VAF (all tumour samples)', col='Mutation class')+
  ggtitle(glue('VAF across the genome'))+
  facet_grid(~Chrom_nr)+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  #  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))
ggsave(glue('Results/20241109_vaf_across_the_genome_1002.pdf'))

######################################################################################################
# Get out mutations that are well mapped in all tumour samples 

for(sample in samples_tumour){
  tumour_sample = twins_filtered_dt[, c('mut_ID', glue('{sample}_MTR'), glue('{sample}_VAF'))]
  setnames(tumour_sample, c(glue('{sample}_MTR'), glue('{sample}_VAF')), c('mtr', 'vaf'))
  tumour_sample[, glue('mapped_{sample}') :=  as.numeric(mtr >= 4) * as.numeric(vaf >= 0.1)]
  twins_filtered_dt = merge(twins_filtered_dt, tumour_sample[, c('mut_ID', glue('mapped_{sample}')), with = FALSE], by = 'mut_ID')
}

mapped_samples = paste0('mapped_', samples_tumour)
twins_filtered_dt[, tumour_mapped_sum := rowSums(.SD == 1), .SDcols = mapped_samples]
twins_filtered_dt[, tumour_mapped_all := as.numeric(tumour_mapped_sum == 10)] 

# identify mutations that are well mapped in all tumour samples
muts_well_mapped_tumour = twins_filtered_dt[tumour_mapped_all==1, mut_ID] %>% unlist() # 641

# require min VAF 0.2 across the tumour 
muts_mapped_tumour = twins_tumour_agg[tumour_vaf > 0.2 & mut_ID %in% muts_well_mapped_tumour, mut_ID] %>% unlist() # 630 
paste('Number of mutations in the tumour:', length(muts_mapped_tumour)) # 630 

# plot the distribution of VAF of those mutations across the tumour 
pdf('Results/20241109_p3_hist_vaf_dist_630muts.pdf')
hist(twins_tumour_agg[mut_ID %in% muts_mapped_tumour, tumour_vaf], breaks = 20,
     xlab = 'VAF (total tumour)', main = '630 mutations (mapped all tumour samples)', xlim = c(0, 1))
abline(v = median(twins_tumour_agg[mut_ID %in% muts_mapped_tumour, tumour_vaf]), col = 'purple')
dev.off()

######################################################################################################
# Patterns of sharing with normal samples 

ggplot(twins_tumour_agg[mut_ID %in% muts_mapped_tumour], aes(x = sum_tumour_dep, y = tumour_vaf, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across the tumour: mutations mapped in all tumour samples'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_630mapped.pdf'))

# Aggregate data for all tumour samples together 
twins_filtered_mtr[, sum_PD62341_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dep[, sum_PD62341_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_filtered_mtr[, sum_PD63383_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_filtered_dep[, sum_PD63383_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]

twins_PD62341_agg = merge(twins_filtered_mtr[, c('mut_ID', 'sum_PD62341_mtr'), with=FALSE], 
                         twins_filtered_dep[, c('mut_ID', 'sum_PD62341_dep'), with=FALSE], by = 'mut_ID')
twins_PD62341_agg[, PD62341_normal_vaf := sum_PD62341_mtr / sum_PD62341_dep]

twins_PD63383_agg = merge(twins_filtered_mtr[, c('mut_ID', 'sum_PD63383_mtr'), with=FALSE], 
                          twins_filtered_dep[, c('mut_ID', 'sum_PD63383_dep'), with=FALSE], by = 'mut_ID')
twins_PD63383_agg[, PD63383_normal_vaf := sum_PD63383_mtr / sum_PD63383_dep]

# Data table with all aggregated data 
twins_agg = merge(merge(twins_tumour_agg, twins_PD62341_agg, by = 'mut_ID'), twins_PD63383_agg, by = 'mut_ID')
twins_agg_filt = twins_agg[mut_ID %in% muts_mapped_tumour]

ggplot(twins_agg_filt, aes(x=tumour_vaf, y=PD62341_vaf))+
  geom_point()+
  theme_classic(base_size=15)+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  labs(x = 'VAF (total tumour)', y = glue('VAF in PD62341'), col = 'Mutation category')+
  ggtitle(glue('630 mutations mapped in the tumour, PD62341'))
ggsave(glue('Results/20241109_vaf_tumour_vs_PD62341_630mapped.pdf'))

ggplot(twins_agg_filt, aes(x=tumour_vaf, y=PD63383_vaf))+
  geom_point()+
  theme_classic(base_size=15)+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  labs(x = 'VAF (total tumour)', y = glue('VAF in PD63383'), col = 'Mutation category')+
  ggtitle(glue('630 mutations mapped in the tumour, PD63383'))
ggsave(glue('Results/20241109_vaf_tumour_vs_PD63383_630mapped.pdf'))

# how many of those mutations are also present in all normal samples?
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD63383 == 6] # 594
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 6] # 593

twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 6] # 588

twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 6] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 4] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 5] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 3] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 2] # 0

twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 1 
# chr14_105458006_C_A our favourite one!
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 2 
# chr3_62055057_C_G aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad
# chr3_62055077_G_C aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 0] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 2 (h and aa only)
# chr10_100754461_C_T # looks okay
# chr6_95827754_C_A not the best mapping quality but fine, aa and h (possibly contamination) 
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 6 (h only)
# chr13_60216504_A_C # looks okay
# chr16_32621521_C_A # mapping quality not excellent (double check because it doesn't look that bad)
# chr17_78256883_C_T # looks okay
# chr22_17774668_G_A # looks okay 
# chr22_30553442_G_T # looks okay
# chr2_199565129_G_A # looks okay

######################################################################################################
# Looking at tumour-clonal mutations

twins_filtered_vaf[sum_tumour > 1 & sum_normal == 0] # 28
twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj == 0] # 52

twins_filtered_vaf[sum_tumour > 1 & sum_normal >= 1 & sum_normal <= 3] # 117
twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj >= 1 & sum_normal_adj <= 3] # 106

# which normal samples are these mutations present in?
colSums(twins_filtered_vaf[sum_tumour > 1 & sum_normal >= 1 & sum_normal <= 3, c(samples_normal_vaf), with=FALSE] >= 0.1)
# PD62341aa = 27
# PD62341h = 75
# PD63383bb = 86

colSums(twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj >= 1 & sum_normal_adj <= 3, c(samples_normal_adj), with=FALSE] >= 0.1)
# PD62341aa = 12
# PD62341h = 52
# PD63383bb = 54

# identify mutations which are present in normal samples other than contaminated ones
contaminated_samples = c('PD62341aa', 'PD62341h', 'PD63383bb')
contaminated_samples_vaf = paste0(contaminated_samples, '_VAF')
contaminated_samples_adj = paste0(contaminated_samples, '_adj')

twins_filtered_vaf[, sum_normal_cont := rowSums(.SD >= 0.1), .SDcols = contaminated_samples_vaf]
twins_filtered_vaf[sum_tumour > 1 & sum_normal <= 3 & sum_normal == sum_normal_cont] # 124

twins_filtered_vaf_adj[, sum_normal_cont_adj := rowSums(.SD >= 0.1), .SDcols = contaminated_samples_adj]
twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj <= 3 & sum_normal_adj == sum_normal_cont_adj] # 136

# check the other mutations (present in normal samples but not contaminated ones)
twins_filtered_vaf[sum_tumour > 1 & sum_normal <= 3 & sum_normal != sum_normal_cont, mut_ID] 
# "chr10_49977949_A_G" # poor mapping 
# "chr13_18565005_G_A" # poor mapping 
# "chr15_30811818_A_G" # poor mapping 
# "chr15_56691722_C_T" # looks real q, bb
# "chr18_71485446_G_A" # looks real q, bb, h
# "chr19_43132229_T_C" # looks real aa, PD63383u ???
# "chr19_7878006_A_G" # not real  
# "chr1_83136026_A_T" # poor mapping 
# "chr20_44114996_C_T" # looks real n, h (n after adjustment)
# "chr2_57814739_A_C" # looks real q, h (q after adjustment)
# "chr5_70784531_C_A" # poor mapping 
# "chr7_120677593_C_T" # looks real n, h, bb (n, h after adjustment)
# "chr7_139659050_G_A" # looks real q, h
# "chr8_131571989_A_T" # looks real q, h, bb (q after adjustment)
# "chr8_46678612_T_C" # looks real q, h, n (q, n after adjustment)
# "chr9_100061581_C_T" # double check, not convinced PD63383ak, w
# "chr9_100061582_G_A" # double check, not convinced PD63383ak, w
# "chr9_12824689_G_A" # looks good q, bb, h
# "chr9_41808224_G_T" # poor mapping  
# "chrX_124531709_C_T" # looks good q, h, bb

muts_tumour_spec = twins_filtered_vaf[sum_tumour > 1 & sum_normal <= 3, mut_ID] %>% unlist() # 144 (NB this is just VAF-based)

# how many tumour samples are these mutations usually in?
pdf('Results/20241109_p4_hist_nr_tumour_samples_144muts.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_spec, sum_tumour], 
     xlab = 'Number of tumour samples with mutation', main = '144 mutations (min 2 tumour, max 3 normal)')
dev.off()

pdf('Results/20241109_p4_hist_vaf_144muts.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_tumour_spec, vaf_tumour_all], breaks = 20,
     xlab = 'VAF (agg tumour)', main = '144 mutations (min 2 tumour, max 3 normal)')
dev.off()

# heatmap so it is easier to think about it 
mut_tumour_144 = twins_agg_vaf[mut_ID %in% muts_tumour_spec, c('mut_ID', samples_vaf), with=FALSE]
mut_tumour_144_mat = as.matrix(mut_tumour_144[, c(samples_vaf), with=FALSE])
rownames(mut_tumour_144_mat) = mut_tumour_144[,1] %>% unlist()  
colnames(mut_tumour_144_mat) = tstrsplit(colnames(mut_tumour_144_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_tumour_144_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_tumour_144muts.pdf')
pheatmap(mut_tumour_144_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="144 mutations: likely tumour samples only", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# more investigation of PD63383 mutations
twins_vaf_sub = twins_agg_vaf[, c('mut_ID', samples_vaf), with=FALSE]
twins_vaf_melt = melt(twins_vaf_sub, id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_vaf_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'
))]
twins_vaf_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'
))]
twins_vaf_melt[, sample_type := paste(status, twin, sep = ', ')]
                 
ggplot(twins_vaf_melt[mut_ID %in% muts_PD63383_tumour_only], aes(x = value, y = mut_ID, col = sample_type))+
  geom_point()+
  theme_classic(base_size=14)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour, col_tumour_PD63383))+
  labs(x = 'VAF', y = glue('Mutation'), col = 'Sample type')+
  ggtitle(glue('PD63383 tumour-restricted mutations'))
ggsave(glue('Results/20241109_vaf_PD63383restircted_tumour.pdf'), height = 8, width = 6)

# investigation of other classes of mutations
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2 & sum_tumour_PD63383 == 2] # 28

# mutations present in 2 tumour samples not in PD63383 tumour 
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2 & sum_tumour_PD63383 != 2] # 13
# chr10_72139813_C_T # ak, ap looks okay
# chr11_43029745_A_T # ae, b
# chr12_19555916_C_T # ae, aj
# chr13_18565005_G_A # PD63383aq, PD62341ap
# chr13_62110424_G_A # am, ap
# chr15_30811818_A_G # ag, ap
# chr15_68707110_G_A # ak, ap
# chr1_14644486_C_A # ae, aj
# chr1_38311094_G_T # ak, ap
# chr2_231591536_G_T # ae, aj
# chr5_58182780_T_A # ak, ap
# chrX_64561555_C_T # ak, ap
# chrX_9804641_G_A # ak, ap

twins_filtered_vaf[sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0] # 30
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0] # 28
setdiff(twins_filtered_vaf[sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist(), twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist())
# "chr22_43290224_C_T" - not buying it 
# "chr5_157248612_A_G"

twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 10] # 27
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 9] # 38
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 8] # 8
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 7] # 1
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 6] # 8
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 5] # 8
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 4] # 10
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 3] # 3
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 2] # 41

# add a column for names of samples which the mutation is missing in:
twins_filtered_vaf[, missing_tumour_sample := apply(twins_filtered_vaf[, ..samples_tumour_vaf], 1, function(row) {
  cols_vaf01 = samples_tumour_vaf[which(row < 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

twins_filtered_vaf[, present_tumour_sample := apply(twins_filtered_vaf[, ..samples_tumour_vaf], 1, function(row) {
  cols_vaf01 = samples_tumour_vaf[which(row >= 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 9, missing_tumour_sample])
# 1 missing in b, 2 in ag, 35 in ak

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 8, missing_tumour_sample])
# 1 missing in ag and u, 2 missing in ak and u, 5 missing in ak and ag

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 7, missing_tumour_sample])
# missing in ak, ag, u

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 6, missing_tumour_sample])
# 1 missing in ak, ag, ae, aj
# 3 missing in ak, ag, u, ae
# 2 missing in ak, ag, u, aj
# 1 missing in aj, ap, b, PD63383ap # chr1 - loss? (mapping isn't great so not 100% convinced)
# 1 missing in am, ap, b, PD63383ap # chr19 - loss?

# mutations in 5 tumour samples or less
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 5]
# chr10_49977949_A_G # not very good 
# chr11_26060345_G_T # looks okay
# chr12_112652964_C_T # looks okay
# chr13_103318073_T_A # looks okay
# chr16_17805832_C_T # looks okay
# chr17_70958446_C_T # looks okay
# chr4_135461656_G_A # looks okay
# chr6_8310684_G_T # looks okay
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 4]
# chr12_114403258_G_A # looks okay
# chr12_82832044_C_A # looks okay
# chr18_1379467_G_T # looks okay
# chr19_7878006_A_G # looks horrible
# chr5_70784531_C_A # terrible mapping 
# chr6_168731130_C_T # looks okay
# chr7_153784080_C_G # only on poorly mapped reads
# chr8_21444725_G_A # looks okay
# chr9_100061581_C_T # looks okay but a bit funny 
# chr9_41808224_G_A # poor mapping 
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 3]
# chr13_28953051_G_T # looks okay
# chr8_141860196_G_A # looks okay
# chr9_100061582_G_A # looks okay but a bit funny 
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2] # checked these and generally look fine

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 5, present_tumour_sample])
table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 4, present_tumour_sample])
table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 3, present_tumour_sample])
table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2, present_tumour_sample])

######################################################################################################
# what about mutations present in all tumour but not all normal samples?
muts_tumour_all_nt[muts_tumour_all_nt %in% muts_tumour_spec] # 24
# makes sense as the other two are on chr3 G>C and C>G

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 0] # 11
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 1] # 21

# likely shared with all normal samples 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 >= 5] # 31

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 %in% c(2, 3, 4)]
# chr15_49480646_T_A # looks okay (probably PD62341 specific and contamination)
# chr17_18402053_G_C # poor mapping 
# chrX_798031_G_A # looks okay 

######################################################################################################

twins_filtered_vaf[sum_normal == 1 & sum_tumour == 0] # 1
# chr19_54746517_T_C Ql not terrible, 4 reads in PD63383t and generally low everywhere, wouldn't buy it 

twins_filtered_vaf[sum_normal == 2 & sum_tumour == 0] # 4
# chr17_20476029_C_T # poor mapping 
# chr4_15905566_C_T # looks fine (but then, PD62341v and PD63383ae?)
# chr6_159851462_T_C # not real 
# chr7_110349267_G_A # mapping to several places (Blat)

twins_filtered_vaf[sum_normal == 0 & sum_tumour == 1] # 4
# chr3_189240074_C_A am looks okay
# chr3_95470911_T_A ak looks okay
# chr5_157248612_A_G  PD63383aq looks okay
# chr8_91252357_A_C ak looks okay

######################################################################################################
# check my mutation classes again

# create lists of mutations: presence defined based on MTR AND VAF 
muts = twins_filtered_vaf[, mut_ID] %>% unlist()
muts_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist())) 
muts_normal = Reduce(intersect, list(twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist())) 
muts_tumour = Reduce(intersect, list(twins_filtered_vaf[sum_tumour>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist()))
muts_normal_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()))
muts_tumour_all = Reduce(intersect, list(twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()))
muts_tumour_only = Reduce(intersect, list(twins_filtered_vaf[sum_normal==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==0, mut_ID] %>% unlist()))
muts_tumour_all_only = Reduce(intersect, list(muts_tumour_all, muts_tumour_only)) %>% unlist()
muts_normal_only = Reduce(intersect, list(twins_filtered_vaf[sum_tumour==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==0, mut_ID] %>% unlist()))
muts_normal_all_only = Reduce(intersect, list(muts_normal_all, muts_normal_only)) %>% unlist()
muts_normal_all_nt = setdiff(muts_normal_all, muts_all)
muts_tumour_all_nt = setdiff(muts_tumour_all, muts_all)

paste('Mutations present in all samples:', length(muts_all)) # 575
paste('Number of muts present in min 1 normal sample:', length(muts_normal)) # 963
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour)) # 981
paste('Number of muts present in all normal samples:', length(muts_normal_all)) # 632
paste('Number of muts present in all tumour samples:', length(muts_tumour_all)) # 641
paste('Number of muts present in only tumour samples:', length(muts_tumour_only)) # 31
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_only)) # 0
paste('Number of muts present in only normal samples:', length(muts_normal_only)) # 16
paste('Number of muts present in all normal samples and max 1 tumour:', length(muts_normal_all_only)) # 0 # makes sense since tumour is from normal

# ONLY VAF
muts_all_vaf = twins_filtered_vaf[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist()
muts_normal_vaf = twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist()
muts_tumour_vaf = twins_filtered_vaf[sum_tumour>=1, mut_ID] %>% unlist()
muts_normal_all_vaf = twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist()
muts_tumour_all_vaf = twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist()
muts_tumour_only_vaf = twins_filtered_vaf[sum_normal==0, mut_ID] %>% unlist()
muts_tumour_all_only_vaf = Reduce(intersect, list(muts_tumour_all_vaf, muts_tumour_only_vaf)) %>% unlist()
muts_normal_only_vaf = twins_filtered_vaf[sum_tumour==0, mut_ID] %>% unlist()
muts_normal_all_only_vaf = Reduce(intersect, list(muts_normal_all_vaf, muts_normal_only_vaf)) %>% unlist()
muts_normal_all_nt_vaf = setdiff(muts_normal_all_vaf, muts_all_vaf)
muts_tumour_all_nt_vaf = setdiff(muts_tumour_all_vaf, muts_all_vaf)

paste('Mutations present in all samples:', length(muts_all_vaf)) # 616
paste('Number of muts present in min 1 normal sample:', length(muts_normal_vaf)) # 971
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour_vaf)) # 983
paste('Number of muts present in all normal samples:', length(muts_normal_all_vaf)) # 683
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_vaf)) # 685
paste('Number of muts present in only normal samples:', length(muts_normal_only_vaf)) # 19
paste('Number of muts present in only tumour samples:', length(muts_tumour_only_vaf)) # 31
paste('Number of muts present in all normal samples and no tumour:', length(muts_normal_all_only_vaf)) # 0
paste('Number of muts present in all tumour samples and no normal:', length(muts_tumour_all_only_vaf)) # 0
paste('Number of muts present in all normal samples but not all tumour:', length(muts_normal_all_nt_vaf)) # 67
paste('Number of muts present in all tumour samples but not all normal:', length(muts_tumour_all_nt_vaf)) # 69

# ONLY MTR
muts_all_mtr = twins_filtered_mtr[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist()
muts_normal_mtr = twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist()
muts_tumour_mtr = twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist()
muts_normal_all_mtr = twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()
muts_tumour_all_mtr = twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()
muts_tumour_only_mtr = twins_filtered_mtr[sum_normal==0, mut_ID] %>% unlist()
muts_tumour_all_only_mtr = Reduce(intersect, list(muts_tumour_all_mtr, muts_tumour_only_mtr)) %>% unlist()
muts_normal_only_mtr = twins_filtered_mtr[sum_tumour==0, mut_ID] %>% unlist()
muts_normal_all_only_mtr = Reduce(intersect, list(muts_normal_all_mtr, muts_normal_only_mtr)) %>% unlist()
muts_normal_all_nt_mtr = setdiff(muts_normal_all_mtr, muts_all_mtr)
muts_tumour_all_nt_mtr = setdiff(muts_tumour_all_mtr, muts_all_mtr)

paste('Mutations present in all samples:', length(muts_all_mtr)) # 581
paste('Number of muts present in min 1 normal sample:', length(muts_normal_mtr)) # 963
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour_mtr)) # 984
paste('Number of muts present in all normal samples:', length(muts_normal_all_mtr)) # 637
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_mtr)) # 649
paste('Number of muts present in only tumour samples:', length(muts_tumour_only_mtr)) # 39
paste('Number of muts present in only normal samples:', length(muts_normal_only_mtr)) # 18
paste('Number of muts present in all normal samples and no tumour:', length(muts_normal_all_only_mtr)) # 0
paste('Number of muts present in all tumour samples and no normal:', length(muts_tumour_all_only_mtr)) # 1
paste('Number of muts present in all normal samples but not all tumour:', length(muts_normal_all_nt_mtr)) # 56
paste('Number of muts present in all tumour samples but not all normal:', length(muts_tumour_all_nt_mtr)) # 68

######################################################################################################
# I think out of the ~60-70 mutations that are in all tumour but not all normal, we only placed 6 on the phylogenetic tree

sum(muts_tumour_all_nt %in% muts_tumour_spec) # 24
# okay so 24 of these are classified as tumour mutations
# therefore, 24 of those are assigned to tumour-specific mutations
# we know that other 8 mutations are shared b/n PD62341 normal and tumour 

# we can try to place the other mutations on this 
muts_tumour_all_nt_unplaced = setdiff(muts_tumour_all_nt, muts_tumour_spec)
muts_tumour_all_nt_unplaced = setdiff(muts_tumour_all_nt_unplaced, muts_PD62341) 
paste('Number of yet-unplaced mutations:', length(muts_tumour_all_nt_unplaced)) # 34

# always present in 10 tumour samples

# usually present in 10-12 tumour samples 
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced, sum_normal],
     xlab = 'Nr of normal samples',  main = '34 muts (unplaced tumour not all normal)')

twins_filtered_vaf[, missing_normal_sample := apply(twins_filtered_vaf[, ..samples_normal_vaf], 1, function(row) {
  cols_vaf01 = samples_normal_vaf[which(row < 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

twins_filtered_vaf[, present_normal_sample := apply(twins_filtered_vaf[, ..samples_normal_vaf], 1, function(row) {
  cols_vaf01 = samples_normal_vaf[which(row >= 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

# present in different nr of normal samples
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==4] # 1
# chr5_44907911_G_A q, aa, bb, h (could be contamination, but also 0.03 in ad and 0.06 in n, so maybe a late one)
# looks very believable, I'd be keen to buy this 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==5] # 0 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==6] # 0 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==7] # 0 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==8] # 1
# chr17_21307477_G_A missing in PD62341ad, h, n, PD63383w (is there at VAF 0.05-0.09)
# mapping doesn't give me a lot of confidence in this one 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==9] # 2
# chr17_18402053_G_C missing in PD63383w, ae, PD62341n (is there at VAF > 0.08 tho) # mapping isn't excellent, plausible
# chrX_798031_G_A missing in PD63383u, bb, PD62341h (is there at VAF 0.05 - 0.09)
# the thing with this one is that it looks believable but is in all normal samples at a really low VAF?

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==10] # 3
# chr11_4220977_T_C missing in PD62341q (0.03), PD63383ak (0.04) (mut only on poorly mapped reads)
# chr19_90678_G_A missing in PD62341n (0.087), PD63383bb (0.09) - likely present in all (looks super bad)
# chr5_176771316_G_A missing in PD62341n (0), PD62341aa (0.085) - I can see reads in PD62341n on Jbrowse tho, looks okay?

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==11] # 12
# chr11_4594537_A_G # poor mapping 
# chr17_36165401_A_G # poor mapping
# chr17_36388753_C_G # poor mapping
# chr1_157363251_C_T # looks real (0.08 in PD63383t)
# chr5_176771327_T_A # looks real (0.08 in PD62341n)
# chr6_32384367_G_A # looks quite plausible (0.04 in PD62341h)
# chr7_102544208_T_C # mapping is questionable
# chr7_55402890_G_T # poor mapping 
# chr7_64978551_T_C # looks okay (0.04 in PD62341h)
# chr7_76798258_T_C # looks okay (0.0976 in PD62341ad, likely everywhere, mapping could be better)
# chr8_7418832_G_C # looks kind of okay (0.087 in PD62341q)
# chr9_42329754_C_T # very questionable mapping 

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==12]
# 15 

######################################################################################################
# how many mutations do I have with VAF 0.1 - 0.2 in agg and 0.1 - 0.2 in both twins?
muts_low_vaf_both = twins_agg_vaf[vaf_normal_all > 0.1 & vaf_normal_all < 0.2 & vaf_normal_PD62341 > 0.1 & vaf_normal_PD62341 < 0.2 & vaf_normal_PD63383 > 0.1 & vaf_normal_PD63383 < 0.2, mut_ID] %>% unlist() # 59 
# wow there are actually 59 of those!

twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==4] # 1
# chr8_124634298_G_T (useless)
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==5] # 1
# chr14_19049856_C_T (useless)
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==6] # 0
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==7] # 6
# chr17_21439338_C_T # poor mapping 
# chr1_146568235_A_G # poor mapping 
# chr5_100144922_C_T # poor mapping 
# chr6_57048152_G_A # poor mapping 
# chr7_56823245_C_G # poor mapping 
# chrX_115889128_C_A # poor mapping 
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==8]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==9]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==10]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==11]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==12]
# chr11_4220746_T_C # mapping could be better?
# chr12_44087006_G_A # double check but possible insertions
# chr19_40853394_C_A # poor mapping 
# chr3_190409217_C_T # looks real! > likely comes from the cluster on chr3
# chr6_94069453_A_T # insertions

# the real stuff seems to map back to clusters seen for mutations present in all samples
# perhaps they are just at lower VAF due to sampling (dropout in some samples etc.)

######################################################################################################
# are there any drivers of clonal hem in the set of normal mutations which would explain uneven transfusion b/n twins?

table(twins_dt[mut_ID %in% muts_normal_all & Gene != '-' & Effect != 'intronic', Gene])
# AGAP6, CCL3L3, CCL4L2, CROCC, CSH2, CYP2A6, DEFB4B, DUOX1, DUSP22 ,EIF5AL1, ERVV-2, GOLGA6L9, GSTT4
# GTF3C5, HERC2, KCNJ18, KIR2DL1, KRTAP10-4, LHFPL5, LILRA2, MAGEA8, MIR3680-1, NQO1, OR2T3, OR4F4, OR4K2
# OR4M1, POMZP3, POTEH, PRAMEF18, RNU1-2, SIMC1, SLC35E2B, TCAF1, THOC3, ZDHHC11B, ZNF705G

# clonal hem driver genes (from https://pmc.ncbi.nlm.nih.gov/articles/PMC11176083/)
# ZBTB33, ZNF318, ZNF234, SPRED2, SH2B3, SRCAP, SIK3, SRSF1, CHEK2, CCDC115, CCL22, BAX, YLPM1, MYD88, MTA2, MAGEC3 and IGLL5
# DNMT3A, TET2, ASXL1, PPM1D, SRSF2, SF3B1, GNB1, IDH1, IDH2, TP53, BRCC2, GNAS, JAK2, KDM6A, CBL, PHIP
# okay so it doesn't seem I have any CH drivers in my dataset 

CH_genes = c('ZBTB33', 'ZNF318', 'ZNF234', 'SPRED2', 'SH2B3', 'SRCAP', 'SIK3', 'SRSF1', 'CHEK2', 'CCDC115', 'CCL22', 'BAX', 'YLPM1', 'MYD88', 'MTA2', 'MAGEC3',
             'IGLL5', 'DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'SRSF2', 'SF3B1', 'GNB1', 'IDH1', 'IDH2', 'TP53', 'BRCC2', 'GNAS', 'JAK2', 'KDM6A', 'CBL', 'PHIP')

sum(twins_dt[Gene != '-' & Effect != 'intronic', Gene] %>% unlist() %in% CH_genes) # 42 apparently 

twins_dt_filters[Effect != 'intronic' & Gene %in% CH_genes & sum_req_filters6 <= 1, c('mut_ID', samples_vaf, 'Effect', 'Gene',
                                                                                      'f7_likelyGermline_PD63383', 'f7_likelyGermline_PD62341',                                                                                'f7_likelyGermline_bothTwins', 'f7_likelyGermline_aggTwins'), with=FALSE]
# mutations in: SH2B3, YLPM1, IDH2, SRCAP2, CCL22, SRSF1, PPM1D, TP53, ZNF234, ASXL1, CHEK2, CCDC115, IDH1, DNMT3A, JAK2, MAGEC3
# all look like they are germline regardless for how you test for this 

######################################################################################################
# Plot mutation signatures for specific sets of mutations 

get_trinucs <- function(mybed, genome) {
  
  mybed$order <- 1:nrow(mybed) # order that rows are supplied in 
  
  gr <- GRanges(seqnames = mybed$Chrom, IRanges(start = mybed$Pos, width=1), ref=mybed$Ref, alt=mybed$Alt, order=mybed$order)
  # create a GRanges object with coordinates for the mutation
  
  if (all(substr(seqlevels(gr), 1, 3) != "chr")) {
    gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
  }
  
  seqinfo(gr) <- seqinfo(genome)[seqlevels(gr)]
  gr <- sort(gr)
  bases <- c("A", "C", "G", "T")
  trinuc_levels <- paste0(rep(bases, each = 16), rep(rep(bases, each = 4), 4), rep(bases, 16))
  
  get_trinuc <- function(seqname) {
    pos <- start(gr[seqnames(gr) == seqname])
    view <- Views(genome[[seqname]], start = pos - 1, end = pos + 1)
    ans <- factor(as.character(view), levels = trinuc_levels, labels = 1:64)
    return(as.numeric(ans))
  }
  
  trinuc <- sapply(seqlevels(gr), get_trinuc)
  gr$trinuc <- factor(unlist(trinuc, use.names = FALSE), levels = 1:64, labels = trinuc_levels)
  remove(trinuc)
  gr$REF <- gr$ref
  gr$ALT <- gr$alt
  gr$context <- gr$trinuc
  torc <- which(gr$ref %in% c("A", "G"))
  gr$REF[torc] <- as.character(reverseComplement(DNAStringSet(gr$REF[torc])))
  gr$ALT[torc] <- as.character(reverseComplement(DNAStringSet(gr$ALT[torc])))
  gr$context[torc] <- as.character(reverseComplement(DNAStringSet(gr$context[torc])))
  gr$class <- paste(gr$REF, gr$ALT, "in", gr$context, sep = ".")
  class_levels <- paste(rep(c("C", "T"), each = 48), rep(c("A", "G", "T", "A", "C", "G"), each = 16), "in", paste0(rep(rep(bases, each = 4), 6), rep(c("C", "T"), each = 48), rep(bases, 24)), sep = ".")
  gr$class <- factor(gr$class, levels = class_levels)
  
  # match them up with one another using the original order that I put in - I think that they may have been reshuffled.
  grdf <- as.data.frame(gr)
  grdf <- grdf[with(grdf, order(order)),]
  return(grdf$class)
} 

# Plot the distribution of mutation types on the unfiltered dataset (362,814 mutations)
mybed_muts_tumour = twins_dt[mut_ID %in% muts_tumour_spec, c('Chrom', 'Pos', 'Ref', 'Alt')]
mybed_muts_putgerm = twins_dt[mut_ID %in% muts_normal_all, c('Chrom', 'Pos', 'Ref', 'Alt')]
mybed_muts_PD63383_tumour = twins_dt[mut_ID %in% muts_PD63383_tumour_only, c('Chrom', 'Pos', 'Ref', 'Alt')]

trins_muts_tumour = get_trinucs(mybed_muts_tumour, BSgenome.Hsapiens.UCSC.hg38)
trins_muts_putgerm = get_trinucs(mybed_muts_putgerm, BSgenome.Hsapiens.UCSC.hg38)
trins_muts_PD63383_tumour = get_trinucs(mybed_muts_PD63383_tumour, BSgenome.Hsapiens.UCSC.hg38)

dt_tumour = twins_dt[mut_ID %in% muts_tumour_spec]
dt_tumour$trins_muts_tumour=trins_muts_tumour

# plot the distribution of different mutations across different contexts 
mut_sign_counts = data.table(table(dt_tumour[, trins_muts_tumour]))
setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = '144 mutations (tumour)')+
  theme_minimal(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241109_p1_mut_trins_allmuts.pdf', width = 10, height = 5.5)

trins_muts_putgerm = get_trinucs(mybed_muts_putgerm, BSgenome.Hsapiens.UCSC.hg38)
trins_muts_PD63383_tumour = get_trinucs(mybed_muts_PD63383_tumour, BSgenome.Hsapiens.UCSC.hg38)

dt_pg = twins_dt[mut_ID %in% muts_normal_all]
dt_pg$trins_muts_pg=trins_muts_putgerm

# plot the distribution of different mutations across different contexts 
mut_sign_counts = data.table(table(dt_pg[, trins_muts_pg]))
setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = '632 mutations (all normal samples)')+
  theme_minimal(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241109_p1_mut_trins_allmuts.pdf', width = 10, height = 5.5)

dt_pd = twins_dt[mut_ID %in% muts_PD63383_tumour_only]
dt_pd$trins_muts_pd=trins_muts_PD63383_tumour

# plot the distribution of different mutations across different contexts 
mut_sign_counts = data.table(table(dt_pd[, trins_muts_pd]))
setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = '28 mutations (PD63383 tumour)')+
  theme_minimal(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241109_p1_mut_trins_PD63383tumour.pdf', width = 10, height = 5.5)

######################################################################################################
# Plot PD62341 and PD63383 specific mutations

muts_PD62341_plausible = c("chr15_56691722_C_T", "chr18_71485446_G_A", "chr19_43132229_T_C",
                           "chr20_44114996_C_T", "chr2_57814739_A_C", "chr7_120677593_C_T",
                           "chr7_139659050_G_A", "chr8_131571989_A_T", "chr9_12824689_G_A", 
                           "chrX_124531709_C_T")

muts_PD63383 = c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                 "chr13_50815806_A_G", "chr4_74625500_G_T", "chr7_73831920_C_T",
                 "chr3_77633967_C_T", "chrX_115066661_C_T", "chr7_149688370_C_T",
                 'chr1_103587565_A_C', 'chr21_40193588_G_A', 'chr3_165901319')

muts_PD62341 = c("chr14_105458006_C_A", "chr17_33422229_C_A", "chr15_49480646_T_A",
                 "chr16_5479739_C_T", "chr2_95662131_G_A", "chr3_50106043_C_T",
                 "chr3_62055057_C_G", "chr3_62055077_G_C") 

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

ggplot(mut_PD63383_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_PD62341, '#C42F0F', col_PD63383, '#F99A49'))+
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  ggtitle(glue('PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_p3_vaf_dist_PD63383_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

ggplot(mut_PD62341_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_PD62341, '#C42F0F', col_PD63383, '#F99A49'))+
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  ggtitle(glue('PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_p3_vaf_dist_PD62341_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

# Is there a way I can quantify this 
# Maybe first, compare mean VAF for 5 non-spleen samples and spleen 
means_PD62341muts_spleen_PD62341 = mut_PD62341_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD62341 = data.table(mut_PD62341_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD62341muts_spleen_PD63383 = mut_PD62341_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD63383 = data.table(mut_PD62341_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD62341muts = merge(merge(means_PD62341muts_spleen_PD62341, means_PD62341muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD62341muts_spleen_PD63383, means_PD62341muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD62341muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))

# I would do the calculations for mutations where VAF in PD63383 non-spleen is 0 for now
# VAF_PD62341_spleen = VAF_PD62341_spleen * fraction_PD62341_spleen + VAF_PD63383 * fraction_PD63383 
# VAF_PD63383_spleen = VAF_PD62341 * fraction_PD62341 + VAF_PD63383_spleen * fraction_PD63383_spleen
means_PD62341muts_clean = means_PD62341muts[PD63383_nonspleen==0]
means_PD62341muts_clean[, PD62341_spleen_fPD62341 := PD62341_spleen / PD62341_nonspleen]
means_PD62341muts_clean[, PD63383_spleen_fPD62341 := PD63383_spleen / PD62341_nonspleen]

ggplot(means_PD62341muts_clean, aes(x=PD62341_nonspleen, y=PD62341_spleen))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'mean VAF PD62341 (non-spleen)', y = 'VAF PD62341 (spleen)')+
  ggtitle(glue('PD62341-specific mutations'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.5))+
  ylim(c(0, 0.5))
ggsave(glue('Results/20241109_spleen_vs_nonspleen_PD62341.pdf'), width=4, height=4)

# do the same for PD63383 specific mutations
means_PD63383muts_spleen_PD62341 = mut_PD63383_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD62341 = data.table(mut_PD63383_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD63383muts_spleen_PD63383 = mut_PD63383_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD63383 = data.table(mut_PD63383_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD63383muts = merge(merge(means_PD63383muts_spleen_PD62341, means_PD63383muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD63383muts_spleen_PD63383, means_PD63383muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD63383muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))
means_PD63383muts[, PD62341_spleen_fPD63383 := PD62341_spleen / PD63383_nonspleen] # contamination
means_PD63383muts[, PD62341_spleen_fPD63383_2 := PD62341_spleen / PD63383_spleen] # contamination
means_PD63383muts[, PD63383_spleen_fPD63383 := PD63383_spleen / PD63383_nonspleen] # purity 

# plots of VAF for twin-specific mutations
ggplot(mut_PD62341_melt[status=='normal'], aes(x=value, y=mut_ID, color=sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD62341-specific mutations'))+
  xlim(c(0, 0.7))
ggsave(glue('Results/20241109_PD62341_specific_vaf.pdf'), width=8, height=5)

ggplot(mut_PD63383_melt[status=='normal'], aes(x=value, y=mut_ID, color=sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD63383-specific mutations'))+
  xlim(c(0, 0.7))
ggsave(glue('Results/20241109_PD63383_specific_vaf.pdf'), width=8, height=5)

mut_PD62341_melt[, sample_type4 := as.factor(fcase(
  status == 'tumour', 'tumour',
  twin == 'PD63383' & sample %in% samples_normal, 'PD63383 normal',
  twin == 'PD62341' & sample %in% samples_normal, 'PD62341 normal'
))]
mut_PD63383_melt[, sample_type4 := as.factor(fcase(
  status == 'tumour', 'tumour',
  twin == 'PD63383' & sample %in% samples_normal, 'PD63383 normal',
  twin == 'PD62341' & sample %in% samples_normal, 'PD62341 normal'
))]

ggplot(mut_PD62341_melt, aes(x=value, y=mut_ID, color=sample_type4))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD62341-specific mutations'))+
  xlim(c(0, 0.7))
ggsave(glue('Results/20241109_PD62341_specific_vaf2.pdf'), width=8, height=5)

ggplot(mut_PD63383_melt, aes(x=value, y=mut_ID, color=sample_type4))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD63383-specific mutations'))+
  xlim(c(0, 0.7))
ggsave(glue('Results/20241109_PD63383_specific_vaf2.pdf'), width=8, height=5)

mut_PD62341_dt2 = twins_vaf[mut_ID %in% muts_PD62341_plausible, 1:23]
mut_PD62341_melt2 = melt(mut_PD62341_dt2, id.vars = 'mut_ID')
mut_PD62341_melt2[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD62341_melt2[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_PD62341_melt2[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD62341_melt2[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD62341_melt2[, sample_type4 := as.factor(fcase(
  status == 'tumour', 'tumour',
  twin == 'PD63383' & sample %in% samples_normal, 'PD63383 normal',
  twin == 'PD62341' & sample %in% samples_normal, 'PD62341 normal'
))]

ggplot(mut_PD62341_melt2[status=='normal'], aes(x=value, y=mut_ID, color=sample_type4))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD62341-specific plausible mutations'))
ggsave(glue('Results/20241109_PD62341_specific_vaf_plausible.pdf'), width=8, height=5)

ggplot(mut_PD62341_melt2, aes(x=value, y=mut_ID, color=sample_type4))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD62341-specific plausible mutations'))
ggsave(glue('Results/20241109_PD62341_specific_vaf_plausible2.pdf'), width=8, height=5)

######################################################################################################
# v quickly is there anything else plausible from the germline filters

muts_germ_PD62341 = twins_dt_filters[sum_req_filters_qual == 0 & f7_likelyGermline_PD62341 == 1 & f7_likelyGermline_PD63383 == 0, mut_ID] %>% unlist()
muts_germ_PD63383 = twins_dt_filters[sum_req_filters_qual == 0 & f7_likelyGermline_PD62341 == 0 & f7_likelyGermline_PD63383 == 1, mut_ID] %>% unlist()

dt_germ_PD62341 = twins_agg_vaf[mut_ID %in% muts_germ_PD62341, c('mut_ID', 'vaf_normal_PD62341', 'vaf_normal_PD63383', 'vaf_tumour_all'), with=FALSE]
dt_germ_PD63383 = twins_agg_vaf[mut_ID %in% muts_germ_PD63383, c('mut_ID', 'vaf_normal_PD62341', 'vaf_normal_PD63383', 'vaf_tumour_all'), with=FALSE]
dt_germ_PD62341[, vaf_ratio := vaf_normal_PD62341 / vaf_normal_PD63383]
dt_germ_PD63383[, vaf_ratio := vaf_normal_PD62341 / vaf_normal_PD63383]
dt_germ_PD62341[, vaf_ratio_log :=log2(vaf_ratio)]
dt_germ_PD63383[, vaf_ratio_log :=log2(vaf_ratio)]

dt_germ_PD62341[abs(vaf_ratio_log) > 0.4, mut_ID]
# "chr10_13474524_G_A" # deletions  
# "chr14_105458006_C_A" # have it 
# "chr15_49480646_T_A" # have it   
# "chr15_82483067_G_A" # poor mapping 
# "chr16_5479739_C_T" # have it   
# "chr17_33422229_C_A" # have it 
# "chr1_16864604_G_A" # poor mapping    
# "chr1_248477453_G_A" # maybe - not very convinced
# "chr2_95662131_G_A"  # have it  
# "chr3_50106043_C_T" # have it   
# "chr6_149763201_T_A" # deletions  
# "chr9_23627511_T_G" # have it 
dt_germ_PD63383[abs(vaf_ratio_log) > 0.4, mut_ID]
# "chr12_129105967_T_C" # poor mapping
# "chr13_69946674_A_G"  # poor mapping
# "chr15_55075682_T_C" # looks good!
# "chr16_22670560_C_T" # also looks good
# "chr17_21644661_A_G"  # poor mapping 
# "chr1_38827952_C_A" # have this one   
# "chr1_76183938_A_G" # oh this one is super interesting!  
# "chr6_55120674_A_G" # deletions 
# "chr7_76547796_C_T" # poor mapping   
# "chr7_94839706_A_T" # looks okay?

muts_PD63383_plausible = c('chr15_55075682_T_C', 'chr7_94839706_A_T', 
                           'chr1_76183938_A_G', 'chr16_22670560_C_T')
mut_PD63383_dt2 = twins_vaf[mut_ID %in% muts_PD63383_plausible, 1:23]
mut_PD63383_melt2 = melt(mut_PD63383_dt2, id.vars = 'mut_ID')
mut_PD63383_melt2[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD63383_melt2[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_PD63383_melt2[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD63383_melt2[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD63383_melt2[, sample_type4 := as.factor(fcase(
  status == 'tumour', 'tumour',
  twin == 'PD63383' & sample %in% samples_normal, 'PD63383 normal',
  twin == 'PD62341' & sample %in% samples_normal, 'PD62341 normal'
))]

ggplot(mut_PD63383_melt2[status=='normal'], aes(x=value, y=mut_ID, color=sample_type4))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD63383-specific plausible mutations'))
ggsave(glue('Results/20241109_PD63383_specific_vaf_plausible.pdf'), width=8, height=5)

ggplot(mut_PD63383_melt2, aes(x=value, y=mut_ID, color=sample_type4))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  ggtitle(glue('PD63383-specific plausible mutations'))
ggsave(glue('Results/20241109_PD63383_specific_vaf_plausible2.pdf'), width=8, height=5)

######################################################################################################
# CHECK: ARE THERE MUTATIONS PRESENT IN A SPECIFIC TISSUE OF BOTH TWINS?

# search through samples which are from the same tissue 
samples_spleen = c('PD62341v_VAF', 'PD63383w_VAF') 
samples_skin = c('PD62341aa_VAF', 'PD63383bb_VAF') 
samples_cerebellum = c('PD62341ad_VAF', 'PD63383ak_VAF') 
samples_liver = c('PD62341h_VAF', 'PD63383t_VAF') 
samples_heart = c('PD62341n_VAF', 'PD63383ae_VAF') 
samples_pancreas = c('PD62341q_VAF', 'PD63383u_VAF') 

twins_filtered_vaf[, mean_vaf_spleen := apply(.SD, 1, mean), .SDcols = samples_spleen]
twins_filtered_vaf[, mean_vaf_skin := apply(.SD, 1, mean), .SDcols = samples_skin]
twins_filtered_vaf[, mean_vaf_cerebellum := apply(.SD, 1, mean), .SDcols = samples_cerebellum]
twins_filtered_vaf[, mean_vaf_liver := apply(.SD, 1, mean), .SDcols = samples_liver]
twins_filtered_vaf[, mean_vaf_heart := apply(.SD, 1, mean), .SDcols = samples_heart]
twins_filtered_vaf[, mean_vaf_pancreas := apply(.SD, 1, mean), .SDcols = samples_pancreas]
twins_filtered_vaf[, mean_vaf_normal := apply(.SD, 1, mean), .SDcols = samples_normal_vaf]

# you can also specify sample order and plot them as a heatmap - that will be easy to see 
mut_tissues_nt = twins_filtered_vaf[, c('mut_ID', samples_normal_vaf), with=FALSE] # only include mutations that passed your filter 
setcolorder(mut_tissues_nt, c('mut_ID', 'PD62341ad_VAF', 'PD63383ak_VAF', 
                              'PD62341aa_VAF', 'PD63383bb_VAF',
                              'PD62341q_VAF', 'PD63383u_VAF',
                              'PD62341h_VAF', 'PD63383t_VAF',
                              'PD62341n_VAF', 'PD63383ae_VAF',
                              'PD62341v_VAF', 'PD63383w_VAF'))
mut_tissues_nt_mat = as.matrix(mut_tissues_nt[, c(samples_normal_vaf), with=FALSE])
rownames(mut_tissues_nt_mat) = mut_tissues_nt[,1] %>% unlist()  
colnames(mut_tissues_nt_mat) = tstrsplit(colnames(mut_tissues_nt_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Twin = rep(c('PD62341', 'PD63383'), 6), 
                            Tissue = rep(c('cerebellum', 'skin', 'pancreas', 'liver', 'heart', 'spleen'), each = 2))
rownames(col_annotation) = colnames(mut_tissues_nt_mat)
annotation_colors = list(Twin = c(PD62341=col_PD62341, PD63383=col_PD63383), 
                         Tissue = c(cerebellum = '#0957cc', skin = '#8accfc', 
                                    pancreas = '#bc7f0c', liver = '#e4bf28', 
                                    heart = '#d1140d', spleen = '#f1732a'))

# heatmap
pdf('Results/20241109_p6_heatmap_compare_same_tissues.pdf', heigh = 20, width = 10)
pheatmap(mut_tissues_nt_mat,
         cellwidth=40, cellheight=1,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="1,002 mutations", 
         legend = T, 
         treeheight_row = 0, 
         cluster_rows = T, cluster_cols = F, 
         show_rownames = F, show_colnames = T,
         fontsize=20, cexCol=2) 
dev.off()

# search for mutations which seem to be enriched in a given tissue over the others 
twins_filtered_vaf[, ratio_spleen_normal := mean_vaf_spleen / mean_vaf_normal]
twins_filtered_vaf[, ratio_skin_normal := mean_vaf_skin / mean_vaf_normal]
twins_filtered_vaf[, ratio_cerebellum_normal := mean_vaf_cerebellum / mean_vaf_normal]
twins_filtered_vaf[, ratio_liver_normal := mean_vaf_liver / mean_vaf_normal]
twins_filtered_vaf[, ratio_heart_normal := mean_vaf_heart / mean_vaf_normal]
twins_filtered_vaf[, ratio_pancreas_normal := mean_vaf_pancreas / mean_vaf_normal]

twins_filtered_vaf[, ratio_spleen_normal_log := log2(ratio_spleen_normal)]
twins_filtered_vaf[, ratio_skin_normal_log := log2(ratio_skin_normal)]
twins_filtered_vaf[, ratio_cerebellum_normal_log := log2(ratio_cerebellum_normal)]
twins_filtered_vaf[, ratio_liver_normal_log := log2(ratio_liver_normal)]
twins_filtered_vaf[, ratio_heart_normal_log := log2(ratio_heart_normal)]
twins_filtered_vaf[, ratio_pancreas_normal_log := log2(ratio_pancreas_normal)]

# can I identify any instances where the log is higher?
# spleen
twins_filtered_vaf[ratio_spleen_normal_log > 1 & PD62341v_VAF > 0.1 & PD63383w_VAF > 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr2_110354557_G_T does not look real (maps everywhere)
twins_filtered_vaf[ratio_spleen_normal_log < -1 & ratio_spleen_normal != 0 & PD62341v_VAF < 0.1 & PD63383w_VAF < 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr14_19049856_C_T poor quality 
# chr14_73648334_G_A looks real, follow it up 
# chr19_54861143_C_A generally looks real, but it does have some split reads
# chr1_25421838_C_T only on badly mapped reads and out of phase with what looks like germline
# chr2_95747706_T_C not ideal, maybe 

# skin
twins_filtered_vaf[ratio_skin_normal_log > 1 & PD62341aa_VAF > 0.1 & PD63383bb_VAF > 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr6_57048152_G_A # poor mapping 
# chr8_124634298_G_T # terrible 
twins_filtered_vaf[ratio_skin_normal_log < -1 & ratio_skin_normal != 0 & PD62341aa_VAF < 0.1 & PD63383bb_VAF < 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr10_54938275_C_A # insertions 
# chr1_1041773_G_T # insertions 
# chr2_124431547_T_C # not convinced - some deletions, very difficult context 
# chr5_175963476_G_C # poor mapping 

# cerebellum
twins_filtered_vaf[ratio_cerebellum_normal_log > 1 & PD62341ad_VAF > 0.1 & PD63383ak_VAF > 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
twins_filtered_vaf[ratio_cerebellum_normal_log < -1 & ratio_cerebellum_normal != 0 & PD62341ad_VAF < 0.1 & PD63383ak_VAF < 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr17_21439338_C_T # poor mapping 
# chr18_56262973_G_A # questionable - deletions + poor mapping 
# chr22_36255758_A_G # poor mapping and very low levels 
# chr2_124431547_T_C # not buying it given deletions and context 
# chr5_100144922_C_T # poor mapping  
# chr6_57048152_G_A # poor mapping
# chr8_124634298_G_T # terrible 
# chrX_115889128_C_A # poor mapping 

# liver 
twins_filtered_vaf[ratio_liver_normal_log > 1 & PD62341h_VAF > 0.1 & PD63383t_VAF > 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
twins_filtered_vaf[ratio_liver_normal_log < -1 & ratio_liver_normal != 0 & PD62341h_VAF < 0.1 & PD63383t_VAF < 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr12_11292597_C_A # poor mapping 
# chr12_18986629_G_A # not sure, maybe?
# chr14_19657092_A_G # poor mapping 
# chr1_1446554_A_G # poor mapping 
# chr1_148153535_T_C # poor mapping 
# chr20_54999345_T_C # quite likely it's deletions 
# chr5_176010320_A_G # poor mapping 
# chrX_106521135_C_A # insertions 

# heart 
twins_filtered_vaf[ratio_heart_normal_log > 1 & PD62341n_VAF > 0.1 & PD63383ae_VAF > 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr18_56262973_G_A # poor mapping 
# chr1_148153535_T_C # poor mapping 
# chr5_175963476_G_C # poor mapping 
twins_filtered_vaf[ratio_heart_normal_log < -1 & ratio_heart_normal != 0 & PD62341n_VAF < 0.1 & PD63383ae_VAF < 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr14_93783931_G_A # insertions
# chr15_22572377_T_C # poor mapping 
# chr5_176010320_A_G # poor mapping 

# pancreas
twins_filtered_vaf[ratio_pancreas_normal_log > 1 & PD62341q_VAF > 0.1 & PD63383u_VAF > 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr2_124431547_T_C # insertions 
twins_filtered_vaf[ratio_pancreas_normal_log < -1 & ratio_pancreas_normal != 0 & PD62341q_VAF < 0.1 & PD63383u_VAF < 0.1 & sum_normal_PD62341 > 1 & sum_normal_PD63383 > 1]
# chr12_11292597_C_A # poor mapping 
# chr17_1199637_C_G # poor mapping 
# chr18_56262973_G_A # poor mapping 
# chr19_56105967_G_A # poor mapping 
# chr2_95747706_T_C # poor mapping (not terrible but not good enough)
# chr5_71157247_C_T # poor mapping 
# chr6_57048152_G_A # poor mapping 

######################################################################################################
# SEARCH FOR PROTEIN-AFFECTING MUTATIONS PRESENT IN THE TUMOUR: POSSIBLE DRIVERS  
muts_coding_subs_tumour = sort(twins_dt_filters[mut_ID %in% muts_tumour_spec & Gene != '-' & Effect %in% c('ess_splice', 'missense', 'nc_variant', 'nonsense', 'start_lost', 'stop_lost'), mut_ID] %>% unlist())
paste('Number of protein-altering mutations specific to tumour samples:', length(muts_coding_subs)) # 0 

muts_ingene_subs_tumour = sort(twins_dt_filters[mut_ID %in% muts_tumour_spec & Gene != '-', mut_ID] %>% unlist())
paste('Number of mutations possibly affecting proteins specific to tumour samples:', length(muts_ingene_subs_tumour)) # 43 

genes_muts_subs_tumour = sort(twins_dt_filters[mut_ID %in% muts_tumour_spec & Gene != '-', Gene] %>% unique() %>% unlist())
paste('Number of genes with protein-altering mutations specific to tumour samples:', length(genes_muts_subs_tumour)) # 43 

paste('Number of driver genes screened for:', length(driver_genes)) # 750 
paste('Number of genes with tumour specific mutations that are plausible drivers:', length(Reduce(intersect, list(driver_genes, genes_muts_subs_tumour)))) # 0

######################################################################################################
# You have 1,002 mutations. For each of them, I want you to be able to say what class it belongs to or if it is an artefact 


