###################################################################################################################################
# SCRIPT 7

# Script to analyse the scRNA-seq genotyping allele integrator output 
# December 2024
# Barbara Walkowiak bw18

# INPUT: 
# 1 .tsv dataframes generated from alleleIntegrator
# allele integrator ran 04/12/2024; jobID 664854
# downloaded files from the farm via rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scRNA_genotyping/output/" /Users/bw18/Desktop/1SB/Data/scRNAseq/
# 2 data frame with final 255 mutations (ID + assignment to informative group on the phylogeny)
# 3 Seurat object with clusters (generated from tumour from either twin separately)

# OUTPUT
# 1 UMAP plots with presence of each mutation (mutant / wt reads) across different clusters

###################################################################################################################################
# LIBRARIES 

# Load needed libraries
library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)
library(gridExtra)

###################################################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
col_mutant = '#ad0505'
col_wt = '#07a7d0'
col_tumour = '#fed4a3'
col_background = '#c8c8c8'

# set colors so mutant / wt always plotted in the same color 
myCols = c('mutant read' = col_mutant, 'wt read' = col_wt, "normal cells" = col_background, 'tumour cells' = col_tumour)

###################################################################################################################################
# INPUT FILES 

setwd('/Users/bw18/Desktop/1SB')

# Read the merged dataframe with the CAVEMAN pileup 
twins_dt = fread('Data/pileup_merged_20241016.tsv') # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Out/F1/F1_mutations_final_20241208_255.txt') 
muts = muts$V1 %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 255

# Create dataframe with mutations that passed current filters 
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# load AlleleIntegrator output from all samples 
samples_scRNAseq = c('NB13652544', 'NB13652545', 'NB13652546',
                     'NB13760628', 'NB13760629', 'NB13760630')
                     # 'NB14599068', 'NB14599069', 'NB14599070') # exclude nuclear RNA-seq - we established it is very poor quality (only 194 cells captured)

ai_counts_dt = data.table()
for (sample in samples_scRNAseq){
  ai_counts = fread(paste0('Data/scRNAseq/', sample, '_somaticVar_scRNA_alleleCounts.tsv'), sep='\t', fill = TRUE, header=TRUE)
  ai_counts[, sample_ID := sample]
  ai_counts_dt = rbind(ai_counts_dt, ai_counts)
}

aiissue = fread('Data/scRNAseq/NB13760629_somaticVar_scRNA_alleleCounts.tsv', sep = '\t')

# load Seurat object with UMAP-based clustering 
tumour_PD62341_clusters = readRDS(file = "Out/F6/F6_tPD62341.agg.rds")
tumour_PD63383_clusters = readRDS(file = "Out/F6/F6_tPD63383.agg.rds")

# load list of final mutations used for the phylogeny
muts_assignment = fread('Out/F3/F3_muts_classes_255_20241208.csv', sep = ',', header = TRUE)

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
# CALCULATED USEFUL VALUES (AGGREGATED VAF ETC) FOR PILEUP FILTERED FOR 255 FINAL MUTATIONS

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
twins_filtered_dt[, sum_all_mtr := rowSums(.SD>=4), .SDcols = c(samples_tumour_mtr, samples_normal_mtr)]
twins_filtered_dt[, sum_tumour_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, sum_normal_mtr := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_dt[, sum_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_dt[, sum_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_dt[, sum_tumour_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, sum_tumour_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, sum_normal_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, sum_normal_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# Add presence / absence based on VAF data
twins_filtered_dt[, sum_all_vaf := rowSums(.SD>=0.1), .SDcols = c(samples_tumour_vaf, samples_normal_vaf)]
twins_filtered_dt[, sum_tumour_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_filtered_dt[, sum_normal_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_vaf]
twins_filtered_dt[, sum_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD62341_vaf]
twins_filtered_dt[, sum_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD63383_vaf]
twins_filtered_dt[, sum_tumour_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_filtered_dt[, sum_tumour_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]
twins_filtered_dt[, sum_normal_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_filtered_dt[, sum_normal_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

twins_filtered_dt[, sum_all_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_all_mtr', 'sum_all_vaf')]
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

twins_filtered_dt[, agg_normal_PD63383_clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_clean_mtr]
twins_filtered_dt[, agg_normal_PD63383_clean_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_clean_dep]
twins_filtered_dt[, agg_normal_PD63383_clean_vaf := agg_normal_PD63383_clean_mtr / agg_normal_PD63383_clean_dep]

# Add a section for clean normal PD62341 samples (exclude spleen which we know contains PD63383 cells)
samples_normal_PD62341_clean = c("PD62341n", "PD62341h", "PD62341ad", "PD62341aa", "PD62341q")
samples_normal_PD62341_clean_dep = paste0(samples_normal_PD62341_clean, '_DEP')
samples_normal_PD62341_clean_mtr = paste0(samples_normal_PD62341_clean, '_MTR')
samples_normal_PD62341_clean_vaf = paste0(samples_normal_PD62341_clean, '_VAF')

twins_filtered_dt[, sum_normal_PD62341_clean_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_clean_mtr]
twins_filtered_dt[, sum_normal_PD62341_clean_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_clean_vaf]
twins_filtered_dt[, sum_normal_PD62341_clean_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD62341_clean_mtr', 'sum_normal_PD62341_clean_vaf')]

twins_filtered_dt[, agg_normal_PD62341_clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_clean_mtr]
twins_filtered_dt[, agg_normal_PD62341_clean_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_clean_dep]
twins_filtered_dt[, agg_normal_PD62341_clean_vaf := agg_normal_PD62341_clean_mtr / agg_normal_PD62341_clean_dep]

# Add a section to identify contaminated normal samples 
samples_normal_contaminated = c("PD62341aa", "PD62341h", "PD63383bb")
samples_normal_contaminated_dep = paste0(samples_normal_contaminated, '_DEP')
samples_normal_contaminated_mtr = paste0(samples_normal_contaminated, '_MTR')
samples_normal_contaminated_vaf = paste0(samples_normal_contaminated, '_VAF')

twins_filtered_dt[, sum_normal_contaminated_mtr := rowSums(.SD>=4), .SDcols = samples_normal_contaminated_mtr]
twins_filtered_dt[, sum_normal_contaminated_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_contaminated_vaf]
twins_filtered_dt[, sum_normal_contaminated_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_contaminated_mtr', 'sum_normal_contaminated_vaf')]

twins_filtered_dt[, agg_normal_contaminated_mtr := rowSums(.SD), .SDcols = samples_normal_contaminated_mtr]
twins_filtered_dt[, agg_normal_contaminated_dep := rowSums(.SD), .SDcols = samples_normal_contaminated_dep]
twins_filtered_dt[, agg_normal_contaminated_vaf := agg_normal_contaminated_mtr / agg_normal_contaminated_dep]

###################################################################################################################################
# Processing the AlleleIntegrator file 

# Number of cells in which a read spanning the mutation position was identified
paste('Number of cells with >= 1 read spanning mutation position:', length(ai_counts_dt[, barcode] %>% unlist() %>% unique())) # 692
paste('Number of reads spanning mutation position (across all cells):', sum(ai_counts_dt[, Tot])) # 748
paste('Number of cells with reads spanning > 1 mutant positions:', sum(duplicated(ai_counts_dt[, barcode]))) # 34 

# Histogram of reads
pdf('Figures/F7/20241208_genotyping_hist_Tot.pdf', height = 3, width = 4)
hist(ai_counts_dt[, Tot], xlab = 'Number of reads over a position per cell',
     main = 'Read number per cell') # only 1 or 2 reads (in 25 cases), 4 reads observed in 1 case
dev.off()

# check the number of reads makes sense given that 
table(ai_counts_dt[, Tot]) # 706 + 38 + 4 = 745

# Add source of the sample (ie, tumour from which twin)
ai_counts_dt[, twin := factor(fcase(
  sample_ID %in% c('NB13652544', 'NB13652545', 'NB13652546'), 'PD62341',
  sample_ID %in% c('NB13760628', 'NB13760629', 'NB13760630'), 'PD63383'
))]

table(ai_counts_dt[, sample_ID]) 
# identified cells from all samples, more in PD62341 than PD63383 (> likely more cells); roughly even between lanes 
table(ai_counts_dt[, twin])
# 486 barcodes from twin 1, 240 from twin 2 (726)

# Add mutation identity and features 

# identify mutation ID 
ai_counts_dt[, Chrom := paste0('chr', chr)] # set the same chromosome name as mut_ID 
ai_counts_dt[, mut_ID0 := paste0(Chrom, '_', pos)]
muts_assignment[, mut_ID0 := substr(mut_ID, 1, nchar(mut_ID)-4)]
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
paste('Number of mutant positions with reads in scRNA-seq:', length(ai_counts_dt[, mut_ID] %>% unlist() %>% unique())) #88 (out of 255)
paste('Number of reads reporting wt variant:', sum(ai_counts_dt[status=='wt', Tot])) # 723 
paste('Number of reads reporting mutation:', sum(ai_counts_dt[status=='mutant', Tot])) # 25 
  
# is each mutation seen only once across all cells, or multiple times?
mut_counts = data.table(table(ai_counts_dt[, mut_ID]))
pdf('Figures/F7/20241208_genotyping_hist_mutID.pdf', height = 5, width = 5)
hist(mut_counts[,N], xlab = 'Number of times mutation seen', main = 'Mutation occurrence', breaks = 40)
dev.off()

# which mutation is seen > 100 times?
muts_10reads = mut_counts[N>10, V1] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_10reads, c('mut_ID', 'Gene', 'Effect'), with=FALSE]
# I checked positions which are not annotated to a gene in the CaVEMan pileup output on Jbrowse (tracks > protein-coding):
# chr10_31333522_G_A - ZEB1 (non-coding region)
# chr1_51714096_G_C - OSBPL9 (oxysterol binding protein like 9, very close to a coding region)
# chr5_54017866_G_T - ARL15 (non-coding region)
# chr7_73831920_C_T - CLDN4 (blue??)
# chrX_70352624_G_A - KIF4A (this seems to be in an exon, too?)

# are those mutations reporting wt or mutant allele?
mut_counts_status = data.table(table(ai_counts_dt[, c('mut_ID', 'status')]))
mut_counts_status[N>10]
mut_counts_status_wide = reshape(mut_counts_status, idvar = c("mut_ID"), timevar = "status", direction = "wide")
mut_counts_status_wide[, dep := N.mutant + N.wt]
mut_counts_status_wide[, vaf := N.mutant / (N.mutant + N.wt)]
mut_counts_status_wide[ dep > 10]

# merge with assignment what type of mutation this is 
mut_counts_status_wide = merge(mut_counts_status_wide, muts_assignment[, c('mut_ID', 'mut_class')], by = 'mut_ID')
mut_counts_status_wide = merge(mut_counts_status_wide, twins_filtered_dt[, c('mut_ID','agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf','vaf_all_tumour_PD62341', 'vaf_all_tumour_PD63383')], by = 'mut_ID')

# which sample were the reads from (PD62341 tumour or PD63383 tumour?)
mut_counts_status_wide = merge(mut_counts_status_wide, ai_counts_dt[, c('mut_ID', 'sample_ID', 'twin')], by = 'mut_ID')
mut_counts_status_wide[!duplicated(mut_counts_status_wide)]

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
table(ai_counts_dt[status == 'mutant', mut_class]) # mostly tumour but you get embryonic ones!

###################################################################################################################################
# Integrate analysis with clusters

# in which clusters are there cells with those reads? what about cells with mutant reads specifically?
meta_PD62341 = tumour_PD62341_clusters@meta.data
meta_PD63383 = tumour_PD63383_clusters@meta.data

# integrate with clusters
ai_counts_dt[, barcodeID :=tstrsplit(barcode, '-', fixed=TRUE, keep=1) ]
ai_counts_dt[, CellID := paste0(barcodeID, '-CG_SB_', sample_ID)]
ai_counts_dt$seurat_clusters_PD62341 = meta_PD62341$seurat_clusters[match(ai_counts_dt$CellID,meta_PD62341$CellID)]
ai_counts_dt$seurat_clusters_PD63383 = meta_PD63383$seurat_clusters[match(ai_counts_dt$CellID,meta_PD63383$CellID)]

# why are some cells apparently not captured? are they removed through QC?
# this affects 226/755 cells so weird if all removed via QC 

# EXAMPLE
# ATCCGAATCCAGTAGT-CG_SB_NB13652545 is the first barcode not assigned to either PD62341 / PD63383 clusters 
tumourPD62341.2 = Read10X(data.dir = "Data/scRNAseq/CG_SB_NB13652545/filtered_feature_bc_matrix/") # lane 2
tPD62341.2 = CreateSeuratObject(counts = tumourPD62341.2, project = "twins.sarcoma", names.field = 1, names.delim = "_", meta.data = NULL)
tPD62341.2@meta.data$CellID = rownames(tPD62341.2@meta.data) 
"ATCCGAATCCAGTAGT-1" %in% tPD62341.2@meta.data$CellID # FALSE 
# so why does my scRNA-seq AlleleIntegrator output say otherwise?

# add information on variant presence in a given cluster 
# create a dataframe with cells assigned to PD62341 tumour clusters, report no read / wt read / mutant read
PD62341_reads = ai_counts_dt[twin=='PD62341' & !is.na(seurat_clusters_PD62341), c('CellID', 'status'), with=FALSE]
PD63383_reads = ai_counts_dt[twin=='PD63383' & !is.na(seurat_clusters_PD63383), c('CellID', 'status'), with=FALSE]

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
  scale_color_manual(values = myCols)+
  theme_classic(base_size = 13)+
  coord_equal(ratio=1)+
  xlim(c(-20, 20))+
  ylim(c(-20, 20))+
  labs(title = 'Reads spanning mutant positions\nPD62341 tumour', x = 'UMAP 1', y = 'UMAP 2', col = 'Read')
ggsave('Figures/F7/20241208_umap_all_mutations_PD62341tumour.pdf', height = 5, width = 5)
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
  scale_color_manual(values = myCols)+
  theme_classic(base_size = 13)+
  coord_equal(ratio=1)+
  xlim(c(-20, 20))+
  ylim(c(-20, 20))+
  labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Tumour mutations positions\nPD62341 tumour', col = 'Read')
ggsave('Figures/F7/20241208_umap_tumour_mutations.pdf', height = 5, width = 5)

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
  scale_color_manual(values = myCols)+
  theme_classic(base_size = 13)+
  coord_equal(ratio=1)+
  xlim(c(-20, 20))+
  ylim(c(-20, 20))+
  labs(x = 'UMAP 1', y = 'UMAP 2', title = 'Normal mutations positions\nPD62341 tumour', col = 'Read')
ggsave('Figures/F7/20241208_umap_normal_mutations.pdf', height = 5, width = 5)
# okay can I check which clusters are those mutations in 

###################################################################################################################################
# Make the plot for each mutation separately

# get all clusters with mean expression of NR5A1 (SF1, tumour marker) > 1
gene_expression_PD62341 = FetchData(tumour_PD62341_clusters, vars = c('NR5A1', "ident"))
tumour_PD62341 = aggregate(gene_expression_PD62341[, 'NR5A1', drop=FALSE  ],
                            by = list(cluster = gene_expression_PD62341$ident), FUN = mean) %>%  
                            filter(NR5A1 > 1) %>% pull(cluster) 

gene_expression_PD63383 = FetchData(tumour_PD63383_clusters, vars = c('NR5A1', "ident"))
tumour_PD63383 = aggregate(gene_expression_PD63383[, 'NR5A1', drop=FALSE ],
                          by = list(cluster = gene_expression_PD63383$ident), FUN = mean) %>%  
                          filter(NR5A1 > 1) %>% pull(cluster)

# list of mutant positions with reads spanning those in scRNA-seq
muts_in_scrna = ai_counts_dt[, mut_ID] %>% unlist() %>% unique()

for (mut in muts_in_scrna){
  
  # get nicely formatted mutation name (for plotting)
  mut_name = gsub("_(.*?)_(.*?)_", ";\\1 \\2>", mut)
  
  # determine which mutation class this mutation was assigned to
  mut_assignment = as.character(muts_assignment[mut_ID == mut, mut_class])
  
  # find which cells (barcodes) report reads with the mutation
  reads = ai_counts_dt[mut_ID == mut]
  reads = reads[!duplicated(reads[,CellID])] # remove duplicate cells
                         
  rownames(reads) = reads$CellID
  
  # only do the plot if the mutant position-spanning reads have been identified in any PD62341 cell
  
  # if reads in PD62341 tumour only 
  if (dim(reads[!is.na(seurat_clusters_PD62341)])[1] > 0 & dim(reads[!is.na(seurat_clusters_PD63383)])[1]==0){
    
    reads = ai_counts_dt[(!is.na(seurat_clusters_PD62341) & is.na(seurat_clusters_PD63383)) & mut_ID == mut]
    reads = reads[!duplicated(reads[,CellID])] # remove duplicate cells
    rownames(reads) = reads$CellID

    # add metadata to clusters 
    tumour_PD62341_clusters = AddMetaData(tumour_PD62341_clusters, metadata = reads)
    
    # show presence in tumour from PD62341 
    umap_PD62341 = FetchData(tumour_PD62341_clusters, vars = c("umap_1", "umap_2", "status", "ident"))
    umap_PD62341 = data.table(umap_PD62341)
    umap_PD62341[, status2 := factor(fcase(
      status == 'mutant', 'mutant read',
      status == 'wt', 'wt read',
      is.na(status) & ident %in% tumour_PD62341, 'tumour cells',
      is.na(status) & !ident %in% tumour_PD62341, 'normal cells'))]
    umap_PD62341[, status2 := factor(status2, levels = c('tumour cells', 'normal cells', 'wt read', 'mutant read'))]
    ggplot(setorder(umap_PD62341, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
      geom_point()+
      scale_color_manual(values = myCols)+
      scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+
      theme_classic(base_size = 13)+
      coord_equal(ratio=1)+
      xlim(c(-20, 20))+
      ylim(c(-20, 20))+ 
      guides(size = "none")+
      labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD62341 tumour\n{mut_name}'), col = 'category')
    ggsave(glue('Figures/F7/20241208_umap_tumour_PD62341only_mutation_{mut}.pdf'), height = 4, width = 4)
    
    ggplot(setorder(umap_PD62341, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
      geom_point()+
      scale_color_manual(values = myCols,  )+
      scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+
      theme_classic(base_size = 13)+
      coord_equal(ratio=1)+
      xlim(c(-20, 20))+
      ylim(c(-20, 20))+ 
      guides(size = "none")+
      labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD62341 tumour\n{mut}\n{mut_assignment}'), col = 'category')
    ggsave(glue('Figures/F7/20241208_umap_tumour_PD62341only_mutation_{mut}_desc.pdf'), height = 4, width = 4)

  }

  # if reads in PD63383 tumour only
  if (dim(reads[!is.na(seurat_clusters_PD62341)])[1] == 0 & dim(reads[!is.na(seurat_clusters_PD63383)])[1] > 0){
    
    reads = ai_counts_dt[(is.na(seurat_clusters_PD62341) & !is.na(seurat_clusters_PD63383)) & mut_ID == mut]
    reads = reads[!duplicated(reads[,CellID])] # remove duplicate cells
    rownames(reads) = reads$CellID
    
    # add metadata to clusters 
    tumour_PD63383_clusters = AddMetaData(tumour_PD63383_clusters, metadata = reads)
    
    # show presence in tumour from PD62341 
    umap_PD63383 = FetchData(tumour_PD63383_clusters, vars = c("umap_1", "umap_2", "status", "ident"))
    umap_PD63383 = data.table(umap_PD63383)
    umap_PD63383[, status2 := factor(fcase(
      status == 'mutant', 'mutant read',
      status == 'wt', 'wt read',
      is.na(status) & ident %in% tumour_PD63383, 'tumour cells',
      is.na(status) & !ident %in% tumour_PD63383, 'normal cells'))]
    umap_PD63383[, status2 := factor(status2, levels = c('tumour cells', 'normal cells', 'wt read', 'mutant read'))]
    
    ggplot(setorder(umap_PD63383, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
      geom_point()+
      scale_color_manual(values = myCols,  )+
      scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+
      theme_classic(base_size = 13)+
      coord_equal(ratio=1)+
      xlim(c(-20, 20))+
      ylim(c(-20, 20))+
      guides(size = "none")+
      labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD63383 tumour\n{mut_name}'), col = 'category')
    ggsave(glue('Figures/F7/20241208_umap_tumour_PD63383only_mutation_{mut}.pdf'), height = 4, width = 4)  
    
    ggplot(setorder(umap_PD63383, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
      geom_point()+
      scale_color_manual(values = myCols,  )+
      scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+
      theme_classic(base_size = 13)+
      coord_equal(ratio=1)+
      xlim(c(-20, 20))+
      ylim(c(-20, 20))+
      guides(size = "none")+
      labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD63383 tumour\n{mut}\n{mut_assignment}'), col = 'category')
    ggsave(glue('Figures/F7/20241208_umap_tumour_PD63383only_mutation_{mut}_desc.pdf'), height = 4, width = 4)  
    
  }
  
  # if reads in both PD62341 and PD63383 
  if (dim(reads[!is.na(seurat_clusters_PD62341)])[1] > 0 & dim(reads[!is.na(seurat_clusters_PD63383)])[1] > 0){
    
    reads_PD62341 = ai_counts_dt[!is.na(seurat_clusters_PD62341) & mut_ID == mut]
    reads_PD63383 = ai_counts_dt[!is.na(seurat_clusters_PD63383) & mut_ID == mut]
    
    reads_PD62341 = reads_PD62341[!duplicated(reads_PD62341[,CellID])] # remove duplicate cells
    reads_PD63383 = reads_PD63383[!duplicated(reads_PD63383[,CellID])]
    
    rownames(reads_PD62341) = reads_PD62341$CellID
    rownames(reads_PD63383) = reads_PD63383$CellID
    
    # add metadata to clusters 
    tumour_PD62341_clusters = AddMetaData(tumour_PD62341_clusters, metadata = reads_PD62341)
    tumour_PD63383_clusters = AddMetaData(tumour_PD63383_clusters, metadata = reads)
      
    # show presence in tumour from PD62341 
    umap_PD62341 = FetchData(tumour_PD62341_clusters, vars = c("umap_1", "umap_2", "status", "ident"))
    umap_PD62341 = data.table(umap_PD62341)
    umap_PD62341[, status2 := factor(fcase(
      status == 'mutant', 'mutant read',
      status == 'wt', 'wt read',
      is.na(status) & ident %in% tumour_PD62341, 'tumour cells',
      is.na(status) & !ident %in% tumour_PD62341, 'normal cells'))]
    umap_PD62341[, status2 := factor(status2, levels = c('tumour cells', 'normal cells', 'wt read', 'mutant read'))]
    
    umap_PD63383 = FetchData(tumour_PD63383_clusters, vars = c("umap_1", "umap_2", "status", "ident"))
    umap_PD63383 = data.table(umap_PD63383)
    umap_PD63383[, status2 := factor(fcase(
      status == 'mutant', 'mutant read',
      status == 'wt', 'wt read',
      is.na(status) & ident %in% tumour_PD63383, 'tumour cells',
      is.na(status) & !ident %in% tumour_PD63383, 'normal cells'))]
    umap_PD63383[, status2 := factor(status2, levels = c('tumour cells', 'normal cells', 'wt read', 'mutant read'))]
    
    # plots with no description, nicely formatted name
    p_PD62341 = ggplot(setorder(umap_PD62341, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
      geom_point()+
      scale_color_manual(values = myCols,  )+
      scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+
      theme_classic(base_size = 13)+
      coord_equal(ratio=1)+
      xlim(c(-20, 20))+
      ylim(c(-20, 20))+ 
      guides(size = "none")+
      labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD62341 tumour\n{mut_name}'), col = 'category')
    
    p_PD63383 = ggplot(setorder(umap_PD63383, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
      geom_point()+
      scale_color_manual(values = myCols,  )+
      scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+
      theme_classic(base_size = 13)+
      coord_equal(ratio=1)+
      xlim(c(-20, 20))+
      ylim(c(-20, 20))+  
      guides(size = "none")+
      labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD63383 tumour\n{mut_name}'), col = 'category')
    
    plots = grid.arrange(p_PD62341, p_PD63383, ncol = 2)
    ggsave(glue('Figures/F7/20241208_umap_tumour_both_mutation_{mut}.pdf'), plots, height = 4, width = 8)
    
    # plots with description of the class of the mutation
    p_PD62341 = ggplot(setorder(umap_PD62341, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
        geom_point()+
        scale_color_manual(values = myCols)+
        scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+
        theme_classic(base_size = 13)+
        coord_equal(ratio=1)+
        xlim(c(-20, 20))+
        ylim(c(-20, 20))+
        guides(size = "none")+
        labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD62341 tumour\n{mut}\n{mut_assignment}'), col = 'category')
    
    p_PD63383 = ggplot(setorder(umap_PD63383, status2), aes(umap_1, umap_2, color = status2, order = status2, size = status2))+
        geom_point()+
        scale_color_manual(values = myCols)+
        scale_size_manual(values = c('tumour cells' = 0.1, 'normal cells' = 0.1, 'wt read' = 0.5, 'mutant read' = 0.5))+  
        theme_classic(base_size = 13)+
        coord_equal(ratio=1)+
        xlim(c(-20, 20))+
        ylim(c(-20, 20))+ 
        guides(size = "none")+
        labs(x = 'UMAP 1', y = 'UMAP 2', title = glue('PD63383 tumour\n{mut}\n{mut_assignment}'), col = 'category')
      
    plots = grid.arrange(p_PD62341, p_PD63383, ncol = 2)
    ggsave(glue('Figures/F7/20241208_umap_tumour_both_mutation_{mut}_desc.pdf'), plots, height = 4, width = 8)  
    
  }

}

###################################################################################################################################
# ALL DONE!
