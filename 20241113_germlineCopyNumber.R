###################################################################################################################################
# SCRIPT 2

# Script to remove germline mutations in copy number-altered regions
# 2024-11-13
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters (from script: 20241113_pileup_filters.R)

# OUTPUT:
# 1 list of mutations which pass required filters AND are not on copy number altered regions
# these mutations will be used for phylogeny reconstruction

# Script to identify regions of copy number alterations (deletions or duplications) in the germline
# These regions can be small (several 100 kb) and would not be picked up by ASCAT 

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

###################################################################################################################################
# INPUT FILES 

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/twins_dt_20241114_1134.csv')) # import high quality pileup with added filters

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Create a dataframe that only includes mutations retained post filtering  
muts_included = read.table('Data/mutations_include_20241114_1134.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts_included)) # 1134

muts_germline = read.table('Data/mutations_putativeGermline_20241114.txt') %>% unlist()
paste('Number of mutations identified as putative germline:', length(muts_germline)) # 333,031

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

col_background = '#c8c8c8'
col_germline = '#d3b28b'
col_other = '#25b6db'
col_normal_all = '#5265c9'
col_normal_all_removed = '#b9dff5'
col_clusters = '#bf068a'
col_coverage = '#bf7a06'
col_clusters_coverage = '#810cb9'
col_excluded = '#a80505'

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
# AGGREGATED MTR, DEP and VAF (by sample category) 

# Aggregated VAF in tumour samples 
twins_dt[, dep_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_dt[, dep_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_dt[, dep_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_dt[, dep_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_dt[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]
twins_dt[, dep_all_tumour := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_dt[, dep_all := rowSums(.SD), .SDcols = samples_dep]

twins_dt[, mtr_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_dt[, mtr_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_dt[, mtr_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_dt[, mtr_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
twins_dt[, mtr_all_normal := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_dt[, mtr_all_tumour := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_dt[, mtr_all := rowSums(.SD), .SDcols = samples_mtr]

twins_dt[, vaf_all_normal_PD62341 := mtr_all_normal_PD62341 / dep_all_normal_PD62341]
twins_dt[, vaf_all_normal_PD63383 := mtr_all_normal_PD63383 / dep_all_normal_PD63383]
twins_dt[, vaf_all_tumour_PD62341 := mtr_all_tumour_PD62341 / dep_all_tumour_PD62341]
twins_dt[, vaf_all_tumour_PD63383 := mtr_all_tumour_PD63383 / dep_all_tumour_PD63383]
twins_dt[, vaf_all_normal := mtr_all_normal / dep_all_normal]
twins_dt[, vaf_all_tumour := mtr_all_tumour / dep_all_tumour]
twins_dt[, vaf_all := mtr_all / dep_all]

######################################################################################################
# MTR, DEP, VAF DATA 

# Create dataframe with mutations that passed current filters 
twins_filtered_dt = twins_dt[mut_ID %in% muts_included]

# Subset MTR, DEP and VAF data
twins_filtered_mtr = twins_filtered_dt[,c('mut_ID', samples_mtr), with=FALSE]
twins_filtered_dep = twins_filtered_dt[,c('mut_ID', samples_dep), with=FALSE]
twins_filtered_vaf = twins_filtered_dt[,c('mut_ID', samples_vaf), with=FALSE]

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
# Plot coverage for all identified mutations across the genome (removing mutations that failed quality filters)

# Identify mutations which are poor quality to exclude from plotting (allow germline mutations of good quality)
columns_req_filters_qual = c('f1_mappedY', 'f2_FailedIndelNearby30','f3_lowDepthNormal', 'f3_highDepthNormal', 
                             'f4_mtrAndVaf', 'f5_strandBiasMutOnly', 'f6_lowQualRatio')
twins_dt[, sum_req_filters_qual := rowSums(.SD), .SDcols = columns_req_filters_qual] 

# Add data on chromosome lengths 
Chrom = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
          "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
          "chr20", "chr21", "chr22", "chrX")
# Source: https://www.ncbi.nlm.nih.gov/grc/human/data (Hg38), accessed 14/11/2024
Chrom_length = as.numeric(c(248956422, 242193529, 198295559, 190214555,
                            181538259, 170805979, 159345973, 145138636,
                            138394717, 133797422, 135086622, 133275309,
                            114364328, 107043718, 101991189, 90338345,
                            83257441, 80373285, 58617616, 64444167,
                            46709983, 50818468, 156040895))
chr_lengths_dt = data.table(cbind(Chrom, Chrom_length))

twins_dt = merge(twins_dt, chr_lengths_dt, by = 'Chrom')
twins_dt[, Chrom := factor(Chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                             "chr20", "chr21", "chr22", "chrX"))]
twins_dt[, pos := tstrsplit(mut_ID, '_', fixed=T, keep=2)]
twins_dt[, pos := as.numeric(pos)]

for (chr in Chrom){
  dt = twins_dt[Chrom == chr & sum_req_filters_qual==0] # exclude mutations removed due to poor quality 
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = dep_all_normal))+
    geom_point(size=1.5, alpha = 0.6, col = col_background)+
    theme_classic(base_size = 12)+
    labs(x = 'Genomic position', y = 'Total coverage (normal samples)')+
    ggtitle(glue('{chr}, {nr} mutations'))+
    xlim(c(0, length))
  ggsave(glue('Results/20241114_p2_dep_across_genome_{chr}_allmuts.pdf'), height=3, width=6.5)
}
# can see clusters of markedly higher coverage across some chromosomes 

######################################################################################################
# Identify mutations present in all normal samples yet not identified as germline 

muts_normal_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist())) 
paste('Number of mutations identified in all normal samples:', length(muts_normal_all)) # 632 

# Determine the VAF of mutations present in all normal samples: why were these not called as germline?
# As seen on the plots below, the VAF of the vast majority of those mutations is clearly below 0.5

pdf('Results/20241114_p2_hist_vaf_632muts_in_all_normal_samples.pdf', width = 4.2, height = 3.2)
hist(twins_dt[mut_ID %in% muts_normal_all, vaf_all_normal] %>% unlist(), 
     xlab = 'VAF (agg normal samples)', main = 'Mutations in all normal samples (632)', xlim = c(0, 1), breaks=30)
abline(v = median(twins_dt[mut_ID %in% muts_normal_all, vaf_all_normal]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241114_p2_hist_vaf_632muts_in_all_normal_samples_PD62341.pdf', width = 4.2, height = 3.2)
hist(twins_dt[mut_ID %in% muts_normal_all, vaf_all_normal_PD62341] %>% unlist(), 
     xlab = 'VAF (agg normal PD62341)', main = 'Mutations in all normal samples (632)', xlim = c(0, 1), breaks=30)
abline(v = median(twins_dt[mut_ID %in% muts_normal_all, vaf_all_normal_PD62341]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241114_p2_hist_vaf_632muts_in_all_normal_samples_PD63383.pdf', width = 4.2, height = 3.2)
hist(twins_dt[mut_ID %in% muts_normal_all, vaf_all_normal_PD63383] %>% unlist(), 
     xlab = 'VAF (agg normal PD63383)', main = 'Mutations in all normal samples (632)', xlim = c(0, 1), breaks=30)
abline(v = median(twins_dt[mut_ID %in% muts_normal_all, vaf_all_normal_PD63383]), col = 'purple', lwd = 2.5)
dev.off()

######################################################################################################
# Plots: VAF vs DEP in high-quality mutations in normal and tumour samples 
# This allows to identify clusters of mutations: germline, germline on copy number changes, somatic (VAF < 0.5)

# Indicate classes of mutations of interest 
twins_dt[, mut_cat := factor(fcase(
  f7_likelyGermline_bothTwins == 1, 'germline',
  mut_ID %in% muts_included == 1 & !mut_ID %in% muts_normal_all, 'included, other',
  mut_ID %in% muts_included == 1 & mut_ID %in% muts_normal_all, 'included, present in all samples' # mutations that could be germline
))]

# there are ~333k germline mutations: to avoid over-plotting, sample 5k germline mutations 
set.seed(10)
muts_germline_sample = sample(muts_germline, 5000, replace = FALSE)

ggplot(twins_dt[sum_req_filters_qual==0 & mut_ID %in% c(muts_germline_sample, muts_included)], aes(x = dep_all, y = vaf_all, col = mut_cat))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (all samples)', y = 'VAF (all samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_germline, col_other, col_normal_all))+
  ggtitle(glue('All samples'))+
  ylim(c(0, 0.75))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_allMuts_all_small.pdf'), height=4, width=6)

ggplot(twins_dt[sum_req_filters_qual==0 & mut_ID %in% c(muts_germline_sample, muts_included)], aes(x = dep_all_normal, y = vaf_all_normal, col = mut_cat))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (normal samples)', y = 'VAF (normal samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_germline, col_other, col_normal_all))+
  ggtitle(glue('Normal samples'))+
  ylim(c(0, 0.75))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_allMuts_normal_small.pdf'), height=4, width=6)

ggplot(twins_dt[sum_req_filters_qual==0 & mut_ID %in% c(muts_germline_sample, muts_included)], aes(x = dep_all_tumour, y = vaf_all_tumour, col = mut_cat))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (tumour samples)', y = 'VAF (tumour samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_germline, col_other, col_normal_all))+
  ggtitle(glue('Tumour samples'))+
  ylim(c(0, 0.75))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_allMuts_tumour_small.pdf'), height=4, width=6)

# Plot VAF vs DEP for all included mutations, indicating all normal vs other mutations 
ggplot(twins_dt[mut_ID %in% muts_included], aes(x = dep_all, y = vaf_all, col = mut_cat))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (all samples)', y = 'VAF (all samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_other, col_normal_all))+
  ggtitle(glue('All samples'))+
  ylim(c(0, 0.6))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_all.pdf'), height=3.5, width=5.5)

ggplot(twins_dt[mut_ID %in% muts_included], aes(x = dep_all_normal, y = vaf_all_normal, col = mut_cat))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (normal samples)', y = 'VAF (normal samples)', col = 'Mutation category')+
  ggtitle(glue('Normal samples'))+
  scale_color_manual(values =c(col_other, col_normal_all))+
  ylim(c(0, 0.6))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_normal.pdf'), height=3.5, width=5.5)

ggplot(twins_dt[mut_ID %in% muts_included], aes(x = dep_all_tumour, y = vaf_all_tumour, col = mut_cat))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (tumour samples)', y = 'VAF (tumour samples)', col = 'Mutation category')+
  ggtitle(glue('Tumour samples'))+
  scale_color_manual(values = c(col_other, col_normal_all))+
  ylim(c(0, 0.6))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_tumour.pdf'), height=3.5, width=5.5)

######################################################################################################
# Possible explanation: mutations present in all samples are germline mutations present in regions of copy number changes 

# 1 plot the coverage of mutations present in all samples against other mutations across the genome

# Add column to indicate mutation category (ignore germline)
twins_dt[, mut_cat2 := factor(fcase(
  mut_ID %in% muts_normal_all, 'mutations present in\nall normal samples',
  !mut_ID %in% muts_normal_all, 'all other mutations'))]

# Specify plotting colors 
my_colors = c(col_background, col_normal)
names(my_colors) = levels(twins_dt[, mut_cat2])

# compare depth across mutations present in all normal samples at lower VAF (< 0.5) to all other mutations in the dataset  
for (chr in Chrom){
  dt = twins_dt[sum_req_filters_qual== 0 & Chrom == chr] # don't plot mutations that fail QC filters 
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(setorder(dt, mut_cat2), aes(x = pos, y = dep_all, col = mut_cat2, alpha = mut_cat2, order = mut_cat2))+
    geom_point(size=1.5)+
    guides(alpha = "none")+ # remove legend for alpha 
    theme_classic(base_size = 12)+
    scale_color_manual(values = c(my_colors))+
    scale_alpha_discrete(range = c(0.1, 0.8))+
    labs(x = 'Genomic position', y = 'Total coverage (all samples)', col = 'Mutation class')+
    ggtitle(glue('{chr}'))+
    xlim(c(0, length))
  ggsave(glue('Results/20241114_p2_dep_all_{chr}_col_mutsAll.pdf'), height=3, width=6.5)
}

for (chr in Chrom){
  dt = twins_dt[sum_req_filters_qual== 0 & Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(setorder(dt, mut_cat2), aes(x = pos, y = dep_all_normal, col = mut_cat2, alpha = mut_cat2, order = mut_cat))+
    geom_point(size=1.5)+
    guides(alpha = "none")+ # remove legend for alpha 
    theme_classic(base_size = 12)+
    scale_color_manual(values = c(my_colors))+
    scale_alpha_discrete(range = c(0.1, 0.8))+
    labs(x = 'Genomic position', y = 'Total coverage (normal samples)', col = 'Mutation class')+
    ggtitle(glue('{chr}'))+
    xlim(c(0, length))+
    ylim(c(0, 900))
  ggsave(glue('Results/20241114_p2_dep_normal_{chr}_col_mutsAllNormal.pdf'), height=3, width=6.5)
}

# for each chromosome, calculate the median depth of all mutations present on the chromosome
chrom_median_depths = twins_dt[sum_req_filters_qual==0, lapply(.SD, function(x) median(x)), .SDcols = 'dep_all_normal', by=Chrom]
setnames(chrom_median_depths, c('Chrom', 'dep_all_normal'), c('Chrom', 'median_dep_chr'))
twins_dt = merge(twins_dt, chrom_median_depths, by = 'Chrom')

# Note: if duplication, then increase from 2 to 3 copies: increase the number of reads 1.5x
# Note: if deletion, then decrease from 2 to 1 copies: increase the number of reads 0.5x
# Therefore, I am setting thresholds lower than this 
twins_dt[, upper_dep := 1.25 * median_dep_chr] 
twins_dt[, lower_dep := 0.75 * median_dep_chr] 

# plot to indicate median depth and threshold depth values to identify mutations to exclude  
for (chr in Chrom){
  dt = twins_dt[sum_req_filters_qual == 0 & Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = dep_all_normal, col = mut_cat2, alpha = mut_cat2, order = mut_cat2))+
      geom_point(size=1.5)+
      guides(alpha = "none")+ # remove legend for alpha 
      theme_classic(base_size = 15)+
      scale_color_manual(values = c(my_colors))+
      scale_alpha_discrete(range = c(0.8, 0.1))+
      labs(x = 'Genomic position', y = 'Coverage (normal samples)', col = 'Mutation class')+
      ggtitle(glue('{chr}'))+
      xlim(c(0, length))+
      ylim(c(0, 850))+
      guides(col="none")+
      theme(axis.text.x = element_text(size=5))+
      geom_hline(yintercept = as.numeric(dt[,median_dep_chr] %>% unique()), col = 'purple', linetype = 'dashed', size = 0.6, alpha = 0.8)+
      geom_hline(yintercept = as.numeric(dt[,upper_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)+
      geom_hline(yintercept = as.numeric(dt[,lower_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)
    ggsave(glue('Results/20241114_p2_dep_normal_{chr}_col_mutsAllNormal_thresholds.pdf'), height=3, width=5.5)
}

ggplot(dt, aes(x = pos, y = dep_all_normal, col = mut_cat2, alpha = mut_cat2, order = mut_cat2))+
  geom_point(size=1.5)+
  guides(alpha = "none")+ # remove legend for alpha 
  theme_classic(base_size = 12)+
  scale_color_manual(values = c(my_colors))+
  scale_alpha_discrete(range = c(0.8, 0.1))+
  labs(x = 'Genomic position', y = 'Total coverage (normal samples)', col = 'Mutation class')+
  ggtitle(glue('{chr}'))+
  xlim(c(0, length))+
  ylim(c(0, 850))+
  geom_hline(yintercept = as.numeric(dt[,median_dep_chr] %>% unique()), col = 'purple', linetype = 'dashed', size = 0.6, alpha = 0.8)+
  geom_hline(yintercept = as.numeric(dt[,upper_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)+
  geom_hline(yintercept = as.numeric(dt[,lower_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)
ggsave(glue('Results/20241114_p2_dep_normal_{chr}_col_mutsAllNormal_thresholds_legend.pdf'), height=3, width=6.5)

# Identify mutations which are above or below the threshold 
twins_dt[, CovThreshold := as.numeric(dep_all_normal < lower_dep | dep_all_normal > upper_dep)]
paste('Number of mutations with unusual coverage compared to the chromosome:', dim(twins_dt[CovThreshold==1])[1]) # 32,835
paste('Fraction of mutations with unusual coverage compared to the chromosome:', round(dim(twins_dt[CovThreshold==1])[1] / dim(twins_dt)[1], 4)) # 0.0906

# Determine the number of mutations present in all samples that are excluded by this criterion
paste('Number of mutations with unusual coverage found in all normal samples at VAF < 0.5:', dim(twins_dt[mut_ID %in% muts_normal_all & CovThreshold==1])[1]) # 476
paste('Fraction of mutations with unusual coverage found in all normal samples at VAF < 0.5:', round(dim(twins_dt[mut_ID %in% muts_normal_all & CovThreshold==1])[1] / length(muts_normal_all), 4)) # 0.7532
# those mutations can be excluded based on the threshold  

twins_dt[, f8_excludeCovThreshold := as.numeric((dep_all_normal < lower_dep | dep_all_normal > upper_dep) & mut_ID %in% muts_normal_all)]
paste('Number of mutations with unusual coverage compared to the chromosome (to exclude):', dim(twins_dt[f8_excludeCovThreshold==1])[1]) # 32,835

# Show mutations excluded by filtering 
twins_dt[, mut_cat3 := factor(fcase(
  f8_excludeCovThreshold == 1, 'all normal samples, excluded',
  mut_ID %in% muts_normal_all & f8_excludeCovThreshold == 0, 'all normal samples, retained',
  !mut_ID %in% muts_normal_all & f7_likelyGermline_bothTwins == 0, 'other', 
  f7_likelyGermline_bothTwins == 1, 'germline'
))]

ggplot(twins_dt[mut_ID %in% c(muts_included, muts_germline_sample)], aes(x = dep_all, y = vaf_all, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (all samples)', y = 'VAF (all samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_normal_all_removed, col_normal_all, col_germline, col_other))+
  ggtitle(glue('All samples'))+
  ylim(c(0, 0.75))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_all_excludeCNGermline_show.pdf'), height=4, width=6)

ggplot(twins_dt[mut_ID %in% c(muts_included, muts_germline_sample)], aes(x = dep_all_normal, y = vaf_all_normal, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (normal samples)', y = 'VAF (normal samples)', col = 'Mutation category')+
  ggtitle(glue('Normal samples'))+
  scale_color_manual(values = c(col_normal_all_removed, col_normal_all, col_germline, col_other))+
  ylim(c(0, 0.75))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_normal_excludeCNGermline_show.pdf'), height=4, width=6)

ggplot(twins_dt[mut_ID %in% c(muts_included, muts_germline_sample)], aes(x = dep_all_tumour, y = vaf_all_tumour, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (tumour samples)', y = 'VAF (tumour samples)', col = 'Mutation category')+
  ggtitle(glue('Tumour samples'))+
  scale_color_manual(values = c(col_normal_all_removed, col_normal_all, col_germline, col_other))+
  ylim(c(0, 0.75))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_tumour_excludeCNGermline_show.pdf'), height=4, width=6)

ggplot(twins_dt[f8_excludeCovThreshold == 0 & mut_ID %in% muts_included], aes(x = dep_all, y = vaf_all, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (all samples)', y = 'VAF (all samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_normal_all, col_other))+
  ggtitle(glue('All samples'))+
  ylim(c(0, 0.6))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_all_excludeCNGermline.pdf'), height=3.5, width=5.5)

ggplot(twins_dt[f8_excludeCovThreshold == 0 &mut_ID %in% muts_included], aes(x = dep_all_normal, y = vaf_all_normal, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (normal samples)', y = 'VAF (normal samples)', col = 'Mutation category')+
  ggtitle(glue('Normal samples'))+
  scale_color_manual(values = c(col_normal_all, col_other))+
  ylim(c(0, 0.6))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_normal_excludeCNGermline.pdf'), height=3.5, width=5.5)

ggplot(twins_dt[f8_excludeCovThreshold == 0 & mut_ID %in% muts_included], aes(x = dep_all_tumour, y = vaf_all_tumour, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (tumour samples)', y = 'VAF (tumour samples)', col = 'Mutation category')+
  ggtitle(glue('Tumour samples'))+
  scale_color_manual(values = c(col_normal_all, col_other))+
  ylim(c(0, 0.6))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_tumour_excludeCNGermline.pdf'), height=3.5, width=5.5)

######################################################################################################
# All mutations in the cluster are likely to be affected by the copy number change
# Therefore, identify and remove all mutations in clusters (even if a given mutation is still within the coverage threshold)
# The threshold is arbitrary - a given mutation may be less well mapped / sequence to a slightly lower coverage and just about meet the threshold

# define a function to find clusters of mutations 
find_clusters = function(numbers, range=20000, min_count = 10){
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

for (chr in twins_dt[, Chrom] %>% unique){
  pos_chr = sort(twins_dt[mut_ID %in% muts_included & Chrom == chr, pos] %>% unlist())
  clusters_chr = find_clusters(pos_chr)
  muts_clusters_chr = twins_dt[mut_ID %in% muts_included & Chrom == chr & pos > min(clusters_chr) & pos < max(clusters_chr), mut_ID] %>% unlist()
  muts_clusters = c(muts_clusters, muts_clusters_chr)
}

muts_clusters = muts_clusters %>% unlist()
paste('Number of mutations presents in clusters:', length(muts_clusters)) # 319

paste('Number of mutations present in clusters with unusually high or low coverage:', length(Reduce(intersect, list(muts_clusters, muts_copynumber_exclude)))) # 229 
paste('Number of mutations present in clusters or with unusually high or low coverage:', length(c(muts_clusters, muts_copynumber_exclude) %>% unique())) # 566 

# Plot the distribution of mutations present in clusters across the genome 
twins_dt[, mut_cat4 := factor(fcase(
  mut_ID %in% muts_clusters & mut_ID %in% muts_copynumber_exclude, 'mutation in clusters with abnormal coverage',
  mut_ID %in% muts_clusters & !mut_ID %in% muts_copynumber_exclude, 'mutation in clusters',
  !mut_ID %in% muts_clusters & mut_ID %in% muts_copynumber_exclude, 'mutation with abnormal coverage',
  !mut_ID %in% muts_clusters & !mut_ID %in% muts_copynumber_exclude, 'other'
))]

my_colors = c(col_clusters, col_clusters_coverage, col_coverage, col_background)
names(my_colors) = levels(twins_dt[, mut_cat4])

for (chr in Chrom){
  dt = twins_dt[sum_req_filters_qual == 0 & Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = dep_all_normal, col = mut_cat4, alpha = mut_cat4, order = mut_cat4))+
    geom_point(size=1.5)+
    guides(alpha = "none")+ # remove legend for alpha 
    theme_classic(base_size = 15)+
    scale_color_manual(values = my_colors)+
    scale_alpha_discrete(range = c(0.8, 0.8, 0.8, 0.1))+
    labs(x = 'Genomic position', y = 'Coverage (normal samples)', col = 'Mutation class')+
    ggtitle(glue('{chr}'))+
    xlim(c(0, length))+
    guides(col="none")+
    ylim(c(0, 850))+
    theme(axis.text.x = element_text(size=5))+
    geom_hline(yintercept = as.numeric(dt[,median_dep_chr] %>% unique()), col = 'purple', linetype = 'dashed', size = 0.6, alpha = 0.8)+
    geom_hline(yintercept = as.numeric(dt[,upper_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)+
    geom_hline(yintercept = as.numeric(dt[,lower_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)
  ggsave(glue('Results/20241114_p2_dep_normal_{chr}_col_mutsAllNormal_clusters.pdf'), height=3, width=5.5)
}

# remove mutations flagged either by coverage or by presence in clusters
twins_dt[, mut_cat5 := factor(fcase(
  mut_ID %in% muts_clusters | mut_ID %in% muts_copynumber_exclude, 'mutation to exclude',
  !mut_ID %in% muts_clusters & !mut_ID %in% muts_copynumber_exclude, 'other'
))]

for (chr in Chrom){
  dt = twins_dt[sum_req_filters_qual == 0 & Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = dep_all_normal, col = mut_cat5, alpha = mut_cat5, order = mut_cat5))+
    geom_point(size=1.5)+
    guides(alpha = "none")+ # remove legend for alpha 
    theme_classic(base_size = 15)+
    scale_color_manual(values = c(col_exclude, col_background))+
    scale_alpha_discrete(range = c(0.8, 0.8, 0.8, 0.1))+
    labs(x = 'Genomic position', y = 'Coverage (normal samples)', col = 'Mutation class')+
    ggtitle(glue('{chr}'))+
    xlim(c(0, length))+
    guides(col="none")+
    ylim(c(0, 850))+
    theme(axis.text.x = element_text(size=5))+
    geom_hline(yintercept = as.numeric(dt[,median_dep_chr] %>% unique()), col = 'purple', linetype = 'dashed', size = 0.6, alpha = 0.8)+
    geom_hline(yintercept = as.numeric(dt[,upper_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)+
    geom_hline(yintercept = as.numeric(dt[,lower_dep] %>% unique()), col = 'black', linetype = 'dashed', size = 0.6, alpha = 0.8)
  ggsave(glue('Results/20241114_p2_dep_normal_{chr}_col_mutsAllNormal_exclude_cov_or_clusters.pdf'), height=3, width=5.5)
}

# Add column on presence in clusters 
twins_dt[, f9_presenceInClusters := as.numeric(mut_ID %in% muts_clusters)]
paste('Number of mutations present in clusters:', length(muts_clusters)) # 319

twins_dt[, mut_cat6 := factor(fcase(
  (mut_ID %in% muts_clusters | mut_ID %in% muts_copynumber_exclude) & mut_ID %in% muts_included, 'present in all normal samples (excluded)',
  (!mut_ID %in% muts_clusters & mut_ID %in% muts_normal_all), 'present in all normal samples (retained)',
  (!mut_ID %in% muts_clusters & !mut_ID %in% muts_copynumber_exclude) & mut_ID %in% muts_included, 'other mutations (retained)',
  mut_ID %in% muts_germline, 'germline', 
  !mut_ID %in% muts_germline & !mut_ID %in% muts_included, 'other'
))]

# Show the set of mutations after removing by coverage or clustering 
ggplot(twins_dt[mut_ID %in% c(muts_included, muts_germline_sample)], aes(x = dep_all, y = vaf_all, col = mut_cat6))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (all samples)', y = 'VAF (all samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_germline, col_other, col_normal_all_removed, col_normal_all))+
  ggtitle(glue('All samples'))+
  ylim(c(0, 0.75))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_all_excludeCNGermlineCluster_show.pdf'), height=4, width=6)

ggplot(twins_dt[mut_ID %in% c(muts_included, muts_germline_sample)], aes(x = dep_all_normal, y = vaf_all_normal, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (normal samples)', y = 'VAF (normal samples)', col = 'Mutation category')+
  ggtitle(glue('Normal samples'))+
  scale_color_manual(values = c(col_germline, col_other, col_normal_all_removed, col_normal_all))+
  ylim(c(0, 0.75))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_normal_excludeCNGermlineCluster_show.pdf'), height=4, width=6)

ggplot(twins_dt[mut_ID %in% c(muts_included, muts_germline_sample)], aes(x = dep_all_tumour, y = vaf_all_tumour, col = mut_cat3))+
  geom_point(size=1.5, alpha = 0.3)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (tumour samples)', y = 'VAF (tumour samples)', col = 'Mutation category')+
  ggtitle(glue('Tumour samples'))+
  scale_color_manual(values = c(col_germline, col_other, col_normal_all_removed, col_normal_all))+
  ylim(c(0, 0.75))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_tumour_excludeCNGermlineCluster_show.pdf'), height=4, width=6)

# Show current set 
ggplot(twins_dt[f8_excludeCovThreshold == 0 & f9_presenceInClusters == 0 & mut_ID %in% muts_included], aes(x = dep_all, y = vaf_all, col = mut_cat6))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (all samples)', y = 'VAF (all samples)', col = 'Mutation category')+
  scale_color_manual(values = c(col_other, col_normal_all))+
  ggtitle(glue('All samples'))+
  ylim(c(0, 0.6))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_all_excludeCNGermlineClusters.pdf'), height=3.5, width=5.5)

ggplot(twins_dt[f8_excludeCovThreshold == 0 & f9_presenceInClusters == 0 & mut_ID %in% muts_included], aes(x = dep_all_normal, y = vaf_all_normal, col = mut_cat6))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (normal samples)', y = 'VAF (normal samples)', col = 'Mutation category')+
  ggtitle(glue('Normal samples'))+
  scale_color_manual(values = c(col_other, col_normal_all))+
  ylim(c(0, 0.6))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_normal_excludeCNGermlineClusters.pdf'), height=3.5, width=5.5)

ggplot(twins_dt[f8_excludeCovThreshold == 0 & f9_presenceInClusters == 0 & mut_ID %in% muts_included], aes(x = dep_all_tumour, y = vaf_all_tumour, col = mut_cat6))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'Coverage (tumour samples)', y = 'VAF (tumour samples)', col = 'Mutation category')+
  ggtitle(glue('Tumour samples'))+
  scale_color_manual(values = c(col_other, col_normal_all))+
  ylim(c(0, 0.6))+
  xlim(c(0, 900))
ggsave(glue('Results/20241114_p2_vaf_vs_dep_1134_tumour_excludeCNGermlineClusters.pdf'), height=3.5, width=5.5)

######################################################################################################
# OUTPUT 1: SAVE FILTERED MUTATIONS TO A TXT FILE
muts_copynumber_exclude = twins_dt[f8_excludeCovThreshold==1, mut_ID] %>% unlist()
muts_retained = setdiff(muts_included, muts_copynumber_exclude)
paste('Number of mutations retained after excluding likely copy number mutations:', length(muts_retained))
write.table(muts_retained, 'Data/mutations_include_20241114_658.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
# OUTPUT 2: SAVE THE DATAFRAME WITH DATA + FILTERS TO A FILE (INCUDING ALL CLASSES OF FILTERS)
write.csv(twins_dt, 'Data/twins_dt_filters_20241114_658_full.csv', quote = FALSE, row.names = F)

######################################################################################################
# OUTPUT 3: SAVE THE DATAFRAME WITH FILTERS TO A FILE (INCUDING ALL CLASSES OF FILTERS)
cols_filters = c('f1_mappedY', 'f2_FailedIndelNearby30','f3_lowDepthNormal', 'f3_highDepthNormal', 
                'f4_mtrAndVaf', 'f5_strandBiasMutOnly', 'f6_lowQualRatio',
                'f7_likelyGermline_bothTwins', 'f8_excludeCovThreshold')
write.csv(twins_dt[, c('mut_ID', cols_filters), with=FALSE], 'Data/twins_dt_filters_20241114_658_filters.csv', quote = FALSE, row.names = F)

######################################################################################################
######################################################################################################

# Added filtering by presence in clusters 

######################################################################################################
# OUTPUT 1: SAVE FILTERED MUTATIONS TO A TXT FILE
muts_copynumber_clusters_exclude = twins_dt[f8_excludeCovThreshold==1 | f9_presenceInClusters==1, mut_ID] %>% unlist()
muts_retained2 = setdiff(muts_included, muts_copynumber_clusters_exclude)
paste('Number of mutations retained after excluding likely copy number mutations and mutations in clusters:', length(muts_retained2)) # 568
write.table(muts_retained2, 'Data/mutations_include_20241114_568.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
# OUTPUT 2: SAVE THE DATAFRAME WITH DATA + FILTERS TO A FILE (INCUDING ALL CLASSES OF FILTERS)
write.csv(twins_dt, 'Data/twins_dt_filters_20241114_568_full.csv', quote = FALSE, row.names = F)

######################################################################################################
# OUTPUT 3: SAVE THE DATAFRAME WITH FILTERS TO A FILE (INCUDING ALL CLASSES OF FILTERS)
cols_filters = c('f1_mappedY', 'f2_FailedIndelNearby30','f3_lowDepthNormal', 'f3_highDepthNormal', 
                 'f4_mtrAndVaf', 'f5_strandBiasMutOnly', 'f6_lowQualRatio',
                 'f7_likelyGermline_bothTwins', 'f8_excludeCovThreshold', 'f9_presenceInClusters')
write.csv(twins_dt[, c('mut_ID', cols_filters), with=FALSE], 'Data/twins_dt_filters_20241114_568_filters.csv', quote = FALSE, row.names = F)

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# ALL DONE  






