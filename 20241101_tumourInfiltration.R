###################################################################################################################################
# SCRIPT 2

# Script to estimate purity of each sample (normal + tumour cell content in each sample)
# November 2024
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters (255, 08/12/2024)

# OUTPUT: 
# 1 list of mutations truncal in the tumour (clonal in the tumour, absent in normal)
# 2 dataframe with estimated normal and tumour cell fraction for each analysed sample 

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

######################################################################################################
# TUMOUR INFILTRATION 

# For true tumour mutations in normal samples, we expect two clusters of mutations present in the tumour
# 1 mutations shared between normal samples and tumour samples (embryonic mosaic)
# 2 mutations private to the tumour and absent from normal samples (expected VAF in normal = 0)

# 1 select mutations which are expected to be clonal in the tumour 
pdf('Figures/F2/20241208_hist_muts255_mean_median.pdf', width = 4, height = 4)
hist(twins_filtered_dt[, vaf_all_tumour], breaks=30, xlab = 'VAF', xlim=c(0,1), main = 'VAF in aggregated tumour samples\n(255 mutations)')
abline(v = median(twins_filtered_dt[, vaf_all_tumour]), col = 'blue', lwd = 2)
abline(v = mean(twins_filtered_dt[, vaf_all_tumour]), col = 'red', lwd = 2)
dev.off()

# indicate the threshold above which mutations can be considered clonal
pdf('Figures/F2/20241208_hist_muts_tumourThreshold.pdf', width = 4, height = 4)
hist(twins_filtered_dt[, vaf_all_tumour], breaks=30, xlab = 'VAF', xlim=c(0,1), main = 'VAF in aggregated tumour samples\n(255 mutations)')
abline(v = 0.2, col = 'purple', lwd = 3)
dev.off()

# require VAF in aggregated tumour samples > 0.2 and mutation present in all tumour samples 
muts_likely_clonal_tumour = twins_filtered_dt[vaf_all_tumour > 0.2 & sum_tumour_vaf==10, mut_ID] %>% unlist()
paste('Number of likely clonal mutations:', length(muts_likely_clonal_tumour)) # 35 

######################################################################################################
# ESTIMATING THE LEVEL OF TUMOUR INFILTRATION

# Plots of tumour VAF against normal VAF
# This identifies pure samples (clear separation of two clusters) vs contaminated samples (no separation)
for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, vaf_all_tumour]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt, aes(x=vaf_tumour, y=vaf_normal))+
    geom_point()+
    theme_classic(base_size=15)+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    labs(x = 'VAF (total tumour)', y = glue('VAF in {sample}'))+
    ggtitle(glue('{sample}'))+
    coord_equal(ratio = 1)
  ggsave(glue('Figures/F2/20241208_vaf_tumour_vs_normal_{sample}_35muts.pdf'), width=4.5, height=4.5)
}

# Add labels for normal PD63383 samples which are clean (show no tumour infiltration) 
samples_normal_PD63383clean = c("PD63383w",  "PD63383t",  "PD63383u",  "PD63383ae", "PD63383ak") 
samples_normal_PD63383clean_mtr = paste0(samples_normal_PD63383clean, '_MTR')
samples_normal_PD63383clean_dep = paste0(samples_normal_PD63383clean, '_DEP')
samples_normal_PD63383clean_vap = paste0(samples_normal_PD63383clean, '_VAF')
twins_filtered_dt[, total_normal_PD63383clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383clean_mtr]
twins_filtered_dt[, total_normal_PD63383clean_dep := rowSums(.SD), .SDcols = samples_normal_PD63383clean_dep]
twins_filtered_dt[, total_normal_PD63383clean_vaf := total_normal_PD63383clean_mtr / total_normal_PD63383clean_dep]

# Approach 1: use the clean split in PD63383 samples

# first, retain mutations w/ clean split in PD63383 (indicates purity)
# remove mutations with higher aggregate in PD62341 (shared early embryonic)

# select mutations < 0.1 in PD62341
twins_filtered_dt[, col_mut_PD62341 := fcase(
  vaf_all_normal_PD62341 < 0.1, 'tumour (not in PD62341)',
  vaf_all_normal_PD62341 >= 0.1, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD62341]) # 28 tumour
table(twins_filtered_dt[, col_mut_PD62341]) # 246 tumour

# select mutations == 0 in PD63383 clean normal samples
twins_filtered_dt[, col_mut_PD63383 := fcase(
  total_normal_PD63383clean_vaf == 0, 'tumour (not in PD63383)',
  total_normal_PD63383clean_vaf > 0, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD63383]) # 27 tumour
table(twins_filtered_dt[, col_mut_PD63383]) # 216 tumour

twins_filtered_dt[, col_mut_both := fcase(
  (col_mut_PD63383 == 'tumour (not in PD63383)' & col_mut_PD62341 == 'tumour (not in PD62341)'), 'tumour-specific',
  (col_mut_PD63383 != 'tumour (not in PD63383)' | col_mut_PD62341 != 'tumour (not in PD62341)'), 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_both]) # 26 tumour
table(twins_filtered_dt[, col_mut_both]) # 214 tumour

# Estimate infiltration using the set of tumour-specific mutations
muts_tumour_specific = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour & col_mut_both=='tumour-specific', mut_ID] %>% unlist()
paste('Number of mutations classified as tumour-specific:', length(muts_tumour_specific)) # 26

# Plot distribution in each sample
for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, c('vaf_all_tumour', 'col_mut_both')]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt %>% arrange(col_mut_both), aes(x=vaf_all_tumour, y=vaf_normal, col=col_mut_both, order=col_mut_both))+
    geom_point(alpha = 0.8)+
    theme_classic(base_size=17)+
    scale_color_manual(values = c(col_normal, col_tumour))+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    labs(x = 'VAF (total tumour)', y = glue('VAF ({sample})'), col = 'Mutation category')+
    ggtitle(glue('{sample}'))+
    coord_equal(ratio=1)
  ggsave(glue('Figures/F2/20241208_vaf_tumour_vs_normal_{sample}_35muts_tumourVsShared.pdf'), width=6, height=4)
}

# Plot histogram of VAF for selected mutations in each tumour sample 
for (sample in samples_tumour){
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_dt[mut_ID %in% muts_tumour_specific, ..s_vaf] %>% unlist()
  
  pdf(glue('Figures/F2/20241208_hist_tumour_vaf_all_{sample}_26muts.pdf'), width=4.5, height=3.5)
  hist(sample_vaf_all, xlab = 'VAF', main = glue('{sample}'), xlim = c(0, 1))
  abline(v=median(sample_vaf_all),col="blue")
  abline(v=mean(sample_vaf_all),col='red')
  dev.off()
}

# Based on median VAF of 26 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
twins_filtered_vaf = twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE]
median_VAFs = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific, 2:23], median) 
median_VAFs_dt = data.frame(median_VAFs)
median_VAFs_dt = data.table(median_VAFs_dt %>% rownames_to_column('sample'))
median_VAFs_dt[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
median_VAFs_dt[, tumour_cell_fraction := 2 * median_VAFs]
median_VAFs_dt[, normal_cell_fraction := 1 - 2 * median_VAFs]
median_VAFs_dt[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
median_VAFs_dt[, purity_est := round(purity, 1)]

# Based on mean VAF of 26 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
mean_VAFs = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific, 2:23], mean) 
mean_VAFs_dt = data.frame(mean_VAFs)
mean_VAFs_dt = data.table(mean_VAFs_dt %>% rownames_to_column('sample'))
mean_VAFs_dt[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
mean_VAFs_dt[, tumour_cell_fraction := 2 * mean_VAFs]
mean_VAFs_dt[, normal_cell_fraction := 1 - 2 * mean_VAFs]
mean_VAFs_dt[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
mean_VAFs_dt[, purity_est := round(purity, 1)]

# Barchart to show tumour and normal cell fraction for each sample
median_VAFs_melt = data.table::melt(median_VAFs_dt[, c('sample', 'tumour_cell_fraction', 'normal_cell_fraction'), with=FALSE], id.vars = 'sample')
median_VAFs_melt[, sample2 := tstrsplit(sample, '_VAF', fixed = TRUE, keep = 1)]
median_VAFs_melt[, status := factor(fcase(
  sample2 %in% samples_normal, 'normal',
  sample2 %in% samples_tumour, 'tumour'))]
median_VAFs_melt[, twin := factor(fcase(
  sample2 %in% samples_PD62341, 'PD62341',
  sample2 %in% samples_PD63383, 'PD63383'))]
median_VAFs_melt[, sample_type := paste(status, twin, sep = ', ')]
median_VAFs_melt[, composition := tstrsplit(variable, '_', fixed = 1, keep = T)]

ggplot(median_VAFs_melt, aes(x=sample2, y=value, fill=composition))+
  geom_bar(stat = 'identity')+
  theme_classic(base_size=14)+
  scale_fill_manual(values = c(col_normal, col_tumour))+
  facet_wrap(~status, scales = 'free', nrow=1)+
  labs(x = 'Sample', y = glue('Cell fraction'), col = 'Cell type')+
  ggtitle(glue('Estimates of tumour and normal cell fraction for each sample'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Figures/F2/20241208_tumour_normal_content_est_26muts.pdf'), width=12, height=4)

######################################################################################################
# OUTPUT 1: SAVE TABLES WITH ESTIMATES OF NORMAL + TUMOUR CELL FRACTION AND PURITY
write.table(median_VAFs_dt, 'Out/F2/F2_stimates_tumour_cont_26muts_median_20241208.csv', sep = ',', quote=F, row.names=F)
write.table(mean_VAFs_dt, 'Out/F2/F2_estimates_tumour_cont_26muts_mean_20241208.csv', sep = ',', quote=F, row.names=F)

# these are the estimates that I will be using for the analysis going forward 

######################################################################################################
# ALL DONE

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# ALTERNATIVE WAYS OF ESTIMATING TUMOUR INFILTRATION

# Approach 2: select mutations at VAF < 0.05 in both PD63383 and PD62341

# select mutations < 0.05 in PD62341
twins_filtered_dt[, col_mut_PD62341_2 := fcase(
  vaf_all_normal_PD62341 < 0.05, 'tumour (not in PD62341)',
  vaf_all_normal_PD62341 >= 0.05, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD62341_2]) # 8 tumour
table(twins_filtered_dt[, col_mut_PD62341_2]) # 189 tumour

# select mutations < 0.05 in PD63383
twins_filtered_dt[, col_mut_PD63383_2 := fcase(
  vaf_all_normal_PD63383 < 0.05, 'tumour (not in PD63383)',
  vaf_all_normal_PD63383 >= 0.05, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD63383_2]) # 34 tumour
table(twins_filtered_dt[, col_mut_PD63383_2]) # 238 tumour

# indicate if mutation fulfills the condition (< 0.05 in PD62341 and < 0.05 in PD63383)
twins_filtered_dt[, col_mut_both_2 := fcase(
  (col_mut_PD63383_2 == 'tumour (not in PD63383)' & col_mut_PD62341_2 == 'tumour (not in PD62341)'), 'tumour-specific',
  (col_mut_PD63383_2 != 'tumour (not in PD63383)' | col_mut_PD62341_2 != 'tumour (not in PD62341)'), 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_both_2]) # 8 tumour
table(twins_filtered_dt[,col_mut_both_2]) # 176 tumour

# Estimate infiltration using the set of tumour-specific mutations
muts_tumour_specific_2 = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour & col_mut_both_2=='tumour-specific', mut_ID] %>% unlist()
paste('Number of mutations classified as tumour-specific:', length(muts_tumour_specific_2)) # 8

# Plot the distribution of tumour-specific mutations (VAF in normal sample vs VAF in tumour)
for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, c('vaf_all_tumour', 'col_mut_both_2')]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt %>% arrange(col_mut_both_2), aes(x=vaf_all_tumour, y=vaf_normal, col=col_mut_both_2, order=col_mut_both_2))+
    geom_point(alpha = 0.8)+
    theme_classic(base_size=15)+
    scale_color_manual(values = c(col_normal, col_tumour))+
    guides(alpha = "none")+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    coord_equal(ratio=1)+
    labs(x = 'VAF (total tumour)', y = glue('VAF ({sample})'), col = 'Mutation category')+
    ggtitle(glue('{sample}'))
  ggsave(glue('Figures/F2/20241208_vaf_tumour_vs_normal_{sample}_col_both_8muts.pdf'), width=6, height=4.5)
}

# Plot histogram of VAF for selected mutations in each tumour sample 
for (sample in samples_tumour){
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_vaf[mut_ID %in% muts_tumour_specific_2, ..s_vaf] %>% unlist()
  pdf(glue('Figures/F2/20241208_hist_tumour_vaf_all_{sample}_8muts.pdf'), width=4.5, height=3.5)
  hist(sample_vaf_all, xlab = 'VAF', main = glue('{sample}'), xlim = c(0, 1))
  abline(v=median(sample_vaf_all),col="blue")
  abline(v=mean(sample_vaf_all),col='red')
  dev.off()
}

# Based on median VAF of 8 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
median_VAFs_2 = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific_2, 2:23], median) 
median_VAFs_dt_2 = data.frame(median_VAFs_2)
median_VAFs_dt_2 = data.table(median_VAFs_dt_2 %>% rownames_to_column('sample'))
median_VAFs_dt_2[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
median_VAFs_dt_2[, tumour_cell_fraction := 2 * median_VAFs_2]
median_VAFs_dt_2[, normal_cell_fraction := 1 - 2 * median_VAFs_2]
median_VAFs_dt_2[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
median_VAFs_dt_2[, purity_est := round(purity, 1)]

# Based on mean VAF of 8 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
mean_VAFs_2 = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific_2, 2:23], mean) 
mean_VAFs_dt_2 = data.frame(mean_VAFs_2)
mean_VAFs_dt_2 = data.table(mean_VAFs_dt_2 %>% rownames_to_column('sample'))
mean_VAFs_dt_2[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
mean_VAFs_dt_2[, tumour_cell_fraction := 2 * mean_VAFs_2]
mean_VAFs_dt_2[, normal_cell_fraction := 1 - 2 * mean_VAFs_2]
mean_VAFs_dt_2[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
mean_VAFs_dt_2[, purity_est := round(purity, 1)]

# OUTPUT 2: SAVE TABLES WITH ESTIMATES OF NORMAL + TUMOUR CELL FRACTION AND PURITY
write.table(median_VAFs_dt_2, 'Out/F2/F2_estimates_tumour_cont_8muts_median_20241208.csv', sep = ',', quote=F, row.names=F)
write.table(mean_VAFs_dt_2, 'Out/F2/F2_estimates_tumour_cont_8muts_mean_20241208.csv', sep = ',', quote=F, row.names=F)

######################################################################################################
# Approach 3: identify deviation from VAF = 1 in tumour samples for reads on retained chr1 / chr18
# Note: makes the assumption that chr1 / chr18 losses are early (truncal in the tumour, not subclonal)

# identify mutations present at very high VAF in chr1 / chr18 in tumour samples (PD63383ap, PD63383aq)
PD63383_tumour_chr1_18_apq = twins_filtered_dt[PD63383ap_VAF > 0.7 | PD63383aq_VAF > 0.7]
# chr18_71485446_G_A # looks okay
# chr1_60839899_G_A # looks okay
# chr1_69984580_G_A # looks okay

muts_chr1_18 = c('chr18_71485446_G_A', 'chr1_60839899_G_A', 'chr1_69984580_G_A')

# Based on median VAF of 7 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
median_VAFs_3 = sapply(twins_filtered_vaf[mut_ID %in% muts_chr1_18, 2:23], median) 
median_VAFs_dt_3 = data.frame(median_VAFs_3)
median_VAFs_dt_3 = data.table(median_VAFs_dt_3 %>% rownames_to_column('sample'))
median_VAFs_dt_3[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
median_VAFs_dt_3[, tumour_cell_fraction := median_VAFs_3]
median_VAFs_dt_3[, normal_cell_fraction := 1 - median_VAFs_3]
median_VAFs_dt_3[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
median_VAFs_dt_3[, purity_est := round(purity, 1)]

# Based on mean VAF of 7 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
mean_VAFs_3 = sapply(twins_filtered_vaf[mut_ID %in% muts_chr1_18, 2:23], mean) 
mean_VAFs_dt_3 = data.frame(mean_VAFs_3)
mean_VAFs_dt_3 = data.table(mean_VAFs_dt_3 %>% rownames_to_column('sample'))
mean_VAFs_dt_3[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
mean_VAFs_dt_3[, tumour_cell_fraction := mean_VAFs_3]
mean_VAFs_dt_3[, normal_cell_fraction := 1 - mean_VAFs_3]
mean_VAFs_dt_3[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
mean_VAFs_dt_3[, purity_est := round(purity, 1)]

# OUTPUT 3: SAVE TABLES WITH ESTIMATES OF NORMAL + TUMOUR CELL FRACTION AND PURITY
write.table(median_VAFs_dt_3, 'Out/F2/F2_estimates_tumour_cont_3muts_median_chr1_18_20241208.csv', sep = ',', quote=F, row.names=F)
write.table(mean_VAFs_dt_3, 'Out/F2/F2_estimates_tumour_cont_3muts_mean_chr1_18_20241208.csv', sep = ',', quote=F, row.names=F)




