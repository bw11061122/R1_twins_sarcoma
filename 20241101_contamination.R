###################################################################################################################################
# SCRIPT 3

# Script to estimate purity of each sample (normal + tumour cell content in each sample)
# 2024-11-01
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters (1,134, 14/11/2024)

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
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241114_568.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 658

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
col_bar = '#e87811'

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
# CONTAMINATION 

# For true tumour mutations in normal samples, we expect two clusters of mutations present in the tumour
# 1 mutations shared between normal samples and tumour samples (embryonic mosaic)
# 2 mutations private to the tumour and absent from normal samples (expected VAF in normal = 0)

# 1 select mutations which are expected to be clonal in the tumour 
pdf('Results/20241114_p3_hist_658muts_mean_median.pdf', width = 4, height = 4)
hist(twins_filtered_dt[, vaf_all_tumour], breaks=30, xlab = 'VAF', xlim=c(0,1), main = 'VAF in aggregated tumour samples\n(658 mutations)')
abline(v = median(twins_filtered_dt[, vaf_all_tumour]), col = 'blue', lwd = 2)
abline(v = mean(twins_filtered_dt[, vaf_all_tumour]), col = 'red', lwd = 2)
dev.off()

# indicate the threshold above which mutations can be considered clonal
pdf('Results/20241114_p3_hist_658muts.pdf', width = 4, height = 4)
hist(twins_filtered_dt[, vaf_all_tumour], breaks=30, xlab = 'VAF', xlim=c(0,1), main = 'VAF in aggregated tumour samples\n(658 mutations)')
abline(v = 0.2, col = 'purple', lwd = 3)
dev.off()

# require VAF in aggregated tumour samples > 0.2 and mutation present in all tumour samples 
muts_likely_clonal_tumour = twins_filtered_dt[vaf_all_tumour > 0.2 & sum_tumour_vaf==10, mut_ID] %>% unlist()
paste('Number of likely clonal mutations:', length(muts_likely_clonal_tumour)) # 227 

######################################################################################################
# CONTAMINATION
######################################################################################################

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
  ggsave(glue('Results/20241114_p3_vaf_tumour_vs_normal_{sample}_227.pdf'), width=4.5, height=4.5)
}

# Add labels for normal PD63383 samples which are clean (show no contamination) 
samples_normal_PD63383clean = c("PD63383w",  "PD63383t",  "PD63383u",  "PD63383ae", "PD63383ak") 
samples_normal_PD63383clean_mtr = paste0(samples_normal_PD63383clean, '_MTR')
samples_normal_PD63383clean_dep = paste0(samples_normal_PD63383clean, '_DEP')
samples_normal_PD63383clean_vap = paste0(samples_normal_PD63383clean, '_VAF')
twins_filtered_dt[, total_normal_PD63383clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383clean_mtr]
twins_filtered_dt[, total_normal_PD63383clean_dep := rowSums(.SD), .SDcols = samples_normal_PD63383clean_dep]
twins_filtered_dt[, total_normal_PD63383clean_vaf := total_normal_PD63383clean_mtr / total_normal_PD63383clean_dep]

######################################################################################################
# Approach 1: use the clean split in PD63383 samples

# first, retain mutations w/ clean split in PD63383 (indicates purity)
# remove mutations with higher aggregate in PD62341 (shared early embryonic)

# select mutations < 0.1 in PD62341
twins_filtered_dt[, col_mut_PD62341 := fcase(
  vaf_all_normal_PD62341 < 0.1, 'tumour (not in PD62341)',
  vaf_all_normal_PD62341 >= 0.1, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD62341]) # 26 tumour
table(twins_filtered_dt[, col_mut_PD62341]) # 295 tumour

# select mutations == 0 in PD63383 clean normal samples
twins_filtered_dt[, col_mut_PD63383 := fcase(
  total_normal_PD63383clean_vaf == 0, 'tumour (not in PD63383)',
  total_normal_PD63383clean_vaf > 0, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD63383]) # 25 tumour
table(twins_filtered_dt[, col_mut_PD63383]) # 221 tumour

twins_filtered_dt[, col_mut_both := fcase(
  (col_mut_PD63383 == 'tumour (not in PD63383)' & col_mut_PD62341 == 'tumour (not in PD62341)'), 'tumour-specific',
  (col_mut_PD63383 != 'tumour (not in PD63383)' | col_mut_PD62341 != 'tumour (not in PD62341)'), 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_both]) # 24 tumour
table(twins_filtered_dt[, col_mut_both]) # 219 tumour

# Estimate contamination using the set of tumour-specific mutations
muts_tumour_specific = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour & col_mut_both=='tumour-specific', mut_ID] %>% unlist()
paste('Number of mutations classified as tumour-specific:', length(muts_tumour_specific)) # 24

# Plot distribution in each sample
for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, c('vaf_all_tumour', 'col_mut_both')]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt %>% arrange(col_mut_both), aes(x=vaf_all_tumour, y=vaf_normal, col=col_mut_both, alpha=col_mut_both, order=col_mut_both))+
    geom_point()+
    theme_classic(base_size=15)+
    scale_color_manual(values = c(col_normal, col_tumour))+
    scale_alpha_discrete(range = c(0.6, 0.9))+
    guides(alpha = "none")+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    labs(x = 'VAF (total tumour)', y = glue('VAF ({sample})'), col = 'Mutation category')+
    ggtitle(glue('{sample}'))+
    coord_equal(ratio=1)
  ggsave(glue('Results/20241114_p3_vaf_tumour_vs_normal_{sample}_col_both_24.pdf'), width=6, height=4.5)
}

# Plot histogram of VAF for selected mutations in each tumour sample 
for (sample in samples_tumour){
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_vaf[mut_ID %in% muts_tumour_specific, ..s_vaf] %>% unlist()
  
  pdf(glue('Results/20241114_p3_hist_tumour_vaf_all_{sample}_24muts.pdf'), width=4.5, height=3.5)
  hist(sample_vaf_all, xlab = 'VAF', main = glue('{sample}'), xlim = c(0, 1))
  abline(v=median(sample_vaf_all),col="blue")
  abline(v=mean(sample_vaf_all),col='red')
  dev.off()
}

# Based on median VAF of 27 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
median_VAFs = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific, 2:23], median) 
median_VAFs_dt = data.frame(median_VAFs)
median_VAFs_dt = data.table(median_VAFs_dt %>% rownames_to_column('sample'))
median_VAFs_dt[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
median_VAFs_dt[, tumour_cell_fraction := 2 * median_VAFs]
median_VAFs_dt[, normal_cell_fraction := 1 - 2 * median_VAFs]
median_VAFs_dt[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
median_VAFs_dt[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
median_VAFs_dt[, purity_est := round(purity, 1)]

# Based on mean VAF of 27 tumour-specific mutations, determine fraction of normal and tumour cells in each sample 
mean_VAFs = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific, 2:23], mean) 
mean_VAFs_dt = data.frame(mean_VAFs)
mean_VAFs_dt = data.table(mean_VAFs_dt %>% rownames_to_column('sample'))
mean_VAFs_dt[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
mean_VAFs_dt[, tumour_cell_fraction := 2 * mean_VAFs]
mean_VAFs_dt[, normal_cell_fraction := 1 - 2 * mean_VAFs]
mean_VAFs_dt[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
mean_VAFs_dt[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
mean_VAFs_dt[, purity_est := round(purity, 1)]

######################################################################################################
# OUTPUT 1: SAVE TABLES WITH ESTIMATES OF NORMAL + TUMOUR CELL FRACTION AND PURITY
write.table(median_VAFs_dt, 'Data/202411014_estimates_tumour_cont_24muts_median.csv', sep = ',', quote=F, row.names=F)
write.table(mean_VAFs_dt, 'Data/20241114_estimates_tumour_cont_24muts_mean.csv', sep = ',', quote=F, row.names=F)

######################################################################################################
# ALTERNATIVE WAYS OF ESTIMATING CONTAMINATION 
######################################################################################################

######################################################################################################
# Approach 2: select mutations at VAF < 0.05 in both PD63383 and PD62341

# select mutations < 0.05 in PD62341
twins_filtered_dt[, col_mut_PD62341_2 := fcase(
  vaf_all_normal_PD62341 < 0.05, 'tumour (not in PD62341)',
  vaf_all_normal_PD62341 >= 0.05, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD62341_2]) # 8 tumour
table(twins_filtered_dt[, col_mut_PD62341_2]) # 219 tumour

# select mutations < 0.05 in PD63383
twins_filtered_dt[, col_mut_PD63383_2 := fcase(
  vaf_all_normal_PD63383 < 0.05, 'tumour (not in PD63383)',
  vaf_all_normal_PD63383 >= 0.05, 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_PD63383_2]) # 35 tumour
table(twins_filtered_dt[, col_mut_PD63383_2]) # 270 tumour

# indicate if mutation fulfills the condition (< 0.05 in PD62341 and < 0.05 in PD63383)
twins_filtered_dt[, col_mut_both_2 := fcase(
  (col_mut_PD63383_2 == 'tumour (not in PD63383)' & col_mut_PD62341_2 == 'tumour (not in PD62341)'), 'tumour-specific',
  (col_mut_PD63383_2 != 'tumour (not in PD63383)' | col_mut_PD62341_2 != 'tumour (not in PD62341)'), 'shared embryonic')]
table(twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, col_mut_both_2]) # 8 tumour
table(twins_filtered_dt[,col_mut_both_2]) # 197 tumour

# Estimate contamination using the set of tumour-specific mutations
muts_tumour_specific_2 = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour & col_mut_both_2=='tumour-specific', mut_ID] %>% unlist()
paste('Number of mutations classified as tumour-specific:', length(muts_tumour_specific_2)) # 8

# Plot the distribution of tumour-specific mutations (VAF in normal sample vs VAF in tumour)
for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, c('vaf_all_tumour', 'col_mut_both_2')]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt %>% arrange(col_mut_both_2), aes(x=vaf_all_tumour, y=vaf_normal, col=col_mut_both_2, alpha=col_mut_both_2, order=col_mut_both_2))+
    geom_point()+
    theme_classic(base_size=15)+
    scale_color_manual(values = c(col_normal, col_tumour))+
    scale_alpha_discrete(range = c(0.6, 0.9))+
    guides(alpha = "none")+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    coord_equal(ratio=1)+
    labs(x = 'VAF (total tumour)', y = glue('VAF ({sample})'), col = 'Mutation category')+
    ggtitle(glue('{sample}'))
  ggsave(glue('Results/20241114_p3_vaf_tumour_vs_normal_{sample}_col_both_8.pdf'), width=6, height=4.5)
}

# Plot histogram of VAF for selected mutations in each tumour sample 
for (sample in samples_tumour){
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_vaf[mut_ID %in% muts_tumour_specific_2, ..s_vaf] %>% unlist()
  pdf(glue('Results/20241114_p3_hist_tumour_vaf_all_{sample}_8muts.pdf'), width=4.5, height=3.5)
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
median_VAFs_dt_2[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
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
mean_VAFs_dt_2[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
mean_VAFs_dt_2[, purity := fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction)]
mean_VAFs_dt_2[, purity_est := round(purity, 1)]

######################################################################################################
# OUTPUT 2: SAVE TABLES WITH ESTIMATES OF NORMAL + TUMOUR CELL FRACTION AND PURITY
write.table(median_VAFs_dt_2, 'Data/20241114_estimates_tumour_cont_8muts_median.csv', sep = ',', quote=F, row.names=F)
write.table(mean_VAFs_dt_2, 'Data/20241114_estimates_tumour_cont_8muts_mean.csv', sep = ',', quote=F, row.names=F)

######################################################################################################
# Approach 3: identify deviation from VAF = 1 in tumour samples for reads on retained chr1 / chr18
# Note: makes the assumption that chr1 / chr18 losses are early (truncal in the tumour, not subclonal)

# identify mutations present at very high VAF in chr1 / chr18 in tumour samples (PD63383ap, PD63383aq)
PD63383_tumour_chr1_18_apq = twins_filtered_dt[PD63383ap_VAF > 0.7 | PD63383aq_VAF > 0.7]
# chr18_64267234_G_A # suspicious mapping (checked in Blat - maps well to many regions)
# chr18_71485446_G_A # looks okay
# chr1_60839899_G_A # looks okay
# chr1_69984580_G_A # looks okay

muts_chr1_18_retained = c('chr18_71485446_G_A', 'chr1_60839899_G_A', 'chr1_69984580_G_A')
PD63383_tumour_chr1_18_apq = PD63383_tumour_chr1_18_apq[mut_ID %in% muts_chr1_18_retained]

# the opposite way is to estimate the proportion of cells which carry a mutation on the lost chr
PD63383_tumour_chr1_18_apq_lost = twins_filtered_dt[(PD63383ap_VAF < 0.1 | PD63383aq_VAF < 0.1) & vaf_all_normal > 0.3]
# "chr13_69946674_A_G" # not 100% sure with these as the read mate doesn't map correctly 
# "chr13_69946688_A_G" # not 100% sure with these as the read mate doesn't map correctly
# "chr13_69946752_A_G" # not 100% sure with these as the read mate doesn't map correctly
# "chr18_40915401_C_T" # looks okay
# "chr18_56506391_G_C" # quite a lot of insertions, wouldn't risk it 
# "chr18_56725490_T_G" # several different issues 
# "chr18_60251278_G_A" # looks okay
# "chr1_111337001_A_G" # not extremely well mapped but I'd take it
# "chr1_111337295_T_C" # variable mapping (double check?)
# "chr1_12731185_G_A" # looks okay
# "chr1_13225701_G_A" # poor mapping   
# "chr1_13269094_A_G" # poor mapping 
# "chr1_1646526_A_C" # mapping not perfect but probably fine  
# "chr1_16884197_C_G" # looks okay
# "chr1_16892573_C_T" # looks okay
# "chr1_16893359_G_A" # looks okay 
# "chr1_16894027_G_A" # very interesting mutation next to this  
# "chr1_16899420_T_C" # looks okay 
# "chr1_16902600_C_A" # looks okay  
# "chr1_16932907_G_A" # looks okay 
# "chr1_1697008_G_T" # reads with deletion but would still think it's probably fine    
# "chr1_21423048_C_A" # poor mapping  
# "chr1_25403253_C_T" # looks okay 
# "chr1_3074825_G_T" # looks okay 
# "chr1_41824923_G_T" # looks okay
# "chr1_46897283_A_T" # mapping isn't great 
# "chr1_47072883_A_T" # mapping could be better but probably still fine  
# "chr1_48695863_C_T" # looks okay # there is a mutation on Jbrowse chr1_48695882_T_C that looks good but not in the pileup
# "chr1_53061735_G_A" # looks okay

muts_chr1_18_lost = c("chr18_40915401_C_T", "chr18_60251278_G_A", "chr1_111337001_A_G", "chr1_12731185_G_A", "chr1_1646526_A_C",
                      "chr1_16884197_C_G", "chr1_16892573_C_T", "chr1_16893359_G_A", "chr1_16894027_G_A", "chr1_16899420_T_C", 
                      "chr1_16902600_C_A", "chr1_16932907_G_A", "chr1_25403253_C_T", "chr1_3074825_G_T", "chr1_41824923_G_T",
                      "chr1_48695863_C_T", "chr1_53061735_G_A")
PD63383_tumour_chr1_18_apq_lost = PD63383_tumour_chr1_18_apq_lost[mut_ID %in% muts_chr1_18_lost]
chr18_lost = c("chr18_40915401_C_T", "chr18_60251278_G_A")

# save this to a list / table
cells_chr1_18_retained = sapply(PD63383_tumour_chr1_18_apq[,2:23], mean) %>% unlist()
cells_chr18_retained = sapply(PD63383_tumour_chr1_18_apq[mut_ID == 'chr18_71485446_G_A',2:23], mean)
cells_chr1_retained = sapply(PD63383_tumour_chr1_18_apq[mut_ID != 'chr18_71485446_G_A',2:23], mean)
cells_chr1_18_lost = sapply(PD63383_tumour_chr1_18_apq_lost[,2:23], mean)
cells_chr18_lost = sapply(PD63383_tumour_chr1_18_apq_lost[mut_ID %in% chr18_lost,2:23], mean)
cells_chr1_lost = sapply(PD63383_tumour_chr1_18_apq_lost[!mut_ID %in% chr18_lost,2:23], mean)

cells_est = cbind(cells_chr1_18_retained, cells_chr18_retained, cells_chr1_retained,
                  cells_chr1_18_lost, cells_chr18_lost, cells_chr1_lost)

cells_est_dt = data.table(cells_est)
cells_est_dt = data.table(cbind(rownames(cells_est), cells_est_dt))

######################################################################################################
# OUTPUT 3: SAVE TABLES WITH ESTIMATES OF NORMAL + TUMOUR CELL FRACTION AND PURITY
write.csv(cells_est_dt, 'Data/20241114_purity_estimates_chr1_18.csv')

######################################################################################################
# Plot fraction of normal and tumour cells in each dataset 





