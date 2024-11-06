# 2024-11-01
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters (413, 30/10/2024)

# Script to estimate purity of each sample (normal + tumour cell content in each sample)

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
muts = read.table('Data/mutations_include_20241105_850.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 413
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
col_PD62341 = "#8909c1"
col_PD63383 = "#bca4f6"
col_tumour_PD62341 = "#980505"
col_tumour_PD63383 = "#eb6767"
col_normal_PD62341 = "#0785a5"
col_normal_PD63383 = "#70c6db"
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

# Aggregated VAF in tumour samples 
twins_filtered_dt[, sum_tumour_mtr := rowSums(.SD >= 4), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, sum_tumour_vaf := rowSums(.SD >= 0.1), .SDcols = samples_tumour_vaf]
twins_filtered_dt[, sum_normal_mtr := rowSums(.SD >= 4), .SDcols = samples_normal_mtr]
twins_filtered_dt[, sum_normal_vaf := rowSums(.SD >= 0.1), .SDcols = samples_normal_vaf]

twins_filtered_dt[, total_tumour_mtr := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, total_tumour_PD62341_mtr := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, total_tumour_PD63383_mtr := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, total_normal_mtr := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_filtered_dt[, total_normal_PD62341_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, total_normal_PD63383_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]

twins_filtered_dt[, total_tumour_dep := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_filtered_dt[, total_tumour_PD62341_dep := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_filtered_dt[, total_tumour_PD63383_dep := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_filtered_dt[, total_normal_dep := rowSums(.SD), .SDcols = samples_normal_dep]
twins_filtered_dt[, total_normal_PD62341_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_filtered_dt[, total_normal_PD63383_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]

twins_filtered_dt[, total_tumour_vaf := total_tumour_mtr / total_tumour_dep]
twins_filtered_dt[, total_tumour_PD62341_vaf := total_tumour_PD62341_mtr / total_tumour_PD62341_dep]
twins_filtered_dt[, total_tumour_PD63383_vaf := total_tumour_PD63383_mtr / total_tumour_PD63383_dep]
twins_filtered_dt[, total_normal_vaf := total_normal_mtr / total_normal_dep]
twins_filtered_dt[, total_normal_PD62341_vaf := total_normal_PD62341_mtr / total_normal_PD62341_dep]
twins_filtered_dt[, total_normal_PD63383_vaf := total_normal_PD63383_mtr / total_normal_PD63383_dep]

# add for contamination (exclude mutations present in PD63383bb - contamination)
samples_normal_PD63383clean = c("PD63383w",  "PD63383t",  "PD63383u",  "PD63383ae", "PD63383ak") 
samples_normal_PD63383clean_mtr = paste0(samples_normal_PD63383clean, '_MTR')
samples_normal_PD63383clean_dep = paste0(samples_normal_PD63383clean, '_DEP')
samples_normal_PD63383clean_vap = paste0(samples_normal_PD63383clean, '_VAF')
twins_filtered_dt[, total_normal_PD63383clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383clean_mtr]
twins_filtered_dt[, total_normal_PD63383clean_dep := rowSums(.SD), .SDcols = samples_normal_PD63383clean_dep]
twins_filtered_dt[, total_normal_PD63383clean_vaf := total_normal_PD63383clean_mtr / total_normal_PD63383clean_dep]

######################################################################################################
# CONTAMINATION
######################################################################################################

######################################################################################################
# We need to have a better set of mutations 

# 1 select mutations which are clonal in the tumour based on VAF NOT number of samples
pdf('Results/20241105_p1_hist_850muts_mean_median.pdf')
hist(twins_filtered_dt[, total_tumour_vaf], breaks=50, xlab = 'VAF (agg tumour)', xlim=c(0,1), main = 'VAF across the tumour (850)')
abline(v = median(twins_filtered_dt[, total_tumour_vaf]), col = 'blue')
abline(v = mean(twins_filtered_dt[, total_tumour_vaf]), col = 'red')
dev.off()

# would it make sense to take the mutations above the median and below 0.6?
pdf('Results/20241105_p1_hist_850muts.pdf')
hist(twins_filtered_dt[, total_tumour_vaf], breaks=50, xlab = 'VAF (agg tumour)', xlim=c(0,1), main = 'VAF across the tumour (850)')
abline(v = 0.2, col = 'purple')
dev.off()

muts_likely_clonal_tumour = twins_filtered_dt[total_tumour_vaf > 0.2, mut_ID] %>% unlist()
length(muts_likely_clonal_tumour) # 617

# 2 plots of tumour VAF against normal VAF
for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, total_tumour_vaf]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt, aes(x=vaf_tumour, y=vaf_normal))+
    geom_point()+
    theme_classic(base_size=15)+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    labs(x = 'VAF (total tumour)', y = glue('VAF in {sample}'))+
    ggtitle(glue('Min 0.2 VAF in tumour (617), {sample}'))
  ggsave(glue('Results/20241105_vaf_tumour_vs_normal_{sample}_617.pdf'), width=6, height=4.5)
}

######################################################################################################
# Way 1 of choosing mutations: < 0.05 in both PD63383 and PD62341

# select mutations < 0.05 in PD63383
twins_filtered_dt[, col_mut_PD63383 := fcase(
  total_normal_PD63383_vaf <= 0.05, 'tumour (not in PD63383)',
  total_normal_PD63383_vaf > 0.05, 'shared embryonic')]
table(twins_filtered_dt[, col_mut_PD63383]) # 148 tumour

# select mutations < 0.05 in PD62341
twins_filtered_dt[, col_mut_PD62341 := fcase(
  total_normal_PD62341_vaf <= 0.05, 'tumour (not in PD62341)',
  total_normal_PD62341_vaf > 0.05, 'shared embryonic')]
table(twins_filtered_dt[, col_mut_PD62341]) # 97 tumour

twins_filtered_dt[, col_mut_both := fcase(
  (col_mut_PD63383 == 'tumour (not in PD63383)' & col_mut_PD62341 == 'tumour (not in PD62341)'), 'tumour-specific',
  (col_mut_PD63383 != 'tumour (not in PD63383)' | col_mut_PD62341 != 'tumour (not in PD62341)'), 'shared embryonic')]
table(twins_filtered_dt[, col_mut_both]) # 79 tumour

for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, c('total_tumour_vaf', 'col_mut_both')]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt, aes(x=total_tumour_vaf, y=vaf_normal, col=col_mut_both))+
    geom_point()+
    theme_classic(base_size=15)+
    scale_color_manual(values = c(col_normal, col_tumour))+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    labs(x = 'VAF (total tumour)', y = glue('VAF in {sample}'), col = 'Mutation category')+
    ggtitle(glue('Min 0.2 VAF in tumour (617), {sample}'))
  ggsave(glue('Results/20241105_vaf_tumour_vs_normal_{sample}_617_col_both.pdf'), width=6, height=4.5)
}

# Estimate contamination using the set of tumour-specific mutations
muts_tumour_specific_18 = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour & col_mut_both=='tumour-specific', mut_ID] %>% unlist()

# Plot 
for (sample in samples_tumour){
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_vaf[mut_ID %in% muts_tumour_specific_18, ..s_vaf] %>% unlist()
  pdf(glue('Results/20241105_p2_hist_tumour_vaf_all_{sample}_18muts.pdf'), width=4.5, height=3.5)
  hist(sample_vaf_all, breaks=20, xlab = 'Variant allele frequency', 
       main = glue('VAF in sample {sample}, 18 mutations'),
       xlim = c(0, 1))
  abline(v=median(sample_vaf_all),col="blue")
  abline(v=mean(sample_vaf_all),col='red')
  dev.off()
}

median_VAFs = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific_18, 2:23], median) 
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
median_VAFs_dt[, purity := factor(fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction))]

write.table(median_VAFs_dt, 'Data/20241105_estimates_tumour_cont_18muts.csv', sep = ',', quote=F, row.names=F)

######################################################################################################
# Way 2 of mutations (use the clean split in PD63383 samples)

# the way one can justify this is by first retaining mutations w/ clean split in PD63383 
# select mutations == 0 in PD63383 clean normal samples
twins_filtered_dt[, col_mut_PD63383_2 := fcase(
  total_normal_PD63383clean_vaf == 0, 'tumour (not in PD63383)',
  total_normal_PD63383clean_vaf > 0, 'shared embryonic')]
table(twins_filtered_dt[, col_mut_PD63383_2]) # 132 tumour

# select mutations < 0.1 in PD62341
twins_filtered_dt[, col_mut_PD62341_2 := fcase(
  total_normal_PD62341_vaf < 0.1, 'tumour (not in PD62341)',
  total_normal_PD62341_vaf >= 0.1, 'shared embryonic')]
table(twins_filtered_dt[, col_mut_PD62341_2]) # 186 tumour

twins_filtered_dt[, col_mut_both_2 := fcase(
  (col_mut_PD63383_2 == 'tumour (not in PD63383)' & col_mut_PD62341_2 == 'tumour (not in PD62341)'), 'tumour-specific',
  (col_mut_PD63383_2 != 'tumour (not in PD63383)' | col_mut_PD62341_2 != 'tumour (not in PD62341)'), 'shared embryonic')]
table(twins_filtered_dt[, col_mut_both_2]) # 130 tumour

for (sample in samples_normal){
  sample_name = paste0(sample, '_VAF')
  vaf_tumour = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, c('total_tumour_vaf', 'col_mut_both_2')]
  vaf_normal = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour, ..sample_name] %>% unlist()
  dt = data.frame(cbind(vaf_tumour, vaf_normal))
  ggplot(dt, aes(x=total_tumour_vaf, y=vaf_normal, col=col_mut_both_2))+
    geom_point()+
    theme_classic(base_size=15)+
    scale_color_manual(values = c(col_normal, col_tumour))+
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    labs(x = 'VAF (total tumour)', y = glue('VAF in {sample}'), col = 'Mutation category')+
    ggtitle(glue('Min 0.2 VAF (70), {sample}'))
  ggsave(glue('Results/20241105_vaf_tumour_vs_normal_{sample}_617_col_both_70.pdf'), width=6, height=4.5)
}

# Estimate contamination using the set of tumour-specific mutations
muts_tumour_specific_70 = twins_filtered_dt[mut_ID %in% muts_likely_clonal_tumour & col_mut_both_2=='tumour-specific', mut_ID] %>% unlist()

median_VAFs_70 = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific_70, 2:23], median) 
median_VAFs_dt_70 = data.frame(median_VAFs_70)
median_VAFs_dt_70 = data.table(median_VAFs_dt_70 %>% rownames_to_column('sample'))
median_VAFs_dt_70[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
median_VAFs_dt_70[, tumour_cell_fraction := 2 * median_VAFs_70]
median_VAFs_dt_70[, normal_cell_fraction := 1 - 2 * median_VAFs_70]
median_VAFs_dt_70[, status := factor(fcase(
  sample %in% samples_normal_vaf, 'normal',
  sample %in% samples_tumour_vaf, 'tumour'))]
median_VAFs_dt_70[, purity := factor(fcase(
  status == 'normal', normal_cell_fraction,
  status == 'tumour', tumour_cell_fraction))]

write.table(median_VAFs_dt_70, 'Data/20241105_estimates_tumour_cont_70muts.csv', sep = ',', quote=F, row.names=F)

# checked on Jbrowse (up until chr20)
# chr15_48353502_C_A issues 
# chr18_71485446_G_A 1 chr copy only 
# chr1_60839899_G_A 1 chr copy only 

######################################################################################################
# Estimate contamination using the set of tumour-specific mutations

# Plot 
for (sample in samples_tumour){
  
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_vaf[mut_ID %in% muts_tumour_specific2, ..s_vaf] %>% unlist()
  
  pdf(glue('Results/20241105_p2_hist_tumour_vaf_all_{sample}_71muts.pdf'), width=4.5, height=3.5)
  hist(sample_vaf_all, breaks=20, xlab = 'Variant allele frequency', 
       main = glue('VAF in sample {sample}, 71 mutations'),
       xlim = c(0, 1))
  abline(v=median(sample_vaf_all),col="blue")
  abline(v=mean(sample_vaf_all),col='red')
  dev.off()
  
}

median_VAFs = sapply(twins_filtered_vaf[mut_ID %in% muts_tumour_specific2, 2:23], median) 
median_VAFs_dt = data.frame(median_VAFs)
write.table(median_VAFs_dt, 'Data/20241105_estimates_tumour_cont_71muts.csv', sep = ',', quote=F, row.names=T)

######################################################################################################
# Contamination option 1
# Estimate contamination of tumour samples with normal using chr1 / chr18 normal reads
# Note: it could be that the sample has several tumour clones, so everything is tumour but loss of chr1 / chr18 is subclonal
# Is there any way I can know that tumour loss is clonal and occurred early?

# identify mutations present at very high VAF in chr1 / chr18 in tumour samples (PD63383ap, PD63383aq)
PD63383_tumour_chr1_18_apq = twins_filtered_vaf[PD63383ap_VAF > 0.7 | PD63383aq_VAF > 0.7]
# chr18_64267234_G_A # suspicious mapping (checked in Blat - maps well to many regions)
# chr18_71485446_G_A # looks okay
# chr1_60839899_G_A # looks okay
# chr1_69984580_G_A # looks okay
# chr1_80734464_C_T # poor mapping 
# chr8_123301519_G_C # should have been excluded
# chr9_129402893_G_T # should have been excluded 

muts_chr1_18_retained = c('chr18_71485446_G_A', 'chr1_60839899_G_A', 'chr1_69984580_G_A')
PD63383_tumour_chr1_18_apq = PD63383_tumour_chr1_18_apq[mut_ID %in% muts_chr1_18_retained]

# the opposite way is to estimate the proportion of cells which carry a mutation on the lost chr
PD63383_tumour_chr1_18_apq_lost = twins_filtered_vaf[(PD63383ap_VAF < 0.1 | PD63383aq_VAF < 0.1) & mean_vaf_normal > 0.3]
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

write.csv(cells_est_dt, 'Data/20241104_purity_estimates_chr1_18.csv')

# alternatively, identify coordinates of the lost segments and search for mutations in these regions

######################################################################################################
# Contamination option 2: Identify mutations shared by all tumour cells 

twins_filtered_vaf[sum_tumour==10] # 125 mutations present in all tumour samples at VAF >= 10
muts_all_tumour_vaf = twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist() 
paste('Mutations present in all tumour samples (VAF):', length(muts_all_tumour_vaf)) # 125

twins_filtered_mtr[sum_tumour==10] # 79 mutations present in all tumour samples at MTR >= 4
muts_all_tumour_mtr = twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()
paste('Mutations present in all tumour samples (MTR):', length(muts_all_tumour_mtr))

muts_all_tumour_vaf_mtr = Reduce(intersect, list(muts_all_tumour_vaf, muts_all_tumour_mtr))
paste('Mutations present in all tumour samples (VAF + MTR):', length(muts_all_tumour_vaf_mtr)) # 76

# Distribution of VAF of those mutations 
hist(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf, c(samples_tumour_vaf), with=FALSE] %>% unlist(), breaks=50)
abline(v = median(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='blue')
abline(v = mean(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='red')

hist(twins_filtered_vaf[mut_ID %in% muts_all_tumour_mtr, c(samples_tumour_vaf), with=FALSE] %>% unlist(), breaks=50)
abline(v = median(twins_filtered_vaf[mut_ID %in% muts_all_tumour_mtr, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='blue')
abline(v = mean(twins_filtered_vaf[mut_ID %in% muts_all_tumour_mtr, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='red')

pdf('Results/20241101_p1_hist_vaf_76muts.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr, c(samples_tumour_vaf), with=FALSE] %>% unlist(), 
     breaks=50, xlab = 'VAF', main = 'Distribution of VAF (tumour samples)\n76 shared mutations',
     xlim = c(0, 1))
abline(v = median(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='blue')
abline(v = mean(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='red')
dev.off()

# run this for each sample separately to estimate contamination per-sample
for (sample in samples_tumour){
  
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr, ..s_vaf] %>% unlist()
  
  pdf(glue('Results/20241101_p1_hist_tumour_vaf_all_{sample}.pdf'), width=4.5, height=3.5)
  hist(sample_vaf_all, breaks=20, xlab = 'Variant allele frequency', 
       main = glue('VAF in sample {sample}, 76 mutations'),
       xlim = c(0, 1))
  abline(v=median(sample_vaf_all),col="blue")
  abline(v=mean(sample_vaf_all),col='red')
  dev.off()
  
}

median_VAFs = sapply(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr, 2:23], median) 
median_VAFs 
# clearly, these are also present in normal tissues 
# these could suggest that 1 there are some mapping errors and so not specific
# the cell that gave rise to the tumour comes from a lineage that also contributed to other tissues 

######################################################################################################
# Restrict to mutations present in tumour cells only (I guess I need to do this if I want to estimate contamination in normal cells)
# allow mutations to be present in up to 3 normal samples (because there is known contamination)

muts_all_tumour_vaf_2 = twins_filtered_vaf[sum_tumour==10 & sum_normal < 4, mut_ID] %>% unlist() 
paste('Mutations present in all tumour samples (VAF):', length(muts_all_tumour_vaf_2)) # 27 

muts_all_tumour_mtr_2 = twins_filtered_mtr[sum_tumour==10 & sum_normal < 4, mut_ID] %>% unlist() 
paste('Mutations present in all tumour samples (VAF):', length(muts_all_tumour_mtr_2)) # 25 

muts_all_tumour_vaf_mtr_2 = Reduce(intersect, list(muts_all_tumour_vaf_2, muts_all_tumour_mtr_2))
paste('Mutations present in all tumour samples (VAF + MTR):', length(muts_all_tumour_vaf_mtr_2)) # 24

# Validate on Jbrowse
muts_all_tumour_vaf_mtr_2
# "chr10_100754461_C_T" # looks okay
# "chr11_134158143_C_T" # looks okay
# "chr13_60216504_A_C"  # looks okay
# "chr14_104500332_C_T" # looks okay
# "chr15_23462705_C_A" # mapping is v funny  
# "chr15_48353502_C_A" # looks okay 
# "chr15_56691722_C_T" # looks okay 
# "chr16_32621521_C_A" # looks okay
# "chr16_4003400_G_A"  # looks okay 
# "chr17_40061856_G_C" # looks okay
# "chr17_78256883_C_T" # looks okay
# "chr1_69984580_G_A" # that seems to be present at VAF ~1 in some tumour samples 
# "chr22_17774668_G_A" # looks okay
# "chr22_30553442_G_T" # looks okay
# "chr2_133311125_G_T" # looks okay
# "chr2_199565129_G_A" # looks okay
# "chr4_147908994_T_C" # looks okay
# "chr4_178731458_A_G" # looks okay
# "chr5_136349883_C_T" # looks okay
# "chr5_37790015_G_A"  # looks okay
# "chr6_95827754_C_A"  # looks okay
# "chr9_11364991_T_A"  # looks okay
# "chrX_66719643_C_G"  # looks okay
# "chrX_68803487_C_T" # looks okay

# Distribution of VAF of those mutations 
hist(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_2, c(samples_tumour_vaf), with=FALSE] %>% unlist(), breaks=50)
abline(v = median(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_2, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='blue')
abline(v = mean(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_2, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='red')

hist(twins_filtered_vaf[mut_ID %in% muts_all_tumour_mtr_2, c(samples_tumour_vaf), with=FALSE] %>% unlist(), breaks=50)
abline(v = median(twins_filtered_vaf[mut_ID %in% muts_all_tumour_mtr_2, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='blue')
abline(v = mean(twins_filtered_vaf[mut_ID %in% muts_all_tumour_mtr_2, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='red')

pdf('Results/20241101_p1_hist_vaf_24muts.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr_2, c(samples_tumour_vaf), with=FALSE] %>% unlist(), 
     breaks=50, xlab = 'VAF', main = 'Distribution of VAF (tumour samples)\n24 shared, tumour-specific mutations',
     xlim = c(0, 1))
abline(v = median(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr_2, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='blue')
abline(v = mean(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr_2, c(samples_tumour_vaf), with=FALSE] %>% unlist()), col='red')
dev.off()

# run this for each sample separately to estimate contamination per-sample
for (sample in samples_tumour){
  
  s_vaf = paste0(sample, '_VAF')
  sample_vaf_all = twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr_2, ..s_vaf] %>% unlist()
  
  pdf(glue('Results/20241101_p1_hist_tumour_vaf_all_{sample}_24muts.pdf'), width=4.5, height=3.5)
  hist(sample_vaf_all, breaks=20, xlab = 'Variant allele frequency', 
       main = glue('VAF in sample {sample}\n24 tumour-restricted mutations'),
       xlim = c(0, 1))
  abline(v=median(sample_vaf_all),col="blue")
  abline(v=mean(sample_vaf_all),col='red')
  dev.off()
  
}

# plot estimated purity 
median_tumour_vaf = sapply(twins_filtered_vaf[mut_ID %in% muts_all_tumour_vaf_mtr_2, 2:23], median) %>% unlist()
purity_est = data.table(data.frame(median_tumour_vaf) %>% rownames_to_column('sample'))
purity_est[, median_tumour_vaf := as.numeric(as.character(median_tumour_vaf))]
purity_est[, sample := tstrsplit(sample, '_', fixed=TRUE, keep=1)]
purity_est[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
purity_est[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'))]
purity_est[, sample_type := as.factor(paste(status, twin, sep = '_'))]
purity_est[, tumour_cells := 2 * median_tumour_vaf]
purity_est[, normal_cells := 1 - 2 * median_tumour_vaf]
purity_est[, purity := fcase(
  status == 'normal', normal_cells, 
  status == 'tumour', tumour_cells)]
purity_est[, sample := factor(sample, levels = c(
  'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v', 'PD62341aa', 'PD62341ad',
  'PD62341b', 'PD62341u', 'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap',
  'PD63383t', 'PD63383u', 'PD63383w', 'PD63383ae', 'PD63383ak', 'PD63383bb',
  'PD63383ap', 'PD63383aq'
))]

ggplot(purity_est, aes(x = sample, y = purity, col = sample_type))+
  geom_point(size=4)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_classic(base_size = 16)+
  ylim(c(0, 1))+
  labs(x = 'Sample', y = 'Estimated purity', color = 'Sample')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('Purity estimates (24 tumour-specific mutations')
ggsave('Results/20241101_purity_estimates.pdf')

purity_melt = melt(purity_est[, c('sample', 'tumour_cells', 'normal_cells'), with=FALSE], id.vars = 'sample')
purity_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
purity_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'))]
purity_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
purity_melt[, sample := factor(sample, levels = c(
  'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v', 'PD62341aa', 'PD62341ad',
  'PD62341b', 'PD62341u', 'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap',
  'PD63383t', 'PD63383u', 'PD63383w', 'PD63383ae', 'PD63383ak', 'PD63383bb',
  'PD63383ap', 'PD63383aq'))]

ggplot(purity_melt, aes(x = sample, y = value, fill = variable))+
  geom_bar(stat='identity')+
  scale_fill_manual(values = c(col_normal_PD62341, col_tumour_PD63383))+
  theme_classic(base_size = 16)+
  ylim(c(0, 1))+
  facet_grid(~status, scales = 'free')+
  labs(x = 'Sample', y = 'Estimated fraction of normal / tumour cells', color = 'Sample')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('Purity estimates (24 tumour-specific mutations')
ggsave('Results/20241101_purity_estimates_stacked_bar.pdf')

######################################################################################################
######################################################################################################
# Analysis of mutation clonality and contamination in the tumour

# SNZ 2012 figure: for each sample, plot coverage vs vaf 
for (sample in samples_names){
  DEP = paste0(sample, '_DEP')
  VAF = paste0(sample, '_VAF')
  df = twins_filtered_dt[, c('mut_ID', VAF, DEP, 'loss'), with=FALSE]
  df[, chr := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
  ggplot(df, aes(x=get(names(df)[3]), y=get(names(df)[2])))+
    geom_point()+
    theme_classic()+
    labs(x = 'Total number of reads covering a base', y = 'Fraction of reads reporting a variant allele', color = 'Chromosome')+
    ggtitle(glue('Nr of variant reads vs coverage: {sample}'))
  ggsave(glue('Results/20241101_p1_vaf_vs_cov_{sample}.pdf'), width=6, height=4.5)
}
# I think because of the number of mutations that we have (400 < ~70k) this is not informative





