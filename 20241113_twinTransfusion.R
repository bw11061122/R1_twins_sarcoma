###################################################################################################################################
# SCRIPT 5

# Script to analyse twin-twin transfusion (spleen contamination)
# 2024-11-13
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters (1,134, 14/11/2024)

# OUTPUT:
# 1 estimates of twin-twin transfusion in spleen samples (fraction of PD62341 and PD63383 cells in each sample)
# 2 analysis of telomere length across samples from normal tissues of twins 

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
library(VGAM)
library(grid)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# INPUT FILES 

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Create a dataframe that only includes mutations retained post filtering  
muts = read.table('Data/mutations_include_20241114_599.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 599
twins_filtered_dt = twins_dt[mut_ID %in% muts]

###################################################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
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
col_PD62341_spleen = '#047247'
col_PD63383_spleen = '#430970'
col_PD63383_skin = '#c7b3d7'

######################################################################################################
# SAMPLES

# create lists of possible samples of interest
samples_names = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb","PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD63383ap", "PD63383aq", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
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
# Accounting for twin-twin transfusion (spleen samples) - I want to have quantitative estimates of how much transfer there is 

# Subset the dataframe to look at twin-specific / twin-enriched mutations 
mut_PD62341_dt = twins_dt[mut_ID %in% muts_PD62341_normal, c('mut_ID', samples_vaf), with=FALSE]
mut_PD62341_melt = data.table::melt(mut_PD62341_dt, id.vars = 'mut_ID')
mut_PD62341_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD62341_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', 
  sample %in% samples_tumour, 'tumour'))]
mut_PD62341_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', 
  sample %in% samples_PD63383, 'PD63383'))]
mut_PD62341_melt[, sample_type := as.factor(paste(twin, status, sep = ' '))]
mut_PD62341_melt[, sample_type2 := as.factor(fcase( 
  sample == 'PD62341v', 'PD62341 normal, spleen',
  sample == 'PD63383w', 'PD63383 normal, spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(twin, status, sep = ' ')))]
mut_PD62341_melt[, sample_type3 := as.factor(fcase(
  sample == 'PD62341v', 'PD62341 normal, spleen',
  sample == 'PD63383w', 'PD63383 normal, spleen',
  sample == 'PD63383bb', 'PD63383 normal, skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(twin, status, sep = ' ')))]
mut_PD62341_melt[, sample_type2 := factor(sample_type2, levels = c('PD62341 normal', 'PD63383 normal',
                                                                   'PD62341 normal, spleen', 'PD63383 normal, spleen'))]

mut_PD63383_dt = twins_dt[mut_ID %in% muts_PD63383_normal, c('mut_ID', samples_vaf), with=FALSE]
mut_PD63383_melt = data.table::melt(mut_PD63383_dt, id.vars = 'mut_ID')
mut_PD63383_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD63383_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', 
  sample %in% samples_tumour, 'tumour'))]
mut_PD63383_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', 
  sample %in% samples_PD63383, 'PD63383'))]
mut_PD63383_melt[, sample_type := as.factor(paste(twin, status, sep = ' '))]
mut_PD63383_melt[, sample_type2 := as.factor(fcase(
  sample == 'PD62341v', 'PD62341 normal, spleen',
  sample == 'PD63383w', 'PD63383 normal, spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(twin, status, sep = ' ')))]
mut_PD63383_melt[, sample_type3 := as.factor(fcase( # identify skin (which is known to be contaminated with PD62341-derived tumour)
  sample == 'PD62341v', 'PD62341 normal, spleen',
  sample == 'PD63383w', 'PD63383 normal, spleen',
  sample == 'PD63383bb', 'PD63383 normal, skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(twin, status, sep = ' ')))]
mut_PD63383_melt[, sample_type2 := factor(sample_type2, levels = c('PD62341 normal', 'PD63383 normal',
                                                                   'PD62341 normal, spleen', 'PD63383 normal, spleen'))]
# Plot # jitter by twin but not by tissue 
ggplot(mut_PD62341_melt[status=='normal'], aes(x=mut_ID, y=value, colour=sample_type2, alpha=sample_type2, order=sample_type2))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_PD62341_spleen, col_PD63383_spleen))+
  theme_classic(base_size = 14)+
  scale_alpha_manual(values = c(0.4, 0.4, 0.9, 0.9))+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  guides(alpha = "none")+
  ggtitle(glue('PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241114_p5_vaf_dist_PD62341_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

ggplot(mut_PD63383_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2, alpha=sample_type2, order=sample_type2))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_PD62341_spleen, col_PD63383_spleen))+
  theme_classic(base_size = 14)+
  scale_alpha_manual(values = c(0.4, 0.4, 0.9, 0.9))+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  guides(alpha = "none")+
  ggtitle(glue('PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241114_p5_vaf_dist_PD63383_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

# Mutation on the y axis
ggplot(mut_PD62341_melt[status=='normal'], aes(x=value, y=mut_ID, colour=sample_type2, alpha=sample_type2, order=sample_type2))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_PD62341_spleen, col_PD63383_spleen))+
  theme_classic(base_size = 14)+
  scale_alpha_manual(values = c(0.4, 0.4, 0.9, 0.9))+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  guides(alpha = "none")+
  ggtitle(glue('PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.8))
ggsave(glue('Results/20241114_p5_vaf_dist_PD62341_muts_samples_normal_labelspleen_yaxis.pdf'), width=5.5, height=2.5)

ggplot(mut_PD63383_melt[status=='normal'], aes(x=value, y=mut_ID, color=sample_type2, alpha=sample_type2, order=sample_type2))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_PD62341_spleen, col_PD63383_spleen))+
  theme_classic(base_size = 14)+
  scale_alpha_manual(values = c(0.4, 0.4, 0.9, 0.9))+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  guides(alpha = "none")+
  ggtitle(glue('PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.8))
ggsave(glue('Results/20241114_p5_vaf_dist_PD63383_muts_samples_normal_labelspleen_yaxis.pdf'), width=5.5, height=3.5)


# show correlation between PD62341v for PD63383 specific mutations 
# perhaps an easier way of showing this is to compare VAFs on a sample-by-sample basis 

PD63383_vafs = twins_filtered_dt[mut_ID %in% muts_PD63383, 'vaf_all_normal_PD63383', with=FALSE] %>% unlist()
for (sample_name in samples_normal_PD62341){
  sample = paste0(sample_name, '_VAF')
  PD62341_vaf = twins_filtered_dt[mut_ID %in% muts_PD63383, ..sample] %>% unlist()
  dt = data.table(cbind(PD62341_vaf, PD63383_vafs))
  ggplot(dt, aes(x=PD63383_vafs, y=PD62341_vaf))+
    geom_point(size=2.5, alpha = 0.8)+
    theme_classic(base_size = 10)+
    labs(x = 'VAF (PD63383)', y = glue('VAF ({sample_name})'))+
    ggtitle(glue('PD63383-specific mutations\n{sample_name}'))+
    coord_equal(ratio=1)+
    xlim(c(0, 0.3))+
    ylim(c(0, 0.3))
  ggsave(glue('Results/20241114_p5_vaf_dist_PD63383_muts_samples_PD63383agg_vs_{sample_name}.pdf'), width=3, height=3)
}

PD62341_vafs = twins_filtered_dt[mut_ID %in% muts_PD62341, 'vaf_all_normal_PD62341', with=FALSE] %>% unlist()
for (sample_name in samples_normal_PD63383){
  sample = paste0(sample_name, '_VAF')
  PD63383_vaf = twins_filtered_dt[mut_ID %in% muts_PD62341, ..sample] %>% unlist()
  dt = data.table(cbind(PD63383_vaf, PD62341_vafs))
  ggplot(dt, aes(x=PD62341_vafs, y=PD63383_vaf))+
    geom_point(size=2.5, alpha = 0.8)+
    theme_classic(base_size = 10)+
    labs(x = 'VAF (PD62341)', y = glue('VAF ({sample_name})'))+
    ggtitle(glue('PD62341-specific mutations\n{sample_name}'))+
    coord_equal(ratio=1)+
    xlim(c(0, 0.6))+
    ylim(c(0, 0.3))
  ggsave(glue('Results/20241114_p5_vaf_dist_PD62341_muts_samples_PD62341agg_vs_{sample_name}.pdf'), width=3, height=2)
}

# Quantification of PD62341 spleen contamination with PD63383 spleen 
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
ggsave(glue('Results/20241114_p5_spleen_vs_nonspleen_PD62341.pdf'), width=4, height=4)

ggplot(means_PD63383muts, aes(x=PD63383_spleen, y=PD62341_spleen))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD63383 spleen)', y = 'VAF (PD62341 spleen)')+
  ggtitle(glue('PD63383-specific mutations'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.5))+
  ylim(c(0, 0.5))
ggsave(glue('Results/20241114_p5_spleen_PD62341_vs_PD63383_PD63383muts.pdf'), width=4, height=4)

# do the same for PD63383 specific mutations
means_PD63383muts_spleen_PD62341 = mut_PD63383_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD62341 = data.table(mut_PD63383_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD63383muts_spleen_PD63383 = mut_PD63383_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD63383 = data.table(mut_PD63383_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD63383muts = merge(merge(means_PD63383muts_spleen_PD62341, means_PD63383muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD63383muts_spleen_PD63383, means_PD63383muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD63383muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))
means_PD63383muts[, PD62341_spleen_fPD63383 := PD62341_spleen / PD63383_nonspleen] # contamination
means_PD63383muts[, PD63383_spleen_fPD63383 := PD63383_spleen / PD63383_nonspleen] # purity 

mut_PD62341_melt[, sample := factor(sample, levels = 
                                      c('PD62341ad', 'PD62341n', 'PD62341q', 'PD62341h', 'PD62341aa', 'PD62341v',
                                        'PD63383w', 'PD63383u', 'PD63383t', 'PD63383ae', 'PD63383ak', 'PD63383bb',
                                        'PD62341b', 'PD62341u', 'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341ap', 'PD62341am',
                                        'PD63383ap', 'PD63383aq'))]
mut_PD63383_melt[, sample := factor(sample, levels = 
                                      c('PD62341ad', 'PD62341n', 'PD62341q', 'PD62341h', 'PD62341aa', 'PD62341v',
                                        'PD63383w', 'PD63383u', 'PD63383t', 'PD63383ae', 'PD63383ak', 'PD63383bb',
                                        'PD62341b', 'PD62341u', 'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341ap', 'PD62341am',
                                        'PD63383ap', 'PD63383aq'))]


# plot VAF for each mutation across samples so we can see samples separately
for (mut in muts_PD62341){
  dt = mut_PD62341_melt[mut_ID == mut]
  ggplot(dt %>% arrange(sample_type), aes(x=sample, y=value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'sample', y = 'VAF')+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour, col_tumour_PD62341))+
    ggtitle(glue('{mut}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241114_p5_PD62341_spec_by_mut_{mut}.pdf'), width=6, height=3.5)
}

for (mut in muts_PD63383){
  dt = mut_PD63383_melt[mut_ID == mut]
  ggplot(dt %>% arrange(sample_type), aes(x=sample, y=value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'sample', y = 'VAF')+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour, col_tumour_PD62341))+
    ggtitle(glue('{mut}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241114_p5_PD63383_spec_by_mut_{mut}.pdf'), width=6, height=3.5)
}

######################################################################################################
# Checking telomere length 

# documentation on the telomerecat module: https://telomerecat.readthedocs.io/en/latest/understanding_output.html

# INPUT: load dataframes of telomere length (output from telomerecat command: 20241113_telomere_length.sh)
list_tel_files = paste0('Data/', list.files(path = "Data/", pattern = "*_telomere_length.csv"))
telomere_dt = data.table(do.call(rbind, lapply(list_tel_files, function(x) read.csv(x, stringsAsFactors = FALSE))))

# systematize column names
setnames(telomere_dt, c('Sample'), 'sample')
telomere_dt[, sample := tstrsplit(sample, '.', fixed = TRUE, keep = 1)]
telomere_dt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v',
  'PD63383ak', 'PD63383bb', 'PD63383t', 'PD63383ae', 'PD63383u', 'PD63383w',
  'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u',
  'PD63383ap', 'PD63383aq'))] # specify order for plotting of samples such that normal / tumour samples plotted next to each other
telomere_dt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
telomere_dt[, status := factor(fcase(sample %in% samples_normal, 'normal', sample %in% samples_tumour, 'tumour'))]
telomere_dt[, sample_type := factor(paste(twin, status, sep = ', '))]
telomere_dt[, sample_type2 := factor(fcase(status == 'tumour', 'tumour', 
                                                 status == 'normal' & twin == 'PD62341', 'PD62341 normal',
                                                 status == 'normal' & twin == 'PD63383', 'PD63383 normal'))]
telomere_dt[, tissue := factor(fcase(
  sample %in% c('PD62341ad', 'PD63383ak'), 'cerebellum',
  sample %in% c('PD62341aa', 'PD63383bb'), 'skin',
  sample %in% c('PD62341v', 'PD63383w'), 'spleen',
  sample %in% c('PD62341h', 'PD63383t'), 'liver',
  sample %in% c('PD62341q', 'PD63383u'), 'pancreas',
  sample %in% c('PD62341n', 'PD63383ae'), 'heart', 
  sample %in% samples_tumour, 'tumour'
))]

# check telomere length across samples
ggplot(telomere_dt, aes(x = sample, y = Length, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'Telomere length (bp)', col = 'Sample category')+
  ggtitle(glue('Telomere lengths across samples'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomere_lengths.pdf'), height = 3, width = 6.5)

ggplot(telomere_dt, aes(x = sample, y = Length, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'Telomere length (bp)', col = 'Sample category')+
  ggtitle(glue('Telomere lengths across samples'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 9000))
ggsave(glue('Results/20241114_p5_telomere_lengths_from0.pdf'), height = 3, width = 6.5)

ggplot(telomere_dt[status=='normal'], aes(x = tissue, y = Length, col = twin))+
  geom_point(size=3)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'Telomere length (bp)', col = 'Twin')+
  ggtitle(glue('Telomere lengths across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomere_lengths_by_tissue.pdf'), height = 3, width = 6.5)

ggplot(telomere_dt[status=='normal'], aes(x = tissue, y = Length, col = twin))+
  geom_point(size=3)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'Telomere length (bp)', col = 'Twin')+
  ggtitle(glue('Telomere lengths across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 9000))
ggsave(glue('Results/20241114_p5_telomere_lengths_by_tissue_from0.pdf'), height = 3, width = 6.5)

# Plot F1, F2 and F4 metrics to better understand the output 
ggplot(telomere_dt, aes(x = sample, y = F1, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'Number of reads', col = 'Sample category')+
  ggtitle(glue('Telomeric reads across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomeres_F1_by_sample.pdf'), height = 3, width = 6.5)

ggplot(telomere_dt, aes(x = sample, y = F2, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'Number of reads', col = 'Sample category')+
  ggtitle(glue('CCCTAA across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomeres_F2_by_sample.pdf'), height = 3, width = 6.5)

ggplot(telomere_dt, aes(x = sample, y = F4, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'Number of reads', col = 'Sample category')+
  ggtitle(glue('TTAGGG across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomeres_F4_by_sample.pdf'), height = 3, width = 6.5)

# Compare F1, F2 and F4 metrics 
ggplot(telomere_dt[status=='normal'], aes(x = tissue, y = F1, col = twin))+
  geom_point(size=3)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'Number of reads', col = 'Twin')+
  ggtitle(glue('Telomeric reads across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomeres_F1_by_tissue.pdf'), height = 3, width = 6.5)

ggplot(telomere_dt[status=='normal'], aes(x = tissue, y = F2, col = twin))+
  geom_point(size=3)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'Number of reads', col = 'Twin')+
  ggtitle(glue('CCCTAA across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomeres_F2_by_tissue.pdf'), height = 3, width = 6.5)

ggplot(telomere_dt[status=='normal'], aes(x = tissue, y = Length, col = twin))+
  geom_point(size=3)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'Number of reads', col = 'Twin')+
  ggtitle(glue('TTAGGG lengths across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p5_telomeres_F4_by_tissue.pdf'), height = 3, width = 6.5)

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
######################################################################################################
######################################################################################################
######################################################################################################
# ALL DONE 

