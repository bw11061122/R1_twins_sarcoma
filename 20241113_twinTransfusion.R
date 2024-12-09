###################################################################################################################################
# SCRIPT 4

# Script to analyse twin-twin transfusion (spleen contamination)
# November 2024 
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters (255, 08/12/2024)
# 3 telomere length file for each sample (telomere.cat)
# 4 lists of twin-specific mutations

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
twins_dt = fread('Data/pileup_merged_20241016.tsv') # import high quality pileup

# Read the dataframe with filtering flags 
twins_dt_filters = fread('Out/F1/F1_twins_dt_20241208_1069.csv', sep = ',') # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Out/F1/F1_mutations_final_20241208_255.txt') 
muts = muts$V1 %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 255

# Create dataframe with mutations that passed current filters 
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Identify twin-specific mutations
muts_assignment = fread('Out/F3/F3_muts_classes_255_20241208.csv')
muts_PD62341_normal = muts_assignment[mut_class=='PD62341_specific_embryonic', mut_ID] %>% unlist()
muts_PD63383_normal = muts_assignment[mut_class=='PD63383_specific_embryonic', mut_ID] %>% unlist()

# Files with telomere length measurements (from telomerecat command)
list_tel_files = paste0('Data/', list.files(path = "Data/", pattern = "*_telomere_length.csv"))
telomere_dt = data.table(do.call(rbind, lapply(list_tel_files, function(x) read.csv(x, stringsAsFactors = FALSE))))

###################################################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#0ac368"
col_PD63383 = "#a249e8"
col_tumour_PD62341 = "#099272"
col_tumour_PD63383 = "#6F09D4"
col_PD62341_spleen = '#047247'
col_PD63383_spleen = '#8c0392'
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
# Accounting for twin-twin transfusion (spleen samples) 
# Obtain a quantitative estimate of how much blood cell transfer there is between twins 

# Subset the dataframe to look at twin-specific / twin-enriched mutations 
# PD62341 mutations 
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
mut_PD62341_melt[, mut_ID := factor(mut_ID, levels = {
  mut_PD62341_melt[, .(mean_col = mean(value)), by = mut_ID][order(mean_col), mut_ID]})]

# PD63383 mutations 
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

mut_PD63383_melt[, mut_ID := factor(mut_ID, levels = {
  mut_PD63383_melt[, .(mean_col = mean(value)), by = mut_ID][order(mean_col), mut_ID]})]

# Plot 
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
ggsave(glue('Figures/F4/20241208_vaf_dist_PD62341_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

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
ggsave(glue('Figures/F4/20241208_vaf_dist_PD63383_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

# Mutation on the y axis
ggplot(mut_PD62341_melt[status=='normal'], aes(x=value, y=mut_ID, colour=sample_type2, alpha=sample_type2, order=sample_type2))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_PD62341_spleen, col_PD63383_spleen))+
  theme_classic(base_size = 12)+
  scale_alpha_manual(values = c(0.2, 0.2, 1, 1))+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  guides(alpha = "none")+
  ggtitle(glue('PD62341-specific mutations'))+
  xlim(c(0, 0.8))
ggsave(glue('Figures/F4/20241208_vaf_dist_PD62341_muts_samples_normal_labelspleen_yaxis.pdf'), width=6.5, height=2.5)

ggplot(mut_PD63383_melt[status=='normal'], aes(x=value, y=mut_ID, color=sample_type2, alpha=sample_type2, order=sample_type2))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_PD62341_spleen, col_PD63383_spleen))+
  theme_classic(base_size = 12)+
  scale_alpha_manual(values = c(0.2, 0.2, 1, 1))+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample category')+
  guides(alpha = "none")+
  ggtitle(glue('PD63383-specific mutations'))+
  xlim(c(0, 0.8))
ggsave(glue('Figures/F4/20241208_vaf_dist_PD63383_muts_samples_normal_labelspleen_yaxis.pdf'), width=6.5, height=3.5)

# Correlation between PD62341v for PD63383 specific mutations 
# PD62341-specific mutations
PD62341_vafs = twins_filtered_dt[mut_ID %in% muts_PD62341_normal, 'vaf_all_normal_PD62341', with=FALSE] %>% unlist()
for (sample_name in samples_normal_PD63383){
  sample = paste0(sample_name, '_VAF')
  PD63383_vaf = twins_filtered_dt[mut_ID %in% muts_PD62341_normal, ..sample] %>% unlist()
  dt = data.table(cbind(PD63383_vaf, PD62341_vafs))
  ggplot(dt, aes(x=PD62341_vafs, y=PD63383_vaf))+
    geom_point(size=2.5, alpha = 0.8, color = col_PD62341)+
    theme_classic(base_size = 12)+
    labs(x = 'VAF (PD62341)', y = glue('VAF ({sample_name})'))+
    ggtitle(glue('PD62341-specific mutations\n{sample_name}'))+
    coord_equal(ratio=1)+
    xlim(c(0, 0.6))+
    ylim(c(0, 0.6))
  ggsave(glue('Figures/F4/20241208_vaf_dist_PD62341_muts_samples_PD62341agg_vs_{sample_name}.pdf'), width=4, height=4)
}

# PD63383-specific mutations
PD63383_vafs = twins_filtered_dt[mut_ID %in% muts_PD63383_normal, 'vaf_all_normal_PD63383', with=FALSE] %>% unlist()
for (sample_name in samples_normal_PD62341){
  sample = paste0(sample_name, '_VAF')
  PD62341_vaf = twins_filtered_dt[mut_ID %in% muts_PD63383_normal, ..sample] %>% unlist()
  dt = data.table(cbind(PD62341_vaf, PD63383_vafs))
  ggplot(dt, aes(x=PD63383_vafs, y=PD62341_vaf))+
    geom_point(size=2.5, alpha = 0.8, color = col_PD63383)+
    theme_classic(base_size = 10)+
    labs(x = 'VAF (PD63383)', y = glue('VAF ({sample_name})'))+
    ggtitle(glue('PD63383-specific mutations\n{sample_name}'))+
    coord_equal(ratio=1)+
    xlim(c(0, 0.3))+
    ylim(c(0, 0.3))
  ggsave(glue('Figures/F4/20241208_vaf_dist_PD63383_muts_samples_PD63383agg_vs_{sample_name}.pdf'), width=3, height=3)
}

# Quantification of PD62341 spleen contamination with PD63383 spleen 
means_PD62341muts_spleen_PD62341 = mut_PD62341_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD62341 = data.table(mut_PD62341_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD62341muts_spleen_PD63383 = mut_PD62341_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD63383 = data.table(mut_PD62341_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD62341muts = merge(merge(means_PD62341muts_spleen_PD62341, means_PD62341muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD62341muts_spleen_PD63383, means_PD62341muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD62341muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))

# observed VAF_PD62341_spleen = VAF_PD62341_spleen * fraction_PD62341_spleen + VAF_PD63383 * fraction_PD63383 
# observed VAF_PD63383_spleen = VAF_PD62341 * fraction_PD62341 + VAF_PD63383_spleen * fraction_PD63383_spleen
means_PD62341muts[, PD62341_spleen_fPD62341 := PD62341_spleen / PD62341_nonspleen] # assume true VAF in PD63383 = 0
means_PD62341muts[, PD63383_spleen_fPD62341 := PD63383_spleen / PD62341_nonspleen] 

ggplot(means_PD62341muts, aes(x=PD62341_nonspleen, y=PD62341_spleen))+
  geom_point(size=2.5, color = col_PD62341)+ 
  theme_classic(base_size = 10)+
  labs(x = 'VAF PD62341 (non-spleen)', y = 'VAF PD62341 (spleen)')+
  ggtitle(glue('PD62341-specific mutations'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.5))+
  ylim(c(0, 0.5))
ggsave(glue('Figures/F4/20241208_spleen_vs_nonspleen_PD62341.pdf'), width=3, height=3)

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

ggplot(means_PD63383muts, aes(x=PD63383_nonspleen, y=PD62341_spleen))+
  geom_point(size=2.5, color = col_PD63383)+ 
  theme_classic(base_size = 10)+
  labs(x = 'VAF PD63383 (non-spleen)', y = 'VAF PD62341 (spleen)')+
  ggtitle(glue('PD63383-specific mutations'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.3))+
  ylim(c(0, 0.3))
ggsave(glue('Figures/F4/20241208_spleen_vs_nonspleen_PD63383.pdf'), width=3, height=3)

######################################################################################################
# Checking telomere length 

# documentation on the telomerecat module: https://telomerecat.readthedocs.io/en/latest/understanding_output.html

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
# note: if we thought this was more interesting, you should use SD to show this data properly 

# plot values for each sample separately
ggplot(telomere_dt, aes(x = sample, y = Length, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'Telomere length (bp)', col = 'Sample category')+
  ggtitle(glue('Telomere lengths across samples'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 9000))
ggsave(glue('Figures/F4/20241208_telomere_lengths_from0.pdf'), height = 3, width = 6.5)

# compare telomere length between tissues of each twin 
ggplot(telomere_dt[status=='normal'], aes(x = tissue, y = Length, col = twin))+
  geom_point(size=3)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'Telomere length (bp)', col = 'Twin')+
  ggtitle(glue('Telomere lengths across normal tissues'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 9000))
ggsave(glue('Figures/F4/20241208_telomere_lengths_by_tissue_from0.pdf'), height = 3, width = 6.5)

######################################################################################################
# Check if there are any possible driver mutations of clonal hem in the set of normal mutations which could affect uneven transfusion b/n twins?
# NB those would have to be twins-specific and present in blood cells; very unlikely 

# clonal hem driver genes (from https://pmc.ncbi.nlm.nih.gov/articles/PMC11176083/)
# ZBTB33, ZNF318, ZNF234, SPRED2, SH2B3, SRCAP, SIK3, SRSF1, CHEK2, CCDC115, CCL22, BAX, YLPM1, MYD88, MTA2, MAGEC3 and IGLL5
# DNMT3A, TET2, ASXL1, PPM1D, SRSF2, SF3B1, GNB1, IDH1, IDH2, T3, BRCC2, GNAS, JAK2, KDM6A, CBL, PHIP
CH_genes = c('ZBTB33', 'ZNF318', 'ZNF234', 'SPRED2', 'SH2B3', 'SRCAP', 'SIK3', 'SRSF1', 'CHEK2', 'CCDC115', 'CCL22', 'BAX', 'YLPM1', 'MYD88', 'MTA2', 'MAGEC3',
             'IGLL5', 'DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'SRSF2', 'SF3B1', 'GNB1', 'IDH1', 'IDH2', 'TP53', 'BRCC2', 'GNAS', 'JAK2', 'KDM6A', 'CBL', 'PHIP')

paste("Number of possible CH-genes with mutations:", length(twins_dt[Gene %in% CH_genes & Effect != 'intronic', Gene] %>% unlist() %>% unique()))
twins_dt_filters[Gene %in% CH_genes & Effect != 'intronic' & f7_likelyGermline_bothTwins==0, c('mut_ID', 'Gene', 'sum_req_filters'), with=FALSE] 
# none of those pass the QC, so likely artefacts; the rest of those mutations are likely germline

######################################################################################################
# ALL DONE 

