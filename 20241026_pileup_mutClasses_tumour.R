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
muts = read.table('Data/mutations_include_20241104_425.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 425
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
# Basic checks: number of mutations

paste('Number of all mutations in the set:', dim(twins_filtered_vaf)[1])

paste('Number of mutations in normal samples:', length(Reduce(intersect, list(twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist())))) # 381 
paste('Number of mutations in tumour samples:',  length(Reduce(intersect, list(twins_filtered_vaf[sum_tumour>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist())))) # 391

paste('Number of mutations in all normal samples:', length(Reduce(intersect, list(twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist())))) # 67
paste('Number of mutations in all tumour samples:',  length(Reduce(intersect, list(twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist())))) # 75 

paste('Number of mutations in max 3 normal samples:', length(Reduce(intersect, list(twins_filtered_vaf[sum_normal<4, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal<4, mut_ID] %>% unlist())))) # 151

# Create lists of mutations 
muts_all = twins_filtered_vaf[, mut_ID] %>% unlist()
muts_normal = Reduce(intersect, list(twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist()))
muts_tumour = Reduce(intersect, list(twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist()))
muts_normal_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()))
muts_tumour_all = Reduce(intersect, list(twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()))
muts_tumour_only = Reduce(intersect, list(twins_filtered_vaf[sum_normal<4, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal<4, mut_ID] %>% unlist()))
muts_tumour_all_only = Reduce(intersect, list(muts_tumour_all, muts_tumour_only)) %>% unlist()

######################################################################################################
# Analysis of tumour evolution: coverage vs VAF plot

twins_filtered_mtr[, sum_tumour_mtr := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_filtered_dep[, sum_tumour_dep := rowSums(.SD), .SDcols = samples_tumour_dep]

# Aggregate data for all tumour samples together 
twins_tumour_agg = merge(twins_filtered_mtr[, c('mut_ID', 'sum_tumour_mtr'), with=FALSE], 
                                            twins_filtered_dep[, c('mut_ID', 'sum_tumour_dep'), with=FALSE])
twins_tumour_agg[, tumour_vaf := sum_tumour_mtr / sum_tumour_dep]
twins_tumour_agg[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_tumour_only), 'other'))]

# Plot coverage vs VAF 
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_bw(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241103_cov_vs_vaf_mut_class.pdf'))

# mutations with max coverage:
# chr5_795141_C_T # questionable  
# chr9_42329754_C_T # seems okay 
# chr11_4594537_A_G # not really 
# chr5_784531_G_A # mapping issues 
# chr19_54861143_C_A # mapping issues, not sure 

# color by type of mutation
twins_tumour_agg[, c('Ref', 'Alt') := tstrsplit(mut_ID, '_', fixed=T, keep=3:4)]
twins_tumour_agg[, mut_type := paste0(Ref, '>', Alt)]
twins_tumour_agg[, mut_type := as.factor(fcase(
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
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = mut_type))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_bw(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241103_cov_vs_vaf_mut_type.pdf'))

# color by chromosome 
twins_tumour_agg[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = Chrom))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_bw(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241103_cov_vs_vaf_chrom.pdf'))

twins_tumour_agg[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'))]
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = loss))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_bw(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241103_cov_vs_vaf_loss.pdf'))

# add number of tumour samples the mutation is present in 
twins_tumour_agg = merge(twins_tumour_agg, twins_filtered_vaf[, c('mut_ID', 'sum_tumour'), with=FALSE], by = 'mut_ID')
twins_tumour_agg[, sum_tumour := as.factor(sum_tumour)]
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = sum_tumour))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_bw(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Number of tumour samples mutation is detected in')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241103_cov_vs_vaf_sum_tumour.pdf'))

######################################################################################################
# Plot VAF and coverage across the genome 
twins_tumour_agg[, pos := tstrsplit(mut_ID, '_', fixed=T, keep=2)]
twins_tumour_agg[, pos := as.numeric(pos)]
twins_tumour_agg[, Chrom_nr := tstrsplit(Chrom, 'chr', fixed=T, keep=2)]
twins_tumour_agg[, Chrom_nr := factor(Chrom_nr, levels = 
                                        c('1', '2', '3', '4', '5',
                                          '6', '7', '8', '9', '10',
                                          '11', '12', '13', '14', '15',
                                          '16', '17', '18', '19', '20',
                                          '21', '22', 'X'))]

ggplot(twins_tumour_agg, aes(x = pos, y = tumour_vaf, col=mut_class))+
  geom_point(size=0.8)+
  labs(x = 'Position', y = 'VAF (all tumour samples)', col='Mutation class')+
  ggtitle(glue('VAF across the genome'))+
  facet_grid(~Chrom_nr)+
  theme_minimal(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())+
#  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))
ggsave(glue('Results/20241103_vaf_across_the_genome.pdf'))

######################################################################################################
# Relationship between tumour clones (PJC 2015 Fig 1)

# let's do a trial plot for PD63383ap and PD63383aq (easy)
sub_PD63383_tumour = twins_filtered_vaf[, c('mut_ID', 'PD63383ap_VAF', 'PD63383aq_VAF'), with=FALSE]

# let's calculate the clonal cell fraction 
est_purity_PD63383ap = 0.9 # not hyper-accurate but should be about right
est_purity_PD63383aq = 0.9

# note: I am not taking into account copy number state 
# the reason is, I am not 100% sure which chr1 / chr18 segments were deleted 
sub_PD63383_tumour[, PD63383ap_ccf := 2 * (PD63383ap_VAF / est_purity_PD63383ap)]
sub_PD63383_tumour[, PD63383aq_ccf := 2 * (PD63383aq_VAF / est_purity_PD63383aq)]

ggplot(sub_PD63383_tumour, aes(x = PD63383ap_ccf, y = PD63383aq_ccf))+
  geom_point(size=.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241103_p2_ccf_PD63383ap_PD63383aq.pdf'))

sub_PD63383_tumour[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_tumour_only), 'other'))]

ggplot(sub_PD63383_tumour, aes(x = PD63383ap_ccf, y = PD63383aq_ccf, col = mut_class))+
  geom_point(size=1.5)+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', y = 'Clonal cell fraction (PD63383aq)', col = 'Mutation class')+
  ggtitle(glue('CCF plot PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241103_p2_ccf_PD63383ap_PD63383aq_mut_class.pdf'))

# create data frame with purity estimates to merge (not extremely accurate atm)
purity_dt = data.table(data.frame(
  sample = c("PD62341b", "PD62341u", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", 'PD63383ap', 'PD63383aq'),
  est_purity  = c(0.76, 0.47, 0.93, 0.31, 0.58, 0.49, 0.78, 0.87, 0.93, 0.91)))

sub_tumour = twins_filtered_vaf[, c('mut_ID', c(samples_tumour_vaf)), with=FALSE]
sub_tumour_melt = melt(sub_tumour, id.vars = 'mut_ID')
sub_tumour_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
sub_tumour_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
sub_tumour_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
sub_tumour_melt[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_tumour_only), 'other'))]
sub_tumour_melt = merge(sub_tumour_melt, purity_dt, by = 'sample')
sub_tumour_melt[, ccf := 2 * (value / est_purity)]

ccf_ap = sub_tumour_melt[sample=='PD63383ap', c('mut_ID', 'ccf', 'mut_class')] 
ccf_aq = sub_tumour_melt[sample=='PD63383aq', c('mut_ID', 'ccf')] 
dt_apq = merge(ccf_ap, ccf_aq, by = 'mut_ID')

ggplot(dt_apq, aes(x = ccf.x, y = ccf.y))+
  geom_point(size=3.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241103_p2_ccf_PD63383ap_PD63383aq.pdf'))

ggplot(dt_apq, aes(x = ccf.x, y = ccf.y, col=mut_class))+
  geom_point(size=2.5)+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', 
       y = 'Clonal cell fraction (PD63383aq)',
       col = 'Mutation class')+
  ggtitle(glue('CCF plot PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241103_p2_ccf_PD63383ap_PD63383aq_mut_class.pdf'))

ccf_am = sub_tumour_melt[sample=='PD62341am', c('mut_ID', 'ccf', 'mut_class')] 
ccf_aq = sub_tumour_melt[sample=='PD63383aq', c('mut_ID', 'ccf')] 
dt_amq = merge(ccf_am, ccf_aq, by = 'mut_ID')
setnames(dt_amq, c('ccf.x', 'ccf.y'), c('ccf_am', 'ccf_aq'))
ggplot(dt_amq, aes(x = ccf_am, y = ccf_aq))+
  geom_point(size=2.5, color = 'lightblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341am)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot: PD62341 vs PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241103_p2_ccf_PD62341am_PD63383aq.pdf'))

ggplot(dt_amq, aes(x = ccf_am, y = ccf_aq, col = mut_class))+
  geom_point(size=2.5)+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341am)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot: PD62341 vs PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241103_p2_ccf_PD62341am_PD63383aq_mut_class.pdf'))

# generate these plots for each pair of samples
for(sample1 in samples_tumour){
  for (sample2 in samples_tumour){
    if (sample1 != sample2){
      ccf1 = sub_tumour_melt[sample==sample1, c('mut_ID', 'ccf', 'mut_class'), with=FALSE] 
      ccf2 = sub_tumour_melt[sample==sample2, c('mut_ID', 'ccf'), with=FALSE] 
      dt = merge(ccf1, ccf2, by = 'mut_ID')
      ggplot(dt, aes(x = ccf.x, y = ccf.y, col = mut_class))+
        geom_point(size=2.5)+
        theme_bw(base_size = 14)+
        labs(x = glue('Clonal cell fraction in {sample1}'), y = glue('Clonal cell fraction in {sample2}'))+
        ggtitle(glue('CCF plot'))+
        xlim(c(0, 2))+
        ylim(c(0, 2))+
        geom_vline(xintercept = 1, color='black', size = 0.5)+
        geom_hline(yintercept = 1, color='black', size = 0.5)
      ggsave(glue('Results/20241103_p2_ccf_{sample1}_{sample2}.pdf'), width=8, height=6.5)
    }
  }
}

######################################################################################################
# CCF plots for different classes of mutations

# all mutations in the set 
sub_tumour_425 = twins_filtered_vaf[, c('mut_ID', c(samples_tumour_vaf)), with=FALSE]
sub_tumour_425_melt = melt(sub_tumour_425, id.vars = 'mut_ID')
sub_tumour_425_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
sub_tumour_425_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
sub_tumour_425_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'))]
sub_tumour_425_melt = merge(sub_tumour_425_melt, purity_dt, by = 'sample')
sub_tumour_425_melt[, ccf := 2 * (value / est_purity)]

ccf_425_ap = sub_tumour_425_melt[sample=='PD63383ap', ccf] 
ccf_425_aq = sub_tumour_425_melt[sample=='PD63383aq', ccf] 
dt_425_apq = data.frame(cbind(ccf_425_ap, ccf_425_aq))

ggplot(dt_425_apq, aes(x = ccf_425_ap, y = ccf_425_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 tumour (425 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD63383ap_PD63383aq_425.pdf'))

# PD62341 vs PD63383 tumour
ccf_425_am = sub_tumour_425_melt[sample=='PD62341am', ccf] 
ccf_425_aq = sub_tumour_425_melt[sample=='PD63383aq', ccf] 
dt_425_amq = data.frame(cbind(ccf_425_am, ccf_425_aq))

ggplot(dt_425_amq, aes(x = ccf_425_am, y = ccf_425_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341am)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 vs PD62341 (425 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD62341am_PD63383aq_425.pdf'))

ccf_425_ag = sub_tumour_425_melt[sample=='PD62341ag', ccf] 
ccf_425_aq = sub_tumour_425_melt[sample=='PD63383aq', ccf] 
dt_425_agq = data.frame(cbind(ccf_425_ag, ccf_425_aq))

ggplot(dt_425_agq, aes(x = ccf_425_ag, y = ccf_425_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341ag)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 vs PD62341 (425 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD62341ag_PD63383aq_425.pdf'))

# mutations present in tumour samples
purity_dt = data.table(data.frame(
  sample = c("PD62341b", "PD62341u", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", 'PD63383ap', 'PD63383aq'),
  est_purity  = c(0.76, 0.47, 0.93, 0.5, 0.58, 0.49, 0.78, 0.87, 0.93, 0.91)))

sub_tumour = twins_filtered_vaf[mut_ID %in% muts_tumour_all_only, c('mut_ID', c(samples_tumour_vaf)), with=FALSE]
sub_tumour_melt = melt(sub_tumour, id.vars = 'mut_ID')
sub_tumour_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
sub_tumour_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
sub_tumour_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]

sub_tumour_melt = merge(sub_tumour_melt, purity_dt, by = 'sample')
sub_tumour_melt[, ccf := 2 * (value / est_purity)]

ccf_ap = sub_tumour_melt[sample=='PD63383ap', ccf] 
ccf_aq = sub_tumour_melt[sample=='PD63383aq', ccf] 
dt_apq = data.frame(cbind(ccf_ap, ccf_aq))

ggplot(dt_apq, aes(x = ccf_ap, y = ccf_aq))+
  geom_point(size=3.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD63383ap_PD63383aq.pdf'))

ccf_am = sub_tumour_melt[sample=='PD62341am', ccf] 
ccf_aq = sub_tumour_melt[sample=='PD63383aq', ccf] 
dt_amq = data.frame(cbind(ccf_am, ccf_aq))

ggplot(dt_apq, aes(x = ccf_am, y = ccf_aq))+
  geom_point(size=2.5, color = 'lightblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341am)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot: PD62341 vs PD63383 tumour'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD62341am_PD63383aq.pdf'))

# mutations present in >= 1 tumour sample 
sub_tumour_393 = twins_filtered_vaf[mut_ID %in% muts_min1_tumour, c('mut_ID', c(samples_tumour_vaf)), with=FALSE]
sub_tumour_393_melt = melt(sub_tumour_393, id.vars = 'mut_ID')
sub_tumour_393_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
sub_tumour_393_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
sub_tumour_393_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'))]
sub_tumour_393_melt = merge(sub_tumour_393_melt, purity_dt, by = 'sample')
sub_tumour_393_melt[, ccf := 2 * (value / est_purity)]

ccf_393_ap = sub_tumour_393_melt[sample=='PD63383ap', ccf] 
ccf_393_aq = sub_tumour_393_melt[sample=='PD63383aq', ccf] 
dt_393_apq = data.frame(cbind(ccf_393_ap, ccf_393_aq))

ggplot(dt_393_apq, aes(x = ccf_393_ap, y = ccf_393_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 tumour (393 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD63383ap_PD63383aq_391.pdf'))

# PD62341 vs PD63383 tumour
ccf_393_am = sub_tumour_393_melt[sample=='PD62341am', ccf] 
ccf_393_aq = sub_tumour_393_melt[sample=='PD63383aq', ccf] 
dt_393_amq = data.frame(cbind(ccf_393_am, ccf_393_aq))

ggplot(dt_393_amq, aes(x = ccf_393_am, y = ccf_393_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341am)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 vs PD62341 (393 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD62341am_PD63383aq_393.pdf'))

ccf_393_ag = sub_tumour_393_melt[sample=='PD62341ag', ccf] 
ccf_393_aq = sub_tumour_393_melt[sample=='PD63383aq', ccf] 
dt_393_agq = data.frame(cbind(ccf_393_ag, ccf_393_aq))

ggplot(dt_393_agq, aes(x = ccf_393_ag, y = ccf_393_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341ag)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 vs PD62341 (393 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD62341ag_PD63383aq_393.pdf'))

# mutations absent from normal samples 
sub_tumour_151 = twins_filtered_vaf[mut_ID %in% muts_tumour_only, c('mut_ID', c(samples_tumour_vaf)), with=FALSE]
sub_tumour_151_melt = melt(sub_tumour_151, id.vars = 'mut_ID')
sub_tumour_151_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
sub_tumour_151_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
sub_tumour_151_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'))]
sub_tumour_151_melt = merge(sub_tumour_151_melt, purity_dt, by = 'sample')
sub_tumour_151_melt[, ccf := 2 * (value / est_purity)]

ccf_151_ap = sub_tumour_151_melt[sample=='PD63383ap', c('ccf', 'mut_ID'), with=F] 
ccf_151_aq = sub_tumour_151_melt[sample=='PD63383aq', c('ccf', 'mut_ID'), with=F] 
dt_151_apq = merge(ccf_151_ap, ccf_151_aq, by='mut_ID')

ggplot(dt_151_apq, aes(x = ccf.x, y = ccf.y))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD63383ap)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 tumour (153 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD63383ap_PD63383aq_151.pdf'))

# can we find something subclonal
dt_151_apq[(ccf.x < 0.5 & ccf.x > 0.1) & (ccf.y < 0.5 & ccf.y > 0.1)] # 11 
# chr10_49977949_A_G # mapping issues 
# chr11_80115755_C_T # this looks okay and CCF ~0.4 - 0.5 # PD63383-specific 
# chr13_18565005_G_A # mapping issues 
# chr19_7878006_A_G # mapping issues 
# chr1_145588334_T_C # mapping issues 
# chr1_83136026_A_T # mapping issues 
# chr2_113525233_T_C # mapping issues 
# chr2_82147172_T_C # looks okay and also low CCF # PD63383-specific 
# chr7_57104457_T_A # mapping issues 
# chr9_100061581_C_T # bit weird
# chr9_100061582_G_A # bit weird 

# PD62341 vs PD63383 tumour
ccf_151_am = sub_tumour_151_melt[sample=='PD62341am', ccf] 
ccf_151_aq = sub_tumour_151_melt[sample=='PD63383aq', ccf] 
dt_151_amq = data.frame(cbind(ccf_151_am, ccf_151_aq))

ggplot(dt_151_amq, aes(x = ccf_151_am, y = ccf_151_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341am)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 vs PD62341 (151 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD62341am_PD63383aq_151.pdf'))

ccf_151_ag = sub_tumour_151_melt[sample=='PD62341ag', ccf] 
ccf_151_aq = sub_tumour_151_melt[sample=='PD63383aq', ccf] 
dt_151_agq = data.frame(cbind(ccf_151_ag, ccf_151_aq))

ggplot(dt_151_agq, aes(x = ccf_151_ag, y = ccf_151_aq))+
  geom_point(size=2.5, color = 'darkblue')+
  theme_bw(base_size = 14)+
  labs(x = 'Clonal cell fraction (PD62341ag)', y = 'Clonal cell fraction (PD63383aq)')+
  ggtitle(glue('CCF plot PD63383 vs PD62341 (153 mutations)'))+
  xlim(c(0, 2))+
  ylim(c(0, 2))+
  geom_vline(xintercept = 1, color='black', size = 0.5)+
  geom_hline(yintercept = 1, color='black', size = 0.5)
ggsave(glue('Results/20241101_p2_ccf_PD62341ag_PD63383aq_151.pdf'))

# are these mutations present mainly in contaminated normal samples?
twins_filtered_vaf[mut_ID %in% muts_tumour_only]
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_only, sum_normal])
dim(twins_filtered_vaf[mut_ID %in% muts_tumour_only & PD63383bb_VAF > 0.1])[1] # 85
dim(twins_filtered_vaf[mut_ID %in% muts_tumour_only & PD62341aa_VAF > 0.1])[1] # 27
dim(twins_filtered_vaf[mut_ID %in% muts_tumour_only & PD62341h_VAF > 0.1])[1] # 72

# are these mutations real?
dim(twins_filtered_vaf[mut_ID %in% muts_tumour_only & PD62341h_VAF < 0.1 & PD62341aa_VAF < 0.1 & PD63383bb_VAF < 0.1 & sum_normal > 0])[1] # 16
# "chr10_119182885_G_C" # mapping issues 
# "chr11_18922568_G_T"  # variable mapping, and on blat actually maps better to sth else 
# "chr13_24392029_C_G"  # mapping issues 
# "chr19_54746517_T_C"  # kind of okay
# "chr19_7878006_A_G"   # clear mapping issues 
# "chr1_103587565_A_C"  # looks okay
# "chr1_145588334_T_C"  # some mapping issues (double check)
# "chr1_83136026_A_T"  # also mapping issues 
# "chr2_113525233_T_C" # mapping issues 
# "chr4_15905566_C_T" # looks okay  
# "chr6_159851462_T_C" # looks fine but there are some deletions so estimates may not be reliable
# "chr7_149688370_C_T"  # looks okay but v low everywhere
# "chr7_63778593_T_G"  # mapping issues  
# "chr9_100061581_C_T" # there are some insertions so not sure - double check 
# "chr9_100061582_G_A" # there are some insertions so not sure - double check
# "chr9_67625770_G_A" # mapping issues   

######################################################################################################
######################################################################################################
# Inspect tumour-specific mutations 

# I will allow up to 3 normal samples due to the possibility of contamination
# I suspect that skin samples (aa, bb) and liver sample (PD62341h) are contaminated with reads from the tumour

# Mutations found in all tumour samples 
twins_filtered_mtr[sum_tumour == 10 & sum_normal == 0] # 1 
# "chr2_199565129_G_A" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 1] # 4 
# "chr16_32621521_C_A" # mapping quality not excellent (double check)
# "chr17_78256883_C_T" # looks okay
# "chr22_17774668_G_A" # looks okay 
# "chr22_30553442_G_T" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 2] # 11
# "chr13_60216504_A_C" # looks okay
# "chr15_48353502_C_A" # looks okay
# "chr15_56691722_C_T" # looks okay
# "chr16_4003400_G_A" # looks okay
# "chr1_69984580_G_A" # can estimate contamination from this one!!
# "chr2_133311125_G_T" # looks okay
# "chr5_136349883_C_T" # looks okay
# "chr5_37790015_G_A" # looks okay
# "chr6_95827754_C_A" # mapping could be better but still decent
# "chr9_11364991_T_A" # looks okay
# "chrX_68803487_C_T" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal == 3] # 9 
# it's all bb, h and aa (skin + liver)
# chr10_100754461_C_T # looks okay
# chr11_134158143_C_T # looks okay
# chr14_104500332_C_T # looks okay
# chr15_23462705_C_A # lots of mis-mapped reads 
# chr17_40061856_G_C # looks okay
# chr4_130207996_C_T # looks okay
# chr4_147908994_T_C # looks okay
# chr4_178731458_A_G # looks okay
# chrX_66719643_C_G # looks okay

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 0] # 2 
# chr12_106354695_C_T # looks okay, missing in PD62341ak
# chrX_115495236_T_C # looks okay to me (double check but I would believe it), missing in PD62341b

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 1] # 7
# chr12_61020950_C_A # looks okay
# chr14_104122986_G_T # looks okay (nb some mis-mapping)
# chr7_115240063_T_C # looks okay
# chr7_120677593_C_T # looks okay
# chr8_46678612_T_C # looks okay
# chrX_46375010_C_T # looks okay
# chrX_56448796_C_A # looks okay

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 2] # 13
# "chr12_10337296_G_T" # variable mapping (double check)
# "chr12_39248376_A_T" # looks okay
# "chr18_71485446_G_A" # looks like del of the other copy - can estimate contamination
# "chr21_45729938_T_A" # looks okay
# "chr22_16694683_C_A" # looks okay
# "chr2_119170376_C_T" # looks okay
# "chr2_187893783_G_T" # looks okay
# "chr3_94699090_G_T" # looks okay
# "chr4_180479634_G_A" # looks okay
# "chr5_112010886_G_A" # looks okay
# "chr7_12297469_T_G" # looks okay
# "chr8_54112201_C_A" # looks okay
# "chr9_94529993_C_T" # quite a lot of mis-mapped reads (double check)

# there are two mutations next to each other that show a very similar pattern and VAF 
# the other is: chr8_54112197_T_C, and is called as present in 9 tumour samples and 3 normal
# normal chr8_54112197_T_C: PD62341aa, h, PD63383 bb
# normal chr8_54112201_C_A: PD62341h, PD63383bb (3 reads in PD62341aa)

twins_filtered_mtr[sum_tumour == 9 & sum_normal == 3] # 21
# "chr13_100405282_G_T" # looks okay
# "chr15_76646556_T_C" # looks okay  
# "chr17_5029328_G_C" # looks okay 
# "chr1_119005653_G_C" # mis-mapped reads (double check cause it really doesn't look too bad) 
# "chr1_60839899_G_A" # looks great
# "chr2_124573116_C_T" # looks great  
# "chr2_1420635_A_T" # looks great # there is a beautiful mutation next to this one 
# "chr2_2023708_C_G"  # looks great  
# "chr2_232306029_G_C" # looks great
# "chr3_49878199_C_A" # looks great  
# "chr4_46067117_A_G" # difficult region but looks okay
# "chr4_93081874_T_A" # looks great  
# "chr5_119778240_G_C" # poor mapping quality (double check)
# "chr7_139659050_G_A" # looks okay
# "chr7_25564729_G_A" # looks okay
# "chr8_131571989_A_T" # looks okay
# "chr8_54112197_T_C" # looks okay
# "chr8_97785584_T_C" # looks okay 
# "chr9_12824689_G_A" # looks okay 
# "chrX_124531709_C_T" # looks okay  
# "chrX_83934262_C_G" # looks okay

# very interesting mutation chr1_60839875_T_C 
# seems to be on the other chromosome 
# can compare reads with different muts to determine contamination in different tumour samples
# also, PD62341q for example has 5 of those reads that belong to the tumour chromosome 
# we don't know tho if the loss of chr1 is subclonal in some tumour samples - what about this?

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 0] # empty  

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 1] # 1
# "chr2_57814739_A_C" # mixed feelings - looks okay but quite a lot of insertions, double check

twins_filtered_mtr[sum_tumour == 8 & sum_normal == 2] # 2
# chr2_78054577_G_A # looks okay
# chr3_153319439_A_G # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 0] # 1  
# chr7_152562130_G_A # looks good 

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 1] # 2
# chr12_84029949_C_A # looks okay
# chr6_77631804_C_A # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 2] # 1
# chr2_158600954_G_A # looks okay

twins_filtered_mtr[sum_tumour == 7 & sum_normal == 3] # 1
# chr8_124634298_G_T # deletions, should have been removed?

# Not sure the below are worth it since we don't know if they tell us anything about the tumour phylogeny or not 
twins_filtered_mtr[sum_tumour == 6 & sum_normal == 0] # 0 

twins_filtered_mtr[sum_tumour == 6 & sum_normal == 1] # 3
# chr10_51734221_A_T # looks okay
# chr12_112652964_C_T # looks okay
# chr15_94939486_G_A # looks okay

twins_filtered_mtr[sum_tumour == 6 & sum_normal == 2] # 1
# chr18_64267234_G_A # looks okay, inspect further - is this on the deleted copy?

twins_filtered_mtr[sum_tumour == 5 & sum_normal == 0] # 3
# chr16_17805832_C_T # looks okay 
# chr7_153784080_C_G # don't think this is real
# chr8_21444725_G_A # looks okay

twins_filtered_mtr[sum_tumour == 5 & sum_normal == 1] # 7
# "chr11_26060345_G_T" # looks okay
# "chr12_129105882_A_G" # poor mapping 
# "chr13_103318073_T_A" # looks okay
# "chr14_54164980_C_G" # looks okay
# "chr17_70958446_C_T" # looks okay
# "chr6_168731130_C_T" # looks okay  
# "chr6_8310684_G_T" # looks okay

twins_filtered_mtr[sum_tumour == 4 & sum_normal == 0] # 1
# chr13_28953051_G_T # looks okay

twins_filtered_mtr[sum_tumour == 4 & sum_normal == 1] # 3
# "chr12_114403258_G_A" # looks okay
# "chr15_23357707_C_T" # really quite poor mapping 
# "chr18_1379467_G_T"  # looks okay
# "chr4_135461656_G_A" # looks okay

twins_filtered_mtr[sum_tumour == 3 & sum_normal == 0] # 5  
# "chr11_43029745_A_T" # looks okay
# "chr12_19555916_C_T" # looks okay
# "chr12_82832044_C_A" # looks okay
# "chr2_231591536_G_T" # looks okay
# "chr5_141433408_A_T" # indels - why not thrown away?

twins_filtered_mtr[sum_tumour == 3 & sum_normal == 1] # 4
# "chr1_145588334_T_C" # poor mapping (not extremely bad)
# "chr1_148952577_G_A" # poor mapping 
# "chr4_179587218_G_A" # looks real!
# "chr8_141860196_G_A" # looks okay

######################################################################################################
# Which PD62341 tumour samples share most with PD63383 tumour samples?
# I want to know which PD62341 tumour sample / clone the PD63383 tumour arose from

twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 1
# "chr12_82832044_C_A" # looks great # shared with am 
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 0
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 1
# chr8_21444725_G_A # looks great # shared with am, b, u
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 0
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 1
# chr7_152562130_G_A # looks great # shared with am, ae, b, aj, ap
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 0

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 3
# chr1_145588334_T_C # looks a bit dodgy 
# chr4_179587218_G_A # looks okay (shared w/ u, 1 read in am)
# chr8_141860196_G_A # looks okay (shared w/ am, 2 reads in u)
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 1
# chr4_135461656_G_A # looks great (b, low levels am, ae)
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 6
# shared w/ b (6), am (5), ap (2), ae (2), u (2)
# chr11_26060345_G_T # looks great
# chr12_129105882_A_G low level in all samples and maps poorly 
# chr13_103318073_T_A # looks great
# chr14_54164980_C_G # looks great
# chr17_70958446_C_T # looks great
# chr6_168731130_C_T # looks great
# seems most sharing is with am and ap
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 3
# typically with am, ap, ae / b / u
# chr10_51734221_A_T # looks okay
# chr12_112652964_C_T # looks okay
# chr15_94939486_G_A # looks okay
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 2
# am, ap, b, u, ae / aj
# chr12_84029949_C_A # looks okay
# chr6_77631804_C_A # looks okay
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 1
# am, ap, b, u, ae, aj
# chr2_57814739_A_C # I think it's fine but double check

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 1
# present in aj, am, ap, b, u
# chr2_158600954_G_A # looks okay
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 2
# present in ae, am, ap, aj, u, b
# chr2_78054577_G_A # looks okay
# chr3_153319439_A_G # looks okay (fun mutation next to this!)

twins_filtered_mtr[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 0
twins_filtered_mtr[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 0
twins_filtered_mtr[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 0
twins_filtered_mtr[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 0
twins_filtered_mtr[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 1
# chr8_124634298_G_T # looks like a lot of insertions 
twins_filtered_mtr[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 2
# chr14_49354012_C_T # looks okay (ae, aj, am, ap, b, u)
# chr5_9088637_A_T # deletions (why not picked up?)

# try doing the same by VAF
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 0
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 1
# chr8_21444725_G_A # looks good (am, b, u)
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 0
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 0
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 1
# chr7_152562130_G_A # looks okay (ae, aj, am, ap, b)
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 0

twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 1
# chr8_141860196_G_A # looks okay (am)
twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 2
# chr12_82832044_C_A (am and b) # looks okay 
# chr6_168731130_C_T (am and b) # looks okay
twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 5
# chr11_26060345_G_T # looks okay
# chr12_112652964_C_T # looks okay
# chr13_103318073_T_A # looks okay
# chr17_70958446_C_T # looks okay
# chr4_135461656_G_A # looks okay
twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 5
# chr10_51734221_A_T # looks okay
# chr12_84029949_C_A # looks okay
# chr14_54164980_C_G # looks okay
# chr15_94939486_G_A # looks okay
# chr6_77631804_C_A # looks okay
twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 0
twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 1
# chr12_106354695_C_T # looks okay

twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 1 
# chr1_145588334_T_C # mapping issues 
twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 1
# chr9_100061581_C_T # double check 
twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 0
twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 1
# chr2_158600954_G_A # looks great
twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 1
twins_filtered_vaf[sum_normal == 2 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 2
# chr14_49354012_C_T # looks fine 
# chr2_57814739_A_C # fine?
# chr2_78054577_G_A # looks okay

twins_filtered_vaf[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 1] # 0
twins_filtered_vaf[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 2] # 0
twins_filtered_vaf[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 3] # 0
twins_filtered_vaf[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 4] # 0
twins_filtered_vaf[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 5] # 0
twins_filtered_vaf[sum_normal == 3 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 6] # 4
# chr2_232306029_G_C # looks okay
# chr5_119778240_G_C # poor mapping, not sure
# chr7_120677593_C_T # looks okay (?)
# chr7_12297469_T_G # looks okay

# get a list of likely subclonal mutations (lower VAF in PD62341) - tho this could be sample purity 
muts_PD63383t_shared = c("chr12_82832044_C_A", "chr8_21444725_G_A", "chr7_152562130_G_A", 'chr2_57814739_A_C', 'chr6_77631804_C_A',
                         'chr12_84029949_C_A', 'chr15_94939486_G_A', 'chr12_112652964_C_T', 'chr11_26060345_G_T', 
                         'chr13_103318073_T_A', 'chr14_54164980_C_G', 'chr17_70958446_C_T', 'chr6_168731130_C_T', 'chr2_158600954_G_A',
                         'chr2_78054577_G_A', 'chr3_153319439_A_G', 'chr14_49354012_C_T', 'chr12_106354695_C_T', 'chr8_141860196_G_A',
                         'chr10_51734221_A_T', 'chr4_135461656_G_A', 'chr2_232306029_G_C', 'chr7_120677593_C_T', 'chr7_12297469_T_G')

muts_tumour_shared = twins_filtered_vaf[mut_ID %in% muts_PD63383t_shared, c(1:23), with=FALSE]

muts_tumour_shared_melt = melt(muts_tumour_shared, id.vars = 'mut_ID')
muts_tumour_shared_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
muts_tumour_shared_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
muts_tumour_shared_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
muts_tumour_shared_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]

ggplot(muts_tumour_shared_melt[status=='tumour'], aes(x=value, y=mut_ID, col=sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF: PD63383 tumour mutations shared with some PD62341'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_muts_shared_PD63383_xy.pdf'), width=9, height=7)

ggplot(muts_tumour_shared_melt[status=='tumour'], aes(x=mut_ID, y=value, col=sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF: PD63383 tumour mutations shared with some PD62341'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))+
  geom_hline(yintercept = 0.5, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_muts_shared_PD63383.pdf'), width=9, height=7)

# add median for each tumour 
muts_tumour_shared_melt2 = data.table(muts_tumour_shared_melt %>% group_by(mut_ID, sample_type) %>% mutate(mean = mean(value)))
muts_tumour_shared_melt2[, mut_ID0 := as.numeric(factor(mut_ID)) - 0.35]
muts_tumour_shared_melt2[, mut_ID1 := as.numeric(factor(mut_ID)) + 0.35]
ggplot(muts_tumour_shared_melt2[status=='tumour'], aes(x=mut_ID, y=value, col=sample_type))+
  geom_point(size=0.8, alpha=0.6)+
  theme_bw(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF')+
  ggtitle(glue('VAF: PD63383 tumour mutations shared with some PD62341'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))+
  geom_hline(yintercept = 0.5, color='black', size = 0.3)+
  geom_segment(aes(x = mut_ID0, xend = mut_ID1, y = mean, yend = mean, color = sample_type), size = 1)+
  scale_color_manual(values = c('tumour_PD62341' = col_tumour_PD62341, 'tumour_PD63383' = col_tumour_PD63383))
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_muts_shared_PD63383.pdf'), width=9, height=7)

######################################################################################################
# Relationship between PD62341 and PD63383 tumours 
# PD63383 is a subclone of PD62341, which likely evolved independetly after initial transmission

######################################################################################################
# Mutations private to PD62341 tumour (subclonal)
# BUT: are these subclonal or are they shared bc the other samples are contaminated?

# searching for mutations acquired after PD63383 tumour migration
# actually, good question: has the tumour migrated once only? or were there several transmissions?

# no muts present in 8/8 PD62341 and none PD63383 tumour
# actually not sure if it makes sense to split PD63383 / PD62341 - if there is nothing in PD63383 tumour, the normal samples couldn't have been contaminated
twins_filtered_mtr[sum_normal_PD62341 == 0 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_mtr[sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0
twins_filtered_mtr[sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 8] # 0

# no muts present in 7/8 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 7] # 0

# muts present in 6/8 PD62341 tumour and none PD63383 tumour (nothing real)
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 6] # 0

# muts present in 5/8 PD62341 and none PD63383 tumour (2 real ones!)
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 2
# "chr16_17805832_C_T" # looks okay 
# ae, ag, aj, b, u (interesting - all non-primary tumour so that would suggest the PD63383 emerged from primary)
# "chr7_153784080_C_G" # not believable, double check
twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 1
# "chr6_8310684_G_T" # looks okay!
# ae, ag, aj, b, u (again!) and present in PD62341h normal - likely contaminated
twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 5] # 0
# empty

# muts present in 4/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4] # 1
# "chr13_28953051_G_T" # looks okay
# present in ae, ag, aj, u (3 reads in b and ap!)

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4] # 2
# "chr12_114403258_G_A" # looks okay
# ae, aj, b, u (2 reads in ag and aa, 1 in am and ap)
# "chr18_1379467_G_T" # looks okay 
# ae, ag, aj, u (3 reads in b, 2 in aa, 1 in v, am, also in h but likely contaminated)

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 4] # 0

# muts present in 3/10 PD62341 and none PD63383 tumour
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3] # 3
# "chr11_43029745_A_T" very low MTR numbers everywhere but looks okay on Jb
# "chr12_19555916_C_T" very low MTR numbers everywhere but looks okay on Jb
# "chr2_231591536_G_T" looks good (present in aj and ae, less in u)

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3] # 3
# "chr1_248574696_C_G" looks okay

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 3, mut_ID] # 0
# "chr10_119182885_G_C" # poor mapping 

# muts present in 2/10 PD62341 and none PD63383 tumour
# not sure these are worth investigating in detail because ultimately each tumour sample would have been seeded by more than one cell probably
# so not very likely that you would identify anything that is truly clonal 
twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2] # 5
# "chr13_62110424_G_A" # looks good 
# "chr1_14644486_C_A" # looks good  
# "chr1_38311094_G_T # looks good 
# "chr1_83136026_A_T" # not great, double check 
# "chrX_9804641_G_A" # looks okay

twins_filtered_mtr[sum_normal == 1 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2] # 4
# "chr11_18922568_G_T" # double check (repetitive region)
# "chr17_1199637_C_G" # poor mapping 
# "chr9_100061581_C_T" # maybe? weird region so may have thrown Illumina off (double check)
# "chr9_67625770_G_A" # some mapping issues 

twins_filtered_mtr[sum_normal == 2 & sum_tumour_PD63383 == 0 & sum_tumour_PD62341 == 2] # 2
# "chr19_50615732_C_T" # double check
# "chr9_39590690_G_A" # poor mapping 

######################################################################################################
# Mutations shared between tumour and one twin 

# Mutations present in ALL tumour samples and all / some normal PD62341 samples 
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 1
# "chr14_105458006_C_A" our favourite one!

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 2
# "chr3_62055057_C_G" aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad
# "chr3_62055077_G_C" aa and h (possibly contamination), n and q (heart and pancreas), missing in v (?) and ad

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 0] # 0

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 1
# chr6_95827754_C_A not the best mapping quality but fine, aa and h (possibly contamination) 

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 4
# "chr16_32621521_C_A" # mapping quality not excellent (double check because it doesn't look that bad)
# "chr17_78256883_C_T" # looks okay
# "chr22_17774668_G_A" # looks okay 
# "chr22_30553442_G_T" # looks okay

# Mutations present in 9 tumour samples and all / some PD62341 normal samples 
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 0
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 0
twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 0] # 1
# chr7_139659050_G_A looks okay
# also there is a very cool mutation on the other chromosome (in all samples?) chr7_139659034_T_G
# missing in ak, which is technically tumour (primary)

twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 0] # 2
# chr18_71485446_G_A # looks okay and on the segment lost in the other chromosome 
# chr21_45729938_T_A # looks okay
# chr2_119170376_C_T # very significant strand bias but technically fine? double check
# all missing in ak 

twins_filtered_mtr[sum_tumour == 9 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 0] # 5
# mostly in h so all of these are likely contamination
# chr12_61020950_C_A # looks okay
# chr7_120677593_C_T # looks okay
# chr8_46678612_T_C # looks okay
# chrX_46375010_C_T # looks okay
# chrX_56448796_C_A # looks okay
# all missing in ak 

# Allow presence in 1 normal sample from PD63383 (because of contamination)
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 1] # 0
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 1] # 3
# chr16_5479739_C_T # looks real (absent from PD62341v, present in PD63383bb)
# chr2_95662131_G_A # looks semi-okay (double check) (absent from PD62341v, present in PD63383bb)
# chr3_50106043_C_T # looks okay (absent from PD62341v, present in PD63383bb)

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 1] # 0
twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 3 & sum_normal_PD63383 == 1] # 1
# chr5_44907911_G_A # looks okay (present in PD62341q, aa, h and PD63383bb)

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 2 & sum_normal_PD63383 == 1] # 9
# always PD62341aa, PD62341h and PD63383bb so very likely contamination 
# "chr10_100754461_C_T" # looks okay
# "chr11_134158143_C_T" # looks okay
# "chr14_104500332_C_T" # looks okay
# "chr15_23462705_C_A" # mis-mapping 
# "chr17_40061856_G_C" # looks okay 
# "chr4_130207996_C_T" # looks okay 
# "chr4_147908994_T_C" # looks okay 
# "chr4_178731458_A_G" # looks okay
# "chrX_66719643_C_G" # looks okay

twins_filtered_mtr[sum_tumour == 10 & sum_normal_PD62341 == 1 & sum_normal_PD63383 == 1] # 10
# always PD62341aa or PD62341h and PD63383bb so very likely contamination 
# "chr13_60216504_A_C" 
# "chr15_48353502_C_A" 
# "chr15_56691722_C_T" 
# "chr16_4003400_G_A" 
# "chr1_69984580_G_A" # looks like it's been deleted in tumour samples 
# "chr2_133311125_G_T" 
# "chr5_136349883_C_T" 
# "chr5_37790015_G_A" 
# "chr9_11364991_T_A"  
# "chrX_68803487_C_T"

######################################################################################################
# Mutations specific to PD63383 tumour ONLY

twins_filtered_mtr[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0] # 12
# "chr10_100308403_C_T" # looks okay
# "chr11_80115755_C_T" # looks okay
# "chr19_43348626_C_A" # not sure, double check (checked on blat and I think it's fine)
# "chr1_51714096_G_C" # looks okay  
# "chr2_72939205_C_T" # looks okay 
# "chr2_82147172_T_C" # looks okay 
# "chr5_157248612_A_G" # looks okay
# "chr5_54017866_G_T" # looks okay
# "chr7_13677323_A_G" # mapping not excellent, double check (I think it's fine but does map to quite a lot of random stuff)
# "chrX_48453603_G_A" # looks okay
# "chrX_72671786_A_T" # poor mapping + strand bias 
# "chrX_77681162_G_A" # looks okay

twins_filtered_mtr[sum_normal_PD63383 == 1  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0] # 15
# "chr10_31333522_G_A" # looks okay # bb
# "chr12_43132316_G_A" # looks okay but close to poly-T # bb
# "chr12_80762031_G_A" # looks okay # bb
# "chr1_160088990_G_A" # looks okay # bb
# "chr22_25281865_C_T" # looks okay # bb
# "chr3_137508691_C_A" # looks okay # bb
# "chr6_86549064_C_T" # looks okay # bb
# "chr6_92179408_A_G" # looks okay # bb
# "chr7_148446425_G_A" # looks okay # bb
# "chr7_49732059_G_A" # looks okay # bb 
# "chr9_98275090_C_T" # looks okay # bb
# "chrX_117418608_G_A" # looks okay # bb
# "chrX_119915655_G_A" # looks okay # bb
# "chrX_36779694_C_T" # looks okay # bb
# "chrX_70352624_G_A" # looks okay # bb

# "chr7_49732059_G_A" note that there is a mutation next to this one 
# chr7_49732077_G_C which is present in more samples
# so can track subclonal evolution in a sense 
# good to check where chr7_49732077_G_C is present (germline?)
# I tried to check but not sure it was called anywhere 

# reads in normal samples are exclusively present in bb
# PD63383bb is skin, so I find it highly plausible that this is contamination with tumour cells
# especially that the kid had lesions in the skin - can I double check with Henry what he thinks?

muts_PD63383_tumour2 = twins_filtered_mtr[sum_normal_PD63383 == 1  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist()

twins_filtered_mtr[sum_normal_PD63383 == 2  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty
twins_filtered_mtr[sum_normal_PD63383 == 3  & sum_normal_PD62341 == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty 
twins_filtered_mtr[sum_normal_PD63383 == 0  & sum_normal_PD62341 == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0, mut_ID] # 0
# empty 

# search by VAF
twins_filtered_vaf[sum_normal == 0 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0] # 10
# "chr10_100308403_C_T" # looks okay
# "chr11_80115755_C_T" # looks okay
# "chr19_43348626_C_A" # not sure, double check (checked on blat and I think it's fine)
# "chr1_51714096_G_C" # looks okay  
# "chr2_72939205_C_T" # looks okay 
# "chr2_82147172_T_C" # looks okay 
# "chr5_54017866_G_T" # looks okay
# "chr7_13677323_A_G" # mapping not excellent, double check (I think it's fine but does map to quite a lot of random stuff)
# "chrX_48453603_G_A" # looks okay
# "chrX_72671786_A_T" # poor mapping + strand bias 

twins_filtered_vaf[sum_normal == 1 & sum_tumour_PD63383 == 2 & sum_tumour_PD62341 == 0] # 17
# "chr10_31333522_G_A" # looks okay # bb
# "chr12_43132316_G_A" # looks okay but close to poly-T # bb
# "chr12_80762031_G_A" # looks okay # bb
# "chr1_160088990_G_A" # looks okay # bb
# "chr22_25281865_C_T" # looks okay # bb
# "chr3_137508691_C_A" # looks okay # bb
# "chr4_179587218_G_A" # looks okay (4 reads in PD62341u, VAF 0.0851 - that's why not called in MTR - ??)
# "chr6_86549064_C_T" # looks okay # bb
# "chr6_92179408_A_G" # looks okay # bb
# "chr7_148446425_G_A" # looks okay # bb
# "chr7_49732059_G_A" # looks okay # bb 
# "chr9_98275090_C_T" # looks okay # bb
# "chrX_117418608_G_A" # looks okay # bb
# "chrX_119915655_G_A" # looks okay # bb
# "chrX_36779694_C_T" # looks okay # bb
# "chrX_70352624_G_A" # looks okay # bb
# "chrX_77681162_G_A" # looks okay

muts_PD63383_tumour_mtr = c("chr10_100308403_C_T", "chr11_80115755_C_T",  "chr19_43348626_C_A",  "chr1_51714096_G_C",   "chr2_72939205_C_T",  
                            "chr2_82147172_T_C",   "chr5_157248612_A_G",  "chr5_54017866_G_T",   "chr7_13677323_A_G",   "chrX_48453603_G_A",  
                            "chrX_77681162_G_A", "chr10_31333522_G_A", "chr12_43132316_G_A", "chr12_80762031_G_A", "chr1_160088990_G_A", 
                            "chr22_25281865_C_T", "chr3_137508691_C_A", "chr6_86549064_C_T",  "chr6_92179408_A_G",  "chr7_148446425_G_A", 
                            "chr7_49732059_G_A", "chr9_98275090_C_T",  "chrX_117418608_G_A", "chrX_119915655_G_A", "chrX_36779694_C_T",  
                            "chrX_70352624_G_A")
muts_PD63383_tumour_vaf = c("chr10_100308403_C_T", "chr11_80115755_C_T",  "chr19_43348626_C_A",  "chr1_51714096_G_C",   "chr2_72939205_C_T",  
                            "chr2_82147172_T_C",   "chr5_157248612_A_G",  "chr7_13677323_A_G",   "chrX_48453603_G_A",  "chr10_31333522_G_A", 
                            "chr12_43132316_G_A", "chr12_80762031_G_A", "chr1_160088990_G_A", "chr22_25281865_C_T", "chr3_137508691_C_A", 
                            "chr4_179587218_G_A", "chr6_86549064_C_T",  "chr6_92179408_A_G",  "chr7_148446425_G_A", "chr7_49732059_G_A", 
                            "chr9_98275090_C_T",  "chrX_117418608_G_A", "chrX_119915655_G_A", "chrX_36779694_C_T", "chrX_70352624_G_A", 
                            "chrX_77681162_G_A")
muts_PD63383_tumour = Reduce(intersect,  list(muts_PD63383_tumour_mtr, muts_PD63383_tumour_vaf))

twins_info = twins_dt[, c('mut_ID', 'Gene', 'Transcript', 'RNA', 'CDS', 'Protein', 'Type', 'Effect'), with=FALSE]
twins_info[mut_ID %in% muts_PD63383_tumour] # nothing in protein coding genes
paste('Number of mutations private to PD63383 tumour:', length(muts_PD63383_tumour)) # 25 

# plots with PD63383-tumour-specific mutations
muts_PD63383t_dt =  twins_filtered_vaf[mut_ID %in% muts_PD63383_tumour, 1:23, with=FALSE]
muts_PD63383t_melt = melt(muts_PD63383t_dt, id.vars = 'mut_ID')
muts_PD63383t_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
muts_PD63383t_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
muts_PD63383t_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
muts_PD63383t_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]

# plot (mutations on the y axis)
ggplot(muts_PD63383t_melt[status=='tumour'], aes(x=value, y=mut_ID, color=sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF: PD63383 tumour-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_muts_tumour.pdf'), width=9, height=7)

# plot (mutations on the y axis)
ggplot(muts_PD63383t_melt[twin=='PD63383'], aes(x=value, y=mut_ID, color=sample_type))+
  geom_point(size=2.5)+
  scale_color_manual(values = c(col_normal_PD63383, col_tumour_PD63383))+
  theme_bw(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF: PD63383 tumour-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_muts_PD63383.pdf'), width=9, height=7)

# try to estimate contamination in PD63383bb (skin sample)
muts_PD63383t_spec = twins_filtered_dt[mut_ID %in% muts_PD63383_tumour, c('mut_ID', samples_mtr, samples_dep), with=FALSE]
muts_PD63383t_spec[, sum_PD63383t_mtr := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
muts_PD63383t_spec[, sum_PD63383t_dep := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
muts_PD63383t_spec[, mean_PD63383t_vaf := sum_PD63383t_mtr / sum_PD63383t_dep] 
muts_PD63383t_spec[, mean_PD63383bb_vaf := PD63383bb_MTR / PD63383bb_DEP] 

paste('Mean VAF of PD63383 tumour-restricted muts in PD63383 tumour:', mean(muts_PD63383t_spec[, mean_PD63383t_vaf] %>% unlist())) # 0.385
paste('Mean VAF of PD63383 tumour-restricted muts in PD63383 skin:', mean(muts_PD63383t_spec[, mean_PD63383bb_vaf] %>% unlist())) # 0.122
paste('Estimated purity of PD63383 tumour:', 2 * mean(muts_PD63383t_spec[, mean_PD63383t_vaf] %>% unlist())) # 0.770
paste('Estimated purity of PD63383 skin:', 1 - 2 * mean(muts_PD63383t_spec[, mean_PD63383bb_vaf] %>% unlist())) # 0.755

# is either PD63383 tumour sample more contaminated?
muts_PD63383t_spec[, mean_PD63383aq_vaf := PD63383aq_MTR / PD63383aq_DEP] 
muts_PD63383t_spec[, mean_PD63383ap_vaf := PD63383ap_MTR / PD63383ap_DEP] 

paste('Mean VAF of PD63383 tumour-restricted muts in PD63383ap:', mean(muts_PD63383t_spec[, mean_PD63383ap_vaf] %>% unlist())) # 0.412
paste('Mean VAF of PD63383 tumour-restricted muts in PD63383aq:', mean(muts_PD63383t_spec[, mean_PD63383aq_vaf] %>% unlist())) # 0.355
paste('Estimated purity of PD63383ap:', 2 * mean(muts_PD63383t_spec[, mean_PD63383ap_vaf] %>% unlist())) # 0.82
paste('Estimated purity of PD63383aq:', 2 * mean(muts_PD63383t_spec[, mean_PD63383aq_vaf] %>% unlist())) # 0.71

# plot (mutations on the y axis): only PD62341ap, PD63383aq and PD63383bb
ggplot(muts_PD63383t_melt[sample %in% c('PD63383aq', 'PD63383ap', 'PD63383bb')], aes(x=value, y=mut_ID, col=sample))+
  geom_point(size=2.5)+
  theme_bw(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle(glue('VAF: PD63383 tumour-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlim(c(0, 0.7))+
  geom_vline(xintercept = 0.5, color='black', size = 0.3)
ggsave(glue('Results/20241030_p3_vaf_dist_PD63383_muts_PD63383tumour_skin.pdf'), width=9, height=7)

######################################################################################################
######################################################################################################
# Another thing I could do is to just treat all tumour samples as a single sample
twins_filtered_mtr[, sum_reads_mtr := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_filtered_dep[, sum_reads_dep := rowSums(.SD), .SDcols = samples_tumour_dep]

tumour_samples_sums = merge(twins_filtered_mtr[, c('mut_ID', 'sum_reads_mtr'), with=FALSE],
                            twins_filtered_dep[, c('mut_ID', 'sum_reads_dep'), with=FALSE], 
                            by = 'mut_ID')
tumour_samples_sums[, tumour_vaf := sum_reads_mtr / sum_reads_dep]

pdf('Results/20241102_p1_hist_vaf_tumour.pdf')
hist(tumour_samples_sums[, tumour_vaf], xlab = 'VAF', main = 'VAF across tumour samples', 
     xlim = c(0, 1))
abline(v = median(tumour_samples_sums[, tumour_vaf]), col = 'blue')
abline(v = mean(tumour_samples_sums[, tumour_vaf]), col = 'red')
dev.off()

pdf('Results/20241102_p1_hist_vaf_tumour_breaks50.pdf')
hist(tumour_samples_sums[, tumour_vaf], xlab = 'VAF', main = 'VAF across tumour samples', 
     xlim = c(0, 1), breaks = 50)
abline(v = median(tumour_samples_sums[, tumour_vaf]), col = 'blue')
abline(v = mean(tumour_samples_sums[, tumour_vaf]), col = 'red')
dev.off()

pdf('Results/20241102_p1_hist_vaf_tumour_24muts.pdf')
hist(tumour_samples_sums[mut_ID %in% muts_tumour_all_only, tumour_vaf] %>% unlist(), 
     xlab = 'VAF', main = 'VAF across tumour samples', xlim = c(0, 1))
abline(v = median(tumour_samples_sums[mut_ID %in% muts_tumour_all_only, tumour_vaf] %>% unlist()), col = 'blue')
abline(v = mean(tumour_samples_sums[mut_ID %in% muts_tumour_all_only, tumour_vaf] %>% unlist()), col = 'red')
dev.off()

tumour_samples_sums[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
tumour_samples_sums[, muts_all_samples := as.factor(as.numeric(mut_ID %in% muts_tumour_all_only))]

ggplot(tumour_samples_sums, aes(x = sum_reads_dep, y = tumour_vaf, col=Chrom))+
  geom_point(size=1.5)+
  theme_bw(base_size = 14)+
  labs(x = 'Coverage (sum across tumour samples)', y = 'VAF (across tumour samples)')+
  ggtitle(glue('VAF vs coverage across tumour samples'))+
  ylim(c(0, 1))+
  geom_hline(yintercept = 0.35, color='black', size = 0.5)
ggsave(glue('Results/20241102_p2_vaf_vs_coverage_tumour.pdf'))

ggplot(tumour_samples_sums, aes(x = sum_reads_dep, y = tumour_vaf, col=muts_all_samples))+
  geom_point(size=1.5)+
  theme_bw(base_size = 14)+
  labs(x = 'Coverage (sum across tumour samples)', y = 'VAF (across tumour samples)',
       col = '24 muts\n(all tumour samples)')+
  ggtitle(glue('VAF vs coverage across tumour samples'))+
  ylim(c(0, 1))+
  geom_hline(yintercept = 0.35, color='black', size = 0.5)
ggsave(glue('Results/20241102_p2_vaf_vs_coverage_tumour_col24.pdf'))


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# Having some fun with driver analysis 

# are there any mutations protein-coding regions identified?
twins_tumour_info = merge(twins_filtered_vaf[sum_tumour >= 6 & sum_normal <= 3, c('mut_ID', samples_vaf, 
                                        'sum_normal_PD62341', 'sum_normal_PD63383', 'sum_tumour_PD62341', 'sum_tumour_PD63383'), with=FALSE],
                          twins_filtered_dt[, c('mut_ID', cols_info), with=FALSE], by = 'mut_ID')
paste('Number of mutations present in most tumour samples:', dim(twins_tumour_info)[1]) # 81

# should I only include mutations present in the tumour and not normal cells?
# I require 6 tumour samples to carry the mutation (in case of drop-outs etc) 
# I require less than 4 normal samples (as there can be contamination, especially bb, aa, h)
# There could be evolution from the somatic pre-cancerous tissue (we don't have adrenal though)
table_effects = table(twins_tumour_info[Effect!='-', Effect])
genes_missense = twins_tumour_info[Effect=='missense', Gene] %>% unlist() # 0
genes_upstream = twins_tumour_info[Effect=='upstream', Gene] %>% unlist() # 3
genes_downstream = twins_tumour_info[Effect=='downstream', Gene] %>% unlist() # 3
genes_intronic = twins_tumour_info[Effect=='intronic', Gene] %>% unlist() # 0
genes_splice = twins_tumour_info[Effect=='splice_region', Gene] %>% unlist() # 0

# check that the mutations are of good quality 
twins_tumour_info[Gene %in% genes_upstream]
# chr12_106354695_C_T # looks great # POLR3B
# chr17_40061856_G_C # looks great # THRA
# chr22_17774668_G_A # looks okay # BID

twins_tumour_info[Gene %in% genes_downstream]
# chr17_5029328_G_C # looks okay # KIF1C
# chr22_30553442_G_T # looks okay # GAL3ST1
# chr9_12824689_G_A # looks okay # LURAP1L

twins_tumour_info[Gene %in% genes_intronic] # 24
# chr10_100754461_C_T # looks okay
# chr10_51734221_A_T # looks okay
# chr11_134158143_C_T # looks okay
# chr13_100405282_G_T # looks okay
# chr14_104500332_C_T # looks okay
# chr15_56691722_C_T # looks okay
# chr15_76646556_T_C # looks okay
# chr16_4003400_G_A # looks okay
# chr1_69984580_G_A # can estimate contamination from this 
# chr20_44114996_C_T # looks okay
# chr21_45729938_T_A # looks okay
# chr2_124573116_C_T # looks okay
# chr2_133311125_G_T # looks okay
# chr2_1420635_A_T # looks okay
# chr2_158600954_G_A # looks okay
# chr2_2023708_C_G # looks okay
# chr2_232306029_G_C # looks okay
# chr4_147908994_T_C # looks okay
# chr4_46067117_A_G # looks okay
# chr4_93081874_T_A # looks okay
# chr5_136349883_C_T # looks okay
# chr7_120677593_C_T # looks okay
# chr7_139659050_G_A # looks okay
# chr8_97785584_T_C # looks okay
# chrX_124531709_C_T # looks okay

######################################################################################################
# Do I have mutations subclonal in all samples?
twins_filtered_vaf[, sum_subclonal := rowSums(.SD < 0.3 & .SD > 0.1), .SDcols = samples_vaf]
twins_filtered_vaf[sum_subclonal >= 19] # 18
# does it matter if these mutations are also present in normal samples?
# chr11_4220549_T_C # variable, double check
# chr11_4220746_T_C # variable, double check
# chr11_4594537_A_G # some mapping issues 
# chr12_44087006_G_A # not convinced
# chr19_40853394_C_A # poor mapping 
# chr19_50082555_A_C # seems fine (nice other mutation next to this) 
# chr19_50082890_A_G # looks okay
# chr3_190409217_C_T # looks okay (very interesting mut next to this)
# chr4_80929402_G_A # okay that looks great (low in w and bb)
# chr5_107405920_C_A # looks fine (interesting mutation next to this)
# chr5_176771316_G_A # looks okay (low in n)
# chr5_176771327_T_A # looks okay (low in n)
# chr6_94069453_A_T # insertions
# chr9_39582303_C_T # mapping issues 
# chr9_42329754_C_T # mapping issues (double check)
# chrX_35787898_C_T # looks okay (fun mutation next to this one)
# chrX_798031_G_A # looks okay (low in u)
# chrX_83728142_A_G # mapping issues 

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# Plots attempt no 1

# are there mutations shared with only a specific normal tissue of PD62341?
twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0] # 20
# "chr12_61020950_C_A"  # looks okay 
# "chr14_105458006_C_A" # looks okay
# "chr16_32621521_C_A"  # poor mapping 
# "chr17_78256883_C_T" # looks okay
# "chr18_64267234_G_A" # some mapping issues but generally seems okay
# "chr18_71485446_G_A" # looks okay 
# "chr21_45729938_T_A" # looks okay
# "chr22_17774668_G_A" # looks quite okay 
# "chr22_30553442_G_T" # generally looks okay
# "chr2_57814739_A_C" # poor mapping quality   
# "chr3_62055057_C_G" # looks okay  
# "chr3_62055077_G_C" # looks okay
# "chr6_58221412_T_C" # muts only on reads with poor mapping   
# "chr6_95827754_C_A" # rather poor mapping quality  
# "chr7_120677593_C_T" # insertions
# "chr7_139659050_G_A" # looks okay 
# "chr8_46678612_T_C" # looks okay
# "chrX_46375010_C_T" # looks okay  
# "chrX_56448796_C_A" # looks okay (but there may be some insertions)

# count samples with mutations (for all mutations)
colSums(twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0,c(samples_normal_PD62341_mtr), with=FALSE] >= 4)
# v: 5, q: 9, aa: 12, ad: 6, h: 20, n: 5
colSums(twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 == 0,c(samples_normal_PD62341_mtr), with=FALSE] >= 10) # less likely to be contamination
# v: 0, q: 3, aa: 4, ad: 1, h: 4, n: 3

# run analysis for mutations with confirmed high quality 
muts_hq = c("chr12_61020950_C_A", "chr14_105458006_C_A", "chr17_78256883_C_T", "chr18_71485446_G_A",
            "chr19_3961910_G_A", "chr21_45729938_T_A", "chr22_17774668_G_A", "chr3_62055057_C_G",
            "chr3_62055077_G_C", "chr7_139659050_G_A", "chr8_46678612_T_C", "chrX_46375010_C_T", "chrX_56448796_C_A")
paste('Number of high-quality mutations shared between tumour and PD62341 normal:', length(muts_hq))

colSums(twins_filtered_mtr[mut_ID %in% muts_hq, c(samples_normal_PD62341_mtr), with=FALSE] >= 4)
# v: 1, q: 6, aa: 5, ad: 1, h: 12, n: 3
# the one that is not 'present' in h is there at 3 reads, and 5 reads in q (contamination?)
# h is at quite low depth there, so VAF is actually 0.1875 

colSums(twins_filtered_mtr[mut_ID %in% muts_hq, c(samples_normal_PD62341_mtr), with=FALSE] >= 10) # less likely to be contamination
# v: 0, q: 2, aa: 3, ad: 1, h: 3, n: 3

# we can also check VAF (> 0.1) to account for differences in coverage
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_normal_PD62341_vaf), with=FALSE] >= 0.1)
# v: 1, q: 6, aa: 4, ad: 1, h: 13, n: 4

# how does this compare to tumour?
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_tumour_vaf), with=FALSE] >= 0.1)
# all are 13 but PD62341ak (maybe contaminated?? some other cells? or duplicated chromosomes so lower VAF?? - but then this wouldn't drop it below 0.1)

# but then it is very rare that any of these mutations seem clonal
# I imagine this can be because of contamination with other cells (which is not known)
# for normal samples, bear in mind they are not monoclonal so really hard to tell
# also if there are really precursors you could have them at lower levels / freq so lower VAF
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_tumour_vaf), with=FALSE] >= 0.4)
colSums(twins_vaf[mut_ID %in% muts_hq, c(samples_normal_PD62341_vaf), with=FALSE] >= 0.4)
# 1 case in PD62341n and PD62341d 

# Examine VAF relationship 
mut_shared_tumour_normal = twins_vaf[mut_ID %in% muts_hq, c('mut_ID', samples_tumour_vaf, samples_normal_PD62341_vaf), with=FALSE]
mut_shared_melt = melt(mut_shared_tumour_normal, id.vars = 'mut_ID')
mut_shared_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_shared_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_shared_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_shared_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_shared_melt[, mut_ID := factor(mut_ID, levels = c("chr14_105458006_C_A", "chr3_62055057_C_G", "chr3_62055077_G_C",
                                                      "chr8_46678612_T_C", "chr7_139659050_G_A", "chr18_71485446_G_A", 
                                                      "chr12_61020950_C_A", "chr17_78256883_C_T", "chr19_3961910_G_A", 
                                                      "chr21_45729938_T_A", "chr22_17774668_G_A", "chrX_46375010_C_T", "chrX_56448796_C_A"))] 

ggplot(mut_shared_melt, aes(x = value, y = mut_ID, col = status))+
  geom_point(size=2, position=position_dodge(.75))+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle('VAF in tumour vs normal')
ggsave('Results/20241020_p2_vaf_dist_for_each_shared_mut13.pdf', width=9, height=7)

# so would we conclude that bc most mutations are shared, the tumour arose in the liver?
# is VAF higher in the liver than other tissues? for mutations which are shared?
mut_shared_melt[, sample_name := factor(fcase(
  sample == 'PD62341h', 'liver',
  sample == 'PD62341q', 'pancreas',
  sample == 'PD62341n', 'heart',
  sample == 'PD62341v', 'spleen',
  sample == 'PD62341ad', 'cerebellum',
  sample == 'PD62341aa', 'skin'
))]
mut_shared_melt[, sample_name := factor(sample_name, 
                                        levels=c('liver', 'pancreas', 'heart', 'spleen', 'cerebellum', 'skin'))]

ggplot(mut_shared_melt[status=='normal'], aes(x = value, y = mut_ID, col = sample_name))+
  geom_point(size=3)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c('#af0101', '#cc7a3d', '#e6c994',
                                '#ffaa6a', '#74a892', '#004367'))+
  labs(x = 'VAF', y = 'Mutation', color = 'Sample')+
  ggtitle('VAF distribution in normal samples')
ggsave('Results/20241020_p2_vaf_dist_for_each_shared_mut13_normalsamples.pdf', width=10, height=7)

# can I see mutations only present in the liver of PD62341 and the tumour?
# there are 6: there is 1 present in PD62341q (5 reads, 3 reads in PD62341h so can be due to coverage differences)
twins_filtered_mtr[mut_ID %in% muts_hq & sum_tumour >= 6 & sum_normal_PD62341 == 1]
colSums(twins_filtered_mtr[mut_ID %in% muts_hq, c(samples_normal_PD62341_mtr), with=FALSE] >= 10) # less likely to be contamination

# note that a real mutation could be present in PD63383 due to contamination
twins_filtered_mtr[sum_tumour >= 6 & sum_normal_PD62341 >= 1 & sum_normal_PD63383 <= 1] # 108 

# better plots for VAF
mut_shared_melt[order(sample, value)]
mut_shared_melt[, sample := as.factor(sample)]

for (s in mut_shared_melt[sample_type=='normal_PD62341',sample] %>% unlist() %>% unique()){
  ggplot(mut_shared_melt[sample==s], aes(x=mut_ID, y=value))+
    geom_line(group = 1)+
    geom_point(size=1.5)+
    geom_area(aes(fill = mut_ID))+
    scale_fill_manual(values = col_normal_PD62341)+
    theme_bw(base_size = 10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x = 'Mutation', y = 'VAF')+
    ggtitle(glue('VAF distribution in sample {s}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ylim(c(0, 0.6))
  ggsave(glue('Results/20241020_p2_vaf_dist_13mut_normal_{s}.pdf'), width = 6, height = 2.5)
}

# look at each mutation one by one 
mut_hq13 = twins_vaf[mut_ID %in% muts_hq, c('mut_ID', samples_tumour_vaf, samples_normal_vaf), with=FALSE]
mut_hq13[, sum_all := rowSums(.SD >= 0.1), .SDcols = 2:23]
mut_hq13[, sum_normal := rowSums(.SD >= 0.1), .SDcols = samples_normal_vaf]
mut_hq13[, sum_tumour := rowSums(.SD >= 0.1), .SDcols = samples_tumour_vaf]
mut_hq13[, sum_normal_PD62341 := rowSums(.SD >= 0.1), .SDcols = samples_normal_PD62341_vaf]
mut_hq13[, sum_tumour_PD62341 := rowSums(.SD >= 0.1), .SDcols = samples_tumour_PD62341_vaf]
mut_hq13[, sum_normal_PD63383 := rowSums(.SD >= 0.1), .SDcols = samples_normal_PD63383_vaf]
mut_hq13[, sum_tumour_PD63383 := rowSums(.SD >= 0.1), .SDcols = samples_tumour_PD63383_vaf]
mut_hq13[order(sum_all)] # note that there are reads in PD63383bb but these are likely contamination (skin!!)

######################################################################################################
# Visualization: heatmaps  

# order the columns such that tumour and normal samples are separated and ordered by # of mutations present
mut_hq13_sub = mut_hq13[, 1:17]
setcolorder(mut_hq13_sub, c('mut_ID', 'PD63383aq_VAF', 'PD63383ap_VAF',
                            'PD62341ae_VAF', 'PD62341ag_VAF', 'PD62341aj_VAF', 'PD62341am_VAF', 'PD62341ap_VAF', 'PD62341b_VAF', 'PD62341u_VAF', 'PD62341ak_VAF',
                            'PD62341h_VAF', 'PD62341q_VAF', 'PD62341n_VAF', 'PD62341aa_VAF', 'PD62341ad_VAF', 'PD62341v_VAF'))
row_order = c(2, 8, 9, 11, 10, 4, 1, 3, 5, 12, 13)
mut_hq13_sub = mut_hq13_sub[row_order,]

# matrix with VAF data (not binary)
mat_vaf = as.matrix(mut_hq13_sub[,2:17])
rownames(mat_vaf) =  mut_hq13_sub[,1] %>% unlist()  
colnames(mat_vaf) = c('PD63383aq', 'PD63383ap',
                      'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                      'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_vaf)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_vaf.pdf')
pheatmap(mat_vaf,
         cellwidth=14, cellheight=14,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations", 
         legend = T,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

# matrix with binary data (VAF > 0.1 - present, else - absent)
mut_hq13_sub_binary_vaf01 = data.table(mut_hq13_sub)

for (j in names(mut_hq13_sub_binary_vaf01)[2:17]){
  
  rows_to_0 = which(mut_hq13_sub_binary_vaf01[[j]] < 0.1)
  set(mut_hq13_sub_binary_vaf01, rows_to_0, j, 0)
  
  rows_to_1 = which(mut_hq13_sub_binary_vaf01[[j]] >= 0.1)
  set(mut_hq13_sub_binary_vaf01, rows_to_1, j, 1)
  
}

mat_vaf01 = as.matrix(mut_hq13_sub_binary_vaf01[,2:17])
rownames(mat_vaf01) =  mut_hq13_sub_binary_vaf01[,1] %>% unlist()  
colnames(mat_vaf01) = c('PD63383aq', 'PD63383ap',
                        'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                        'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_vaf01)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_vaf_binary.pdf')
pheatmap(mat_vaf01,
         cellwidth=14, cellheight=14,
         color = c('#597ad3', '#e76f04'),
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations (VAF > 0.1)", 
         legend = F, legend_breaks = c(-0.5, 0.5, 1),
         cluster_rows=F, cluster_cols=F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

# repeat for MTR >= 4
mut_hq13_mtr = twins_mtr[mut_ID %in% muts_hq, c('mut_ID', samples_tumour_mtr, samples_normal_PD62341_mtr), with=FALSE]
setcolorder(mut_hq13_mtr, c('mut_ID', 'PD63383aq_MTR', 'PD63383ap_MTR',
                            'PD62341ae_MTR', 'PD62341ag_MTR', 'PD62341aj_MTR', 'PD62341am_MTR', 'PD62341ap_MTR', 'PD62341b_MTR', 'PD62341u_MTR', 'PD62341ak_MTR',
                            'PD62341h_MTR', 'PD62341q_MTR', 'PD62341n_MTR', 'PD62341aa_MTR', 'PD62341ad_MTR', 'PD62341v_MTR'))
row_order = c(2, 8, 9, 11, 10, 4, 1, 3, 5, 12, 13)
mut_hq13_mtr = mut_hq13_mtr[row_order,]

# matrix with VAF data (not binary)
mat_mtr = as.matrix(mut_hq13_mtr[,2:17])
rownames(mat_mtr) =  mut_hq13_mtr[,1] %>% unlist()  
colnames(mat_mtr) = c('PD63383aq', 'PD63383ap',
                      'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                      'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_mtr)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_mtr.pdf')
pheatmap(mat_mtr,
         cellwidth=14, cellheight=14,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations", 
         legend = T,
         cluster_rows=F, cluster_cols = F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

# matrix with binary data (MTR >= 4)
mut_hq13_binary_mtr4 = data.table(mut_hq13_mtr)

for (j in names(mut_hq13_binary_mtr4)[2:17]){
  
  rows_to_0 = which(mut_hq13_binary_mtr4[[j]] < 4)
  set(mut_hq13_binary_mtr4, rows_to_0, j, 0)
  
  rows_to_1 = which(mut_hq13_binary_mtr4[[j]] >= 4)
  set(mut_hq13_binary_mtr4, rows_to_1, j, 1)
  
}

mat_mtr4 = as.matrix(mut_hq13_binary_mtr4[,2:17])
rownames(mat_mtr4) =  mut_hq13_binary_mtr4[,1] %>% unlist()  
colnames(mat_mtr4) = c('PD63383aq', 'PD63383ap',
                       'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u', 'PD62341ak',
                       'PD62341h', 'PD62341q', 'PD62341n', 'PD62341aa', 'PD62341ad', 'PD62341v')

# annotation of column names
col_annotation = data.frame(Status = c(rep('tumour', 10), rep('normal', 6)))
rownames(col_annotation) = colnames(mat_mtr4)
annotation_colors = list(
  Status = c(normal=col_normal, tumour=col_tumour))

pdf('Results/20241020_p2_heatmap_status_13mut_mtr_binary.pdf')
pheatmap(mat_mtr4,
         cellwidth=14, cellheight=14,
         color = c('#597ad3', '#e76f04'),
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Status of shared mutations (MTR >= 4)", 
         legend = F, legend_breaks = c(-0.5, 0.5, 1),
         cluster_rows=F, cluster_cols=F, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=6) 
dev.off()

