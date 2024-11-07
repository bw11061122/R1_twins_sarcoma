# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-11-06
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
muts = read.table('Data/mutations_include_20241106_1002.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 1002
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Add column to indicate chromosomes lost in the tumour
twins_filtered_dt[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'
))]

# Import dataframe with purity estimates
purity_dt = data.table(read.csv('Data/20241106_estimates_tumour_cont_68muts.csv'))
purity_dt

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

######################################################################################################
# Basic checks: number of mutations

# Create lists of mutations 
muts_all = twins_filtered_vaf[, mut_ID] %>% unlist()
muts_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist())) 
muts_normal = Reduce(intersect, list(twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist())) 
muts_tumour = Reduce(intersect, list(twins_filtered_vaf[sum_tumour>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist()))
muts_normal_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()))
muts_tumour_all = Reduce(intersect, list(twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()))
muts_tumour_only = Reduce(intersect, list(twins_filtered_vaf[sum_normal<4, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal<4, mut_ID] %>% unlist()))
muts_tumour_all_only = Reduce(intersect, list(muts_tumour_all, muts_tumour_only)) %>% unlist()
muts_normal_only = Reduce(intersect, list(twins_filtered_vaf[sum_tumour<2, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour<2, mut_ID] %>% unlist()))
muts_normal_all_only = Reduce(intersect, list(muts_normal_all, muts_normal_only)) %>% unlist()

paste('Mutations present in all samples:', length(muts_normal)) # 575
paste('Number of muts present in min 1 normal sample:', length(muts_normal)) # 963
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour)) # 981
paste('Number of muts present in all normal samples:', length(muts_normal_all)) # 632
paste('Number of muts present in all tumour samples:', length(muts_tumour_all)) # 641
paste('Number of muts present in only tumour samples:', length(muts_tumour_only)) # 155
paste('Number of muts present in all tumour samples and max 3 normal:', length(muts_tumour_all_only)) # 24
paste('Number of muts present in only normal samples:', length(muts_normal_only)) # 19
paste('Number of muts present in all normal samples and max 1 tumour:', length(muts_normal_all_only)) # 0 # makes sense since tumour is from normal

######################################################################################################
# Analysis of tumour evolution: coverage vs VAF plot

twins_filtered_mtr[, sum_tumour_mtr := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_filtered_dep[, sum_tumour_dep := rowSums(.SD), .SDcols = samples_tumour_dep]

# Aggregate data for all tumour samples together 
twins_tumour_agg = merge(twins_filtered_mtr[, c('mut_ID', 'sum_tumour_mtr'), with=FALSE], 
                         twins_filtered_dep[, c('mut_ID', 'sum_tumour_dep'), with=FALSE])
twins_tumour_agg[, tumour_vaf := sum_tumour_mtr / sum_tumour_dep]
twins_tumour_agg[, mut_class := as.factor(fcase(
  mut_ID %in% muts_all, 'all samples',
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  mut_ID %in% muts_normal_only, 'only normal samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_all, muts_normal_only, muts_tumour_all_only, muts_tumour_only), 'other'))]

# Plot coverage vs VAF 
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_1002.pdf'))

ggplot(twins_tumour_agg[mut_class!='all samples'], aes(x = sum_tumour_dep, y = tumour_vaf, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_not_all_1002.pdf'))

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
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241106_cov_vs_vaf_mut_type.pdf'))

# color by chromosome 
twins_tumour_agg[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = Chrom))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241106_cov_vs_vaf_chrom.pdf'))

twins_tumour_agg[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'))]
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = loss))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241106_cov_vs_vaf_loss.pdf'))

# add number of tumour samples the mutation is present in 
twins_tumour_agg = merge(twins_tumour_agg, twins_filtered_vaf[, c('mut_ID', 'sum_tumour'), with=FALSE], by = 'mut_ID')
twins_tumour_agg[, sum_tumour := as.factor(sum_tumour)]
ggplot(twins_tumour_agg, aes(x = sum_tumour_dep, y = tumour_vaf, col = sum_tumour))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = '# tumour samples\nmut detected in')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241106_cov_vs_vaf_sum_tumour.pdf'))

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
ggsave(glue('Results/20241106_vaf_across_the_genome_1002.pdf'))

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
pdf('Results/20241106_p3_hist_vaf_dist_630muts.pdf')
hist(twins_tumour_agg[mut_ID %in% muts_mapped_tumour, tumour_vaf], breaks = 20,
     xlab = 'VAF (total tumour)', main = '630 mutations', xlim = c(0, 1))
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
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_630mapped.pdf'))

# Aggregate data for all tumour samples together 
twins_filtered_mtr[, sum_PD62341_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dep[, sum_PD62341_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_filtered_mtr[, sum_PD63383_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_filtered_dep[, sum_PD63383_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]

twins_PD62341_agg = merge(twins_filtered_mtr[, c('mut_ID', 'sum_PD62341_mtr'), with=FALSE], 
                         twins_filtered_dep[, c('mut_ID', 'sum_PD62341_dep'), with=FALSE], by = 'mut_ID')
twins_PD62341_agg[, PD62341_vaf := sum_PD62341_mtr / sum_PD62341_dep]

twins_PD63383_agg = merge(twins_filtered_mtr[, c('mut_ID', 'sum_PD63383_mtr'), with=FALSE], 
                          twins_filtered_dep[, c('mut_ID', 'sum_PD63383_dep'), with=FALSE], by = 'mut_ID')
twins_PD63383_agg[, PD63383_vaf := sum_PD63383_mtr / sum_PD63383_dep]

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
ggsave(glue('Results/20241106_vaf_tumour_vs_PD62341_630mapped.pdf'))

ggplot(twins_agg_filt, aes(x=tumour_vaf, y=PD63383_vaf))+
  geom_point()+
  theme_classic(base_size=15)+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  labs(x = 'VAF (total tumour)', y = glue('VAF in PD63383'), col = 'Mutation category')+
  ggtitle(glue('630 mutations mapped in the tumour, PD63383'))
ggsave(glue('Results/20241106_vaf_tumour_vs_PD63383_630mapped.pdf'))

# how many of those mutations are also present in all normal samples?
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD63383 == 6] # 594
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 6] # 593

twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 6] # 588

twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 6] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 4] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 5] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 3] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 0 & sum_normal_PD63383 == 2] # 0

twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 6 & sum_normal_PD63383 == 0] # 1 chr14
# chr14_105458006_C_A our favourite one!
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 5 & sum_normal_PD63383 == 0] # 0
twins_filtered_vaf[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341 == 4 & sum_normal_PD63383 == 0] # 2 chr3
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
# Accounting for contamination (use purity estimates)

# calculate the VAF in the normal tissue assuming 0.5 VAF in the tumour 
# true VAF (normal sample) = observed VAF (normal sample) - 0 * fraction_normal + 0.5 * fraction_tumour
# fraction_normal = fraction of normal cells, fraction_tumour = fraction of tumour cells 

twins_filtered_vaf_adj = data.table(twins_filtered_vaf)
twins_filtered_vaf_adj = merge(twins_filtered_vaf_adj, twins_tumour_agg[, c('mut_ID', 'tumour_vaf'), with=FALSE], by = 'mut_ID')

for (sample in samples_normal){
  sample_vaf = paste0(sample, '_VAF')
  sample_diff = paste0(sample, '_adj')
  sample_purity = as.numeric(purity_dt[sample==sample_vaf, purity_est])
  twins_filtered_vaf_adj[, (sample_diff) := (get(sample_vaf) - tumour_vaf * (1-sample_purity)) / sample_purity]
}

for (sample in samples_tumour){
  sample_vaf = paste0(sample, '_VAF')
  sample_diff = paste0(sample, '_adj')
  sample_purity = as.numeric(purity_dt[sample==sample_vaf, purity_est])
  twins_filtered_vaf_adj[, (sample_diff) := get(sample_vaf) / sample_purity]
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
# everything in h at lower VAF - how do I know this is contamination or real?
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

# okay that's what we expected that looks fine
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 6] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 5] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 4] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 0 & sum_normal_PD63383_adj == 1] # 0

# sounds like it really cleaned up your stuff then 
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 6 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 5 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 4 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 3 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 2 & sum_normal_PD63383_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_mapped_tumour & sum_normal_PD62341_adj == 1 & sum_normal_PD63383_adj == 1] # 0

######################################################################################################
# Accounting for twin-twin transfusion (spleen samples)





