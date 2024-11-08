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
col_PD62341 = "#0ac368"
col_PD63383 = "#a249e8"
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
# Inspect mutations that are present in all samples
# Why were these not called as germline?

twins_filtered_vaf_adj[mut_ID %in% muts_normal_all] # 632 of those 

twins_filtered_vaf_adj[, mean_vaf_adj_normal := apply(.SD, 1, mean), .SDcols = samples_normal_adj]
twins_filtered_vaf_adj[, mean_vaf_adj_normal_PD62341 := apply(.SD, 1, mean), .SDcols = samples_normal_PD62341_adj]
twins_filtered_vaf_adj[, mean_vaf_adj_normal_PD63383 := apply(.SD, 1, mean), .SDcols = samples_normal_PD63383_adj]

pdf('Results/20241106_p4_hist_vaf_632muts_in_all_normal_samples.pdf')
hist(twins_filtered_vaf_adj[mut_ID %in% muts_normal_all, mean_vaf_adj_normal] %>% unlist(), 
     xlab = 'Adjusted VAF (agg normal samples)', main = '632 muts: all normal samples', xlim = c(0, 1), breaks=30)
dev.off()

pdf('Results/20241106_p4_hist_vaf_632muts_in_all_normal_samples_PD62341.pdf')
hist(twins_filtered_vaf_adj[mut_ID %in% muts_normal_all, mean_vaf_adj_normal_PD62341] %>% unlist(), 
     xlab = 'Adjusted VAF (agg normal PD62341)', main = '632 muts: PD62341 normal samples', xlim = c(0, 1), breaks=30)
dev.off()

pdf('Results/20241106_p4_hist_vaf_632muts_in_all_normal_samples_PD63383.pdf')
hist(twins_filtered_vaf_adj[mut_ID %in% muts_normal_all, mean_vaf_adj_normal_PD63383] %>% unlist(), 
     xlab = 'Adjusted VAF (agg normal PD63383)', main = '632 muts: PD63383 normal samples', xlim = c(0, 1), breaks=30)
dev.off()

# conclusion: these are clearly not germline based on their VAF 
# are these artefacts or real?

# plot VAF in PD62341 vs PD63383 
# if real / informative, likely enriched just in one twin
ggplot(twins_filtered_vaf_adj[mut_ID %in% muts_normal_all], aes(x = mean_vaf_adj_normal_PD62341, y = mean_vaf_adj_normal_PD63383))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF (agg PD62341 normal)', y = 'Adjusted VAF (agg PD63383 normal)')+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  ggtitle(glue('632 mutations present in all normal samples'))+
  coord_equal(ratio=1)
ggsave(glue('Results/20241106_adj_vaf_normal_PD62341_PD63383.pdf'))

# color by different mutation categories that are potentially interesting 
twins_filtered_vaf_adj[, vaf_ratio := mean_vaf_adj_normal_PD62341 / mean_vaf_adj_normal_PD63383]
pdf('Results/20241106_p4_hist_vaf_632muts_in_all_normal_samples_ratio_PD62341toPD63383.pdf')
hist(twins_filtered_vaf_adj[mut_ID %in% muts_normal_all, vaf_ratio] %>% unlist(), 
     xlab = 'Adjusted VAF ratio PD62341 / PD63383', main = '632 muts: VAF ratio', 
     breaks=30)
dev.off()

twins_filtered_vaf_adj[, vaf_ratio_log := log2(vaf_ratio)]
pdf('Results/20241106_p4_hist_vaf_632muts_in_all_normal_samples_logratio_PD62341toPD63383.pdf')
hist(twins_filtered_vaf_adj[mut_ID %in% muts_normal_all, vaf_ratio_log] %>% unlist(), 
     xlab = 'log2(Adjusted VAF ratio PD62341 / PD63383)', main = '632 muts: log2(VAF ratio)', 
     breaks=30)
abline(v = -0.5, col = 'red')
abline(v = 0.5, col = 'red')
dev.off()

twins_filtered_vaf_adj[, mut_cat := factor(fcase(
  mean_vaf_adj_normal_PD62341 > 0.45 & mean_vaf_adj_normal_PD63383 > 0.45, 'likely germline', 
  vaf_ratio_log < -0.5, 'PD63383 enriched',
  vaf_ratio_log > 0.5, 'PD62341 enriched',
  vaf_ratio_log >= -0.5 & vaf_ratio_log <= 0.5, 'no difference'
))]

ggplot(twins_filtered_vaf_adj[mut_ID %in% muts_normal_all], aes(x = mean_vaf_adj_normal_PD62341, y = mean_vaf_adj_normal_PD63383, col = mut_cat))+
  geom_point(size=1.5, alpha = 1)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF (agg PD62341 normal)', y = 'Adjusted VAF (agg PD63383 normal)', col = 'Mutation category')+
  scale_color_manual(values = c('darkred', 'grey', col_PD62341, col_PD63383))+
  xlim(c(0, 0.6))+
  ylim(c(0, 0.6))+
  ggtitle(glue('632 mutations present in all normal samples'))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241106_adj_vaf_normal_PD62341_PD63383_col_cat.pdf'), width = 6, height = 4)

# adjusted VAF is nice but not exactly great so let's check we get a similar result with aggregate VAF (non-corrected)
purity_dt[, twin := tstrsplit(sample, '_VAF', fixed = TRUE, keep = 1)]
purity_dt[, twin := factor(fcase(twin %in% samples_PD62341, 'PD62341', twin %in% samples_PD63383, 'PD63383'))]
purity_normal_PD62341 = mean(purity_dt[status=='normal' & twin == 'PD62341', purity_est])
purity_normal_PD63383 = mean(purity_dt[status=='normal' & twin == 'PD63383', purity_est])
purity_tumour_PD62341 = mean(purity_dt[status=='tumour' & twin == 'PD62341', purity_est])
purity_tumour_PD63383 = mean(purity_dt[status=='tumour' & twin == 'PD63383', purity_est])
# not sure if I should even adjust this?

twins_agg[, vaf_ratio := PD62341_normal_vaf / PD63383_normal_vaf]
twins_agg[, vaf_ratio_log := log2(vaf_ratio)]
twins_agg[, mut_cat := factor(fcase(
  PD62341_normal_vaf > 0.45 & PD63383_normal_vaf > 0.45, 'likely germline', 
  vaf_ratio_log < -0.5, 'PD63383 enriched',
  vaf_ratio_log > 0.5, 'PD62341 enriched',
  vaf_ratio_log >= -0.5 & vaf_ratio_log <= 0.5, 'no difference'
))]

ggplot(twins_agg[mut_ID %in% muts_normal_all], aes(x = PD62341_normal_vaf, y = PD63383_normal_vaf, col = mut_cat))+
  geom_point(size=1.5, alpha = 1)+
  theme_classic(base_size = 14)+
  labs(x = 'Non-adj VAF (agg PD62341 normal)', y = 'Non-adj VAF (agg PD63383 normal)', col = 'Mutation category')+
  scale_color_manual(values = c('darkred', 'grey', col_PD62341, col_PD63383))+
  xlim(c(0, 0.6))+
  ylim(c(0, 0.6))+
  ggtitle(glue('632 mutations present in all normal samples'))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241106_agg_vaf_normal_PD62341_PD63383_col_cat.pdf'), width = 6, height = 4)

# inspect enriched mutations
twins_filtered_vaf_adj[mut_ID %in% muts_normal_all & mut_cat == 'likely germline', mut_ID]
# "chr10_92674047_G_A" # doesn't look very good 
# "chr12_131428827_C_A" # mapping is funny (checked on blat - not great)
# "chr15_85186362_A_G"  # again, mapping is funny
# "chr17_21783609_C_A" # mapping issues 
# "chr5_127735154_A_G" # mapping issues 
# "chr7_4912517_G_A" # mapping issues 
twins_filtered_vaf_adj[mut_ID %in% muts_normal_all & mut_cat == 'PD62341 enriched', mut_ID]
# "chr15_82483067_G_A" # kind of okay? double check
# "chr16_18539858_T_C" # mapping is not perfect but on blat looks okay
# "chr1_16864604_G_A"  # poor mapping (reads map somewhere else)
# "chr8_123201027_T_C" # does indeed look real
twins_filtered_vaf_adj[mut_ID %in% muts_normal_all & mut_cat == 'PD63383 enriched', mut_ID]
# "chr10_79789115_C_T" # probably okay, but there is some mismapping 
# "chr15_20262122_C_T" # looks okay
# "chr17_21754695_C_T" # phasing is fine but mapping is really not 
# "chr17_36229878_A_G" # rather poor mapping (double check?)
# "chr1_38827952_C_A" # looks okay 
# "chr9_135524032_G_A" # double check but I think too many insertions to be reliable

# plot the distribution of VAFs across samples for the good looking mutations
twins_vaf_adj_sub = twins_filtered_vaf_adj[, c('mut_ID', samples_normal_adj), with=FALSE]
twins_vaf_adj_melt = melt(twins_vaf_adj_sub, id.vars = 'mut_ID')
twins_vaf_adj_melt[, sample := tstrsplit(variable, '_adj', fixed = TRUE, keep = 1)]
twins_vaf_adj_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341n', 'PD62341v', 'PD62341q', 'PD62341h',
  'PD63383ae', 'PD63383ak', 'PD63383u', 'PD63383bb', 'PD63383w', 'PD63383t'
))]
twins_vaf_adj_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
ggplot(twins_vaf_adj_melt[mut_ID=='chr8_123201027_T_C'], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF (adjusted)', col = 'Twin')+
  ggtitle(glue('Mutation: chr8:123201027, T>C'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241106_dist_samples_mut_chr8_PD62341_enriched.pdf'))

ggplot(twins_vaf_adj_melt[mut_ID=='chr1_38827952_C_A'], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF (adjusted)', col = 'Twin')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241106_dist_samples_mut_chr1_PD63383_enriched.pdf'))

ggplot(twins_vaf_adj_melt[mut_ID=='chr15_20262122_C_T'], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF (adjusted)', col = 'Twin')+
  ggtitle(glue('Mutation: chr15:20262122, C>T'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241106_dist_samples_mut_chr15_PD63383_enriched.pdf'))

twins_agg[mut_ID %in% muts_normal_all & mut_cat == 'likely germline', mut_ID]
# "chr12_131428827_C_A" # mapping issues 
# "chr15_85186362_A_G" # mapping issues 
# "chr17_21783609_C_A" # mapping issues
# "chr5_127735154_A_G" # mapping issues 
# "chr7_4912517_G_A" # mapping issues
twins_agg[mut_ID %in% muts_normal_all & mut_cat == 'PD62341 enriched', mut_ID]
# "chr15_82483067_G_A" # poor mapping 
# "chr1_16864604_G_A" # poor mapping  
# "chr8_123201027_T_C" # looks okay
twins_agg[mut_ID %in% muts_normal_all & mut_cat == 'PD63383 enriched', mut_ID]
# "chr1_38827952_C_A" # looks okay

# Are all these mutations artefacts?
# Found mutations on chr12 that all look great - what is going on?
mut_chr12_cluster = c("chr12_31851749_C_T",  "chr12_31851985_G_A",  "chr12_31895265_G_A",
                      "chr12_31852349_C_G",  "chr12_31852949_G_A",  "chr12_31856122_C_T", 
                      "chr12_31867266_G_A",  "chr12_31869715_A_G",  "chr12_31872355_C_T", 
                      "chr12_31879431_C_T",  "chr12_31880061_G_A",  "chr12_31883693_C_T", 
                      "chr12_31885515_G_A",  "chr12_31885689_A_G",  "chr12_31887208_C_T" )

ggplot(twins_vaf_adj_melt[mut_ID %in% mut_chr12_cluster], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF (adjusted)', col = 'Twin')+
  ggtitle(glue('Mutation on chr12 31,850,000-31,900,000'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.6))
ggsave(glue('Results/20241106_dist_samples_mut_chr12.pdf'))

pdf('Results/20241106_p5_hist_chr12_PD62341_PD63383.pdf')
par(mfrow = c(2, 1))
hist(twins_agg[mut_ID %in% mut_chr12_cluster, PD63383_normal_vaf], xlab = 'VAF', xlim = c(0, 1), main = 'PD63383 chr12', breaks = 6)
hist(twins_agg[mut_ID %in% mut_chr12_cluster, PD62341_normal_vaf], xlab = 'VAF', xlim = c(0, 1), main = 'PD62341 chr12', breaks = 6)
dev.off()

# plot the coverage vs VAF for these mutations
twins_agg[, mut_chr12 := factor(fcase(
  mut_ID %in% mut_chr12_cluster, 'segment chr12',
  !mut_ID %in% mut_chr12_cluster, 'other'))]

ggplot(twins_agg[mut_ID %in% muts_normal_all] %>% arrange(mut_chr12), aes(x = sum_PD62341_dep, y = PD62341_normal_vaf, col = mut_chr12))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred'))+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD62341 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_632_PD62341.pdf'), height=5, width=7)

ggplot(twins_agg[mut_ID %in% muts_normal_all] %>% arrange(mut_chr12), aes(x = sum_PD63383_dep, y = PD63383_normal_vaf, col = mut_chr12))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred'))+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD62341 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_632_PD63383.pdf'), height=5, width=7)

# similar thing on chr 14?
mut_chr14_cluster = c("chr14_19725920_G_T",  "chr14_19735467_A_C",  "chr14_19770060_C_T",  "chr14_19778921_G_T", 
                      "chr14_19803896_T_A",  "chr14_19814466_C_T",  "chr14_19823326_A_T", "chr14_19880026_C_G", 
                      "chr14_19881983_C_T",  "chr14_19887911_C_T",  "chr14_19895402_G_A",  "chr14_19899009_A_G", 
                      "chr14_19899033_A_G",  "chr14_19903456_A_G")


twins_agg[, mut_cluster := factor(fcase(
  mut_ID %in% mut_chr12_cluster, 'segment chr12',
  mut_ID %in% mut_chr14_cluster, 'segment chr14',
  !mut_ID %in% c(mut_chr12_cluster, mut_chr14_cluster), 'other'))]

ggplot(twins_agg[mut_ID %in% muts_normal_all] %>% arrange(mut_cluster), aes(x = sum_PD62341_dep, y = PD62341_normal_vaf, col = mut_cluster))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred', '#0db793'))+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD62341 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_632_PD62341_clusters.pdf'), height=5, width=7)

ggplot(twins_agg[mut_ID %in% muts_normal_all] %>% arrange(mut_cluster), aes(x = sum_PD63383_dep, y = PD63383_normal_vaf, col = mut_cluster))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred', '#0db793'))+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD63383 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_632_PD63383_clusters.pdf'), height=5, width=7)

# can I plot the distribution of these mutations for each chromosome
twins_agg[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
twins_agg[, pos := tstrsplit(mut_ID, '_', fixed=TRUE, keep=2)]
twins_agg[, pos := as.numeric(pos)]

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

twins_agg = merge(twins_agg, chr_lengths_dt, by = 'Chrom')
twins_agg[, Chrom := factor(Chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                                "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                                "chr20", "chr21", "chr22", "chrX"))]

for (chr in Chrom){
  dt = twins_agg[Chrom == chr & mut_ID %in% muts_normal_all]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = PD63383_normal_vaf))+
    geom_point(size=2.5)+
    theme_classic(base_size = 12)+
    scale_color_manual(values = '#7b13d4')+
    labs(x = 'Genomic position', y = 'VAF (PD63383 normal)')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241106_dist_across_genome_632_{chr}_PD63383.pdf'), height=3, width=4.5)
}

# compare to coverage across all mutations in the set
twins_agg[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  mut_ID %in% muts_normal_only, 'only normal samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_all, muts_normal_only, muts_tumour_all_only, muts_tumour_only), 'other'))]

ggplot(twins_agg, aes(x = sum_PD62341_dep, y = PD62341_normal_vaf, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Mutation class')+
  ggtitle(glue('PD62341 normal samples (1002 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_1002_PD62341.pdf'), height=5, width=7)

ggplot(twins_agg, aes(x = sum_PD63383_dep, y = PD63383_normal_vaf, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Mutation class')+
  ggtitle(glue('PD63383 normal samples (1002 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241106_cov_vs_vaf_mut_class_1002_PD63383.pdf'), height=5, width=7)

twins_dt[, sum_normal_PD63383_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_dt[, sum_normal_PD63383_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_dt[, sum_normal_PD62341_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_dt[, sum_normal_PD62341_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_dt[, sum_normal_PD63383_vaf := sum_normal_PD63383_mtr / sum_normal_PD63383_dep]
twins_dt[, sum_normal_PD62341_vaf := sum_normal_PD62341_mtr / sum_normal_PD62341_dep]

par(mfrow = c(1,1))

pdf('Results/20241106_p6_hist_PD62341_allmuts.pdf')
hist(twins_dt[sum_normal_dep < 600, sum_normal_PD62341_dep], 
     xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, all mutations (362,814)', breaks = 50)
abline(v = median(twins_dt[sum_normal_dep < 600, sum_normal_PD62341_dep]), col = 'red', lwd = 2)
dev.off()

pdf('Results/20241106_p6_hist_PD63383_allmuts.pdf')
hist(twins_dt[sum_normal_dep < 600, sum_normal_PD63383_dep], 
     xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, all mutations (362,814)', breaks = 50)
abline(v = median(twins_dt[sum_normal_dep < 600, sum_normal_PD63383_dep]), col = 'red', lwd = 2)
dev.off()

pdf('Results/20241106_p6_hist_PD62341_632muts.pdf')
hist(twins_agg[mut_ID %in% muts_normal_all,sum_PD62341_dep] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, 632 normal')
abline(v = median(twins_agg[mut_ID %in% muts_normal_all, sum_PD62341_dep]), col = 'blue', lwd = 2) # 264.5
dev.off()

pdf('Results/20241106_p6_hist_PD63383_632muts.pdf')
hist(twins_agg[mut_ID %in% muts_normal_all,sum_PD63383_dep] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, 632 normal')
abline(v = median(twins_agg[mut_ID %in% muts_normal_all, sum_PD63383_dep]), col = 'blue', lwd = 2) # 264.5
dev.off()

pdf('Results/20241106_p6_hist_PD62341_1002muts.pdf')
hist(twins_agg[,sum_PD62341_dep] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, 1002 normal')
abline(v = median(twins_agg[,sum_PD62341_dep]), col = 'purple', lwd = 2) # 264.5
dev.off()

pdf('Results/20241106_p6_hist_PD63383_1002muts.pdf')
hist(twins_agg[,sum_PD63383_dep] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, 1002 normal')
abline(v = median(twins_agg[,sum_PD63383_dep]), col = 'purple', lwd = 2) # 264.5
dev.off()

######################################################################################################
# how do mutations present in all normal samples look like in tumour samples?

######################################################################################################
# plot distribution of mutations present in all tumour samples across the genome 
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
  ggsave(glue('Results/20241106_dist_across_genome_641_tumour_{chr}.pdf'), height=3, width=4.5)
}

twins_tumour_agg[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_all, 'all samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  mut_ID %in% muts_normal_only, 'only normal samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_all, muts_normal_only, muts_tumour_all_only, muts_tumour_only), 'other'))]

for (chr in Chrom){
  dt = twins_agg[Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = tumour_vaf, col = mut_class))+
    geom_point(size=2.5)+
    theme_classic(base_size = 10)+
    labs(x = 'Genomic position', y = 'VAF (agg tumour)', col = 'Mutation category')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241106_dist_across_genome_1002_tumour_{chr}_cat.pdf'), height=3, width=5)
}


######################################################################################################
# Accounting for twin-twin transfusion (spleen samples) - I want to have quantitative estimates of how much transfer there is 





