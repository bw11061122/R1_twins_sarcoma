# Script to analyse the pileup (run 15/10/2024)
# 2024-10-21
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# with filters 
# NOTE that for the moment, I am investigating if there is sth wrong in how the pileup was run
# I may alter the input df once everything is resolved

# NOTE ON ADDING COLUMNS WITH FILTERS 
# I want to add a column which says 0 if sample is fine and 1 if sample is wrong 

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

###################################################################################################################################
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB/Analysis')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241018.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 2,559
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

###################################################################################################################################
# Specify settings for plotting 
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
samples_names = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb",
                  "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq",
                  "PD62341b", "PD62341h", "PD62341n", "PD62341u")
samples_normal = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341h", "PD62341n")
samples_tumour = c("PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq", "PD62341b", "PD62341u")
samples_PD62341 = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
samples_PD63383 = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383ap", "PD63383aq", "PD63383bb")
samples_normal_PD62341 = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD62341h", "PD62341n")
samples_normal_PD63383 = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb")
samples_tumour_PD62341 = c("PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap", "PD62341b", "PD62341u")
samples_tumour_PD63383 = c( "PD63383ap", "PD63383aq")

######################################################################################################
# Add column to indicate chromosomes lost in the tumour
twins_dt[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'
))]

######################################################################################################
# Is there contamination of PD62341h sample with tumour cells?

# Names of samples with correct descriptor
samples_normal_mtr = paste(samples_normal, 'MTR', sep='_')
samples_tumour_mtr = paste(samples_tumour, 'MTR', sep='_')
samples_PD62341_mtr = paste(samples_PD62341, 'MTR', sep='_')
samples_PD63383_mtr = paste(samples_PD63383, 'MTR', sep='_')
samples_normal_PD62341_mtr = paste(samples_normal_PD62341, 'MTR', sep='_')
samples_normal_PD63383_mtr = paste(samples_normal_PD63383, 'MTR', sep='_')
samples_tumour_PD62341_mtr = paste(samples_tumour_PD62341, 'MTR', sep='_')
samples_tumour_PD63383_mtr = paste(samples_tumour_PD63383, 'MTR', sep='_')

samples_normal_vaf = paste(samples_normal, 'VAF', sep='_')
samples_tumour_vaf = paste(samples_tumour, 'VAF', sep='_')
samples_PD62341_vaf = paste(samples_PD62341, 'VAF', sep='_')
samples_PD63383_vaf = paste(samples_PD63383, 'VAF', sep='_')
samples_normal_PD62341_vaf = paste(samples_normal_PD62341, 'VAF', sep='_')
samples_normal_PD63383_vaf = paste(samples_normal_PD63383, 'VAF', sep='_')
samples_tumour_PD62341_vaf = paste(samples_tumour_PD62341, 'VAF', sep='_')
samples_tumour_PD63383_vaf = paste(samples_tumour_PD63383, 'VAF', sep='_')

samples_normal_dep = paste(samples_normal, 'DEP', sep='_')
samples_tumour_dep = paste(samples_tumour, 'DEP', sep='_')
samples_PD62341_dep = paste(samples_PD62341, 'DEP', sep='_')
samples_PD63383_dep = paste(samples_PD63383, 'DEP', sep='_')
samples_normal_PD62341_dep = paste(samples_normal_PD62341, 'DEP', sep='_')
samples_normal_PD63383_dep = paste(samples_normal_PD63383, 'DEP', sep='_')
samples_tumour_PD62341_dep = paste(samples_tumour_PD62341, 'DEP', sep='_')
samples_tumour_PD63383_dep = paste(samples_tumour_PD63383, 'DEP', sep='_')

# For each mutation, create subsetted dataframes VAF, DEP, MTR
cols_mtr = grep("MTR", names(twins_filtered_dt), value = TRUE)
cols_dep = grep("DEP", names(twins_filtered_dt), value = TRUE)
cols_vaf = grep("VAF", names(twins_filtered_dt), value = TRUE)

cols_normal_mtr = paste(samples_normal, 'MTR', sep='_')
cols_normal_dep = paste(samples_normal, 'DEP', sep='_')
cols_normal_vaf = paste(samples_normal, 'VAF', sep='_')

# select VAF values of PD62341h and all tumour samples
# tumour sample from the liver: PD62341u (compare to this cf other tumour samples)
twins_dt_liver_mtr = twins_dt[, c('mut_ID', 'PD62341h_MTR', samples_tumour_mtr), with=FALSE]
twins_dt_liver_vaf = twins_dt[, c('mut_ID', 'PD62341h_VAF', samples_tumour_vaf), with=FALSE]

# create a sample, we cannot plot all mutations
twins_dt_liver_mtr_sample = data.table(sample_n(twins_dt_liver_mtr, 1000, replace = FALSE))
twins_dt_liver_vaf_sample = data.table(sample_n(twins_dt_liver_vaf, 1000, replace = FALSE))

# add mutations which are shared
twins_dt_liver_mtr_hq = twins_dt_liver_mtr[mut_ID %in% muts_hq]
twins_dt_liver_vaf_hq = twins_dt_liver_vaf[mut_ID %in% muts_hq] 
  
twins_dt_liver_mtr1 = rbind(twins_dt_liver_mtr_sample, twins_dt_liver_mtr_hq)
twins_dt_liver_vaf1 = rbind(twins_dt_liver_vaf_sample, twins_dt_liver_vaf_hq)

twins_dt_liver_mtr_melt = melt(twins_dt_liver_mtr1, id.vars = 'mut_ID')
twins_dt_liver_vaf_melt = melt(twins_dt_liver_vaf1, id.vars = 'mut_ID')

twins_dt_liver_mtr_melt[, sample := factor(fcase(
  variable == 'PD62341h_MTR', 'normal: liver',
  variable == 'PD62341u_MTR', 'tumour: liver',
  !variable %in% c('PD62341h_MTR', 'PD62341u_MTR'), 'tumour: other'
))]
twins_dt_liver_vaf_melt[, sample := factor(fcase(
  variable == 'PD62341h_VAF', 'normal: liver',
  variable == 'PD62341u_VAF', 'tumour: liver',
  !variable %in% c('PD62341h_VAF', 'PD62341u_VAF'), 'tumour: other'
))]

ggplot(twins_dt_liver_mtr_melt, aes(x = value, y = mut_ID, col = sample))+
  geom_point(size=2, position=position_dodge(.75))+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour, '#decfa4'))+
  labs(x = 'MTR', y = 'Mutation')+
  ggtitle('comparison of PD62341h to tumour and normal samples')+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 

ggplot(twins_dt_liver_vaf_melt, aes(x = value, y = mut_ID, col = sample))+
  geom_point(size=2, position=position_dodge(.75))+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour, '#decfa4'))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle('comparison of PD62341h to tumour and normal samples')+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 

# only the 13 shared mutations + paired 13 random 

# sample 13 mutations at random
twins_dt_liver_mtr_sample2 = data.table(sample_n(twins_dt_liver_mtr, 13, replace = FALSE))
twins_dt_liver_vaf_sample2 = data.table(sample_n(twins_dt_liver_vaf, 13, replace = FALSE))

twins_dt_liver_mtr_sample2[, mut := factor('random')]
twins_dt_liver_mtr_hq[, mut := factor('shared')]

twins_dt_liver_vaf_sample2[, mut := factor('random')]
twins_dt_liver_vaf_hq[, mut := factor('shared')]

twins_dt_liver_mtr2 = rbind(twins_dt_liver_mtr_sample2, twins_dt_liver_mtr_hq)
twins_dt_liver_vaf2 = rbind(twins_dt_liver_vaf_sample2, twins_dt_liver_vaf_hq)

twins_dt_liver_mtr_melt2 = melt(twins_dt_liver_mtr2, id.vars = c('mut_ID', 'mut'))
twins_dt_liver_vaf_melt2 = melt(twins_dt_liver_vaf2, id.vars = c('mut_ID', 'mut'))

twins_dt_liver_mtr_melt2[, sample := factor(fcase(
  variable == 'PD62341h_MTR', 'normal: liver',
  variable == 'PD62341u_MTR', 'tumour: liver',
  !variable %in% c('PD62341h_MTR', 'PD62341u_MTR'), 'tumour: other'))]
twins_dt_liver_vaf_melt2[, sample := factor(fcase(
  variable == 'PD62341h_VAF', 'normal: liver',
  variable == 'PD62341u_VAF', 'tumour: liver',
  !variable %in% c('PD62341h_VAF', 'PD62341u_VAF'), 'tumour: other'))]

ggplot(twins_dt_liver_mtr_melt2, aes(x = value, y = mut_ID, col = sample))+
  geom_point(size=1.5)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour, '#decfa4'))+
  labs(x = 'MTR', y = 'Mutation')+
  ggtitle('comparison of PD62341h to tumour and normal samples')+
  facet_grid(mut ~ .)+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 

ggplot(twins_dt_liver_vaf_melt2, aes(x = value, y = mut_ID, col = sample))+
  geom_point(size=2)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour, '#decfa4'))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle('comparison of PD62341h to tumour and normal samples')+
  theme(axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 

######################################################################################################
# Is there a way I can test contamination?

twins_dt_liver_mtr_melt_hq = melt(twins_dt_liver_mtr_hq, id.vars = c('mut_ID', 'mut'))
twins_dt_liver_vaf_melt_hq = melt(twins_dt_liver_vaf_hq, id.vars = c('mut_ID', 'mut'))

twins_dt_liver_mtr_melt_hq[, sample := factor(fcase(
  variable == 'PD62341h_MTR', 'normal: liver',
  variable == 'PD62341u_MTR', 'tumour: liver',
  !variable %in% c('PD62341h_MTR', 'PD62341u_MTR'), 'tumour: other'))]
twins_dt_liver_vaf_melt_hq[, sample := factor(fcase(
  variable == 'PD62341h_VAF', 'normal: liver',
  variable == 'PD62341u_VAF', 'tumour: liver',
  !variable %in% c('PD62341h_VAF', 'PD62341u_VAF'), 'tumour: other'))]

ggplot(twins_dt_liver_mtr_melt_hq, aes(x = value, y = mut_ID, col = sample))+
  geom_point(size=1.5)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour, '#decfa4'))+
  labs(x = 'MTR', y = 'Mutation')+
  ggtitle('comparison of PD62341h to tumour and normal samples')

ggplot(twins_dt_liver_vaf_melt_hq, aes(x = value, y = mut_ID, col = sample))+
  geom_point(size=2)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal, col_tumour, '#decfa4'))+
  labs(x = 'VAF', y = 'Mutation')+
  ggtitle('comparison of PD62341h to tumour and normal samples')

######################################################################################################
######################################################################################################

# could I do something like clustering and see if the sample is closer to PD62341u?
# tried to do this but can't figure out how to do this in a way that is actually useful 
library(tibble)
library(ggrepel)
library("FactoMineR") ## PCA
library("factoextra") ## PCA 

twins_dt_liver_mtr_sample3 = data.table(sample_n(twins_dt_liver_mtr, 100, replace = FALSE))
twins_dt_liver_vaf_sample3 = data.table(sample_n(twins_dt_liver_vaf, 100, replace = FALSE))

twins_dt_liver_mtr3 = rbind(twins_dt_liver_mtr_sample3, twins_dt_liver_mtr_hq[,1:12])
twins_dt_liver_vaf3 = rbind(twins_dt_liver_vaf_sample3, twins_dt_liver_vaf_hq[,1:12])

twins_dt_liver_mtr_melt3 = melt(twins_dt_liver_mtr3, id.vars = 'mut_ID')
twins_dt_liver_vaf_melt3 = melt(twins_dt_liver_vaf3, id.vars = 'mut_ID')

twins_dt_liver_mtr_melt3[, sample := factor(fcase(
  variable == 'PD62341h_MTR', 'normal: liver',
  variable == 'PD62341u_MTR', 'tumour: liver',
  !variable %in% c('PD62341h_MTR', 'PD62341u_MTR'), 'tumour: other'))]
twins_dt_liver_vaf_melt3[, sample := factor(fcase(
  variable == 'PD62341h_VAF', 'normal: liver',
  variable == 'PD62341u_VAF', 'tumour: liver',
  !variable %in% c('PD62341h_VAF', 'PD62341u_VAF'), 'tumour: other'))]

twins_dt_liver_mtr3_t = data.table(t(twins_dt_liver_mtr3))
twins_dt_liver_vaf3_t = data.table(t(twins_dt_liver_vaf3))

twins_dt_liver_vaf3_t = twins_dt_liver_vaf3_t[2:12,]
twins_dt_liver_vaf3_t = data.table(sapply(twins_dt_liver_vaf3_t, as.numeric))
twins_dt_liver_vaf3_t[, sample := factor(c('normal liver', 
                                           'tumour liver', 
                                           rep('tumour other', 9)))]

res.pca = PCA(twins_dt_liver_vaf3_t[,1:113],  graph = FALSE) # PCA on JUST numerical
pca <- fviz_pca_ind(res.pca,
                             label = "none", # hide individual labels
                             habillage = twins_dt_liver_vaf3_t[,'sample'] %>% unlist(),
                             title="PCA")

data(iris)
head(iris)
res.pca <- prcomp(iris[, -5],  scale = TRUE)
fviz_pca_ind(res.pca, label="none", habillage=iris$Species)
