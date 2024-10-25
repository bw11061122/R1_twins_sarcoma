# Script to analyse the pileup (run 15/10/2024)
# 2024-10-11 - 2024-10-14
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples, high quality (-mq 30 -bq 30)
# 2 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples, low quality (-mq 5 -bq 5)

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
setwd('/Users/bw18/Desktop/1SB')
twins_hq = data.table(read.csv('Data/pileup_merged_20241016.tsv'))
twins_lq = data.table(read.csv('Data/pileup_merged_20241025.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_hq), value = TRUE)
twins_hq[, c(twins_PDv38is) := NULL]

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

# KEY TO SAMPLES PRESENT IN THE DATAFRAME 

# PD62341 NORMAL
# "PD62341n" - liver and nodules
# "PD62341h" - heart (left ventricle)
# "PD62341v" - spleen 
# "PD62341q" - pancreas
# "PD62341aa" - skin 
# "PD62341ad" - cerebellum

# PD62341 TUMOUR 
# "PD62341u" - tumour (liver hilum)
# "PD62341b" - tumour (left adrenal)
# "PD62341ae" - tumour (brain)
# "PD62341ag" - tumour (brain)
# "PD62341aj" - tumour (brain)
# "PD62341ak" - tumour 
# "PD62341am" - tumour 
# "PD62341ap" - tumour 

# PD63383 NORMAL
# "PD63383w" - spleen
# "PD63383t" - liver 
# "PD63383u" - pancreas
# "PD63383ae" - left ventricle
# "PD63383ak" - cerebellum 
# "PD63383bb" - skin 

# PD63383 TUMOUR
# "PD63383ap" - tumour
# "PD63383aq" - tumour

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

######################################################################################################
# CHECKS 

# Check that the dimensions of both pileup dataframes make sense 
paste('Number of columns of the high-quality pileup dataframe:', dim(twins_hq)[2])
paste('Number of columns of the low-quality pileup dataframe:', dim(twins_lq)[2])

# Determine the number of unique mutations identified in at least one sample (ie how many mutations do I have?)
paste('Number of unique mutations identified across samples (high quality pileup):', dim(twins_hq)[1]) # 362,814
paste('Number of unique mutations identified across samples (low quality pileup):', dim(twins_lq)[1]) # 362,814

# Check that class of each column makes sense
sapply(twins_hq, class) # all fine 
sapply(twins_lq, class) # all fine 

######################################################################################################
# COMPARE RESULTS 

samples_mtr = grep("MTR", names(twins_hq), value = TRUE)
samples_dep = grep("DEP", names(twins_hq), value = TRUE)
samples_vaf = grep("VAF", names(twins_hq), value = TRUE)

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
# MAYBE CHECK FILTERS FOR EACH DATAFRAME

twins_mtr_hq = twins_hq[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep_hq = twins_hq[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf_hq = twins_hq[,c('mut_ID', samples_vaf), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf
twins_mtr_hq[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_mtr]
twins_mtr_hq[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_mtr]
twins_mtr_hq[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_mtr]
twins_mtr_hq[,max_mtr := apply(.SD, 1, max), .SDcols = samples_mtr]

twins_dep_hq[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_dep_hq[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_dep_hq[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_dep_hq[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]

twins_vaf_hq[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_vaf]
twins_vaf_hq[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_vaf]
twins_vaf_hq[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_vaf]
twins_vaf_hq[,max_vaf := apply(.SD, 1, max), .SDcols = samples_vaf]

# determine this specifically in normal samples (in case of ploidy changes in the tumour)
twins_mtr_hq[,mean_mtr_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr_hq[,median_mtr_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr_hq[,min_mtr_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr_hq[,max_mtr_normal := apply(.SD, 1, max), .SDcols = samples_normal_mtr]

twins_dep_hq[,mean_dep_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_dep]
twins_dep_hq[,median_dep_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_dep]
twins_dep_hq[,min_dep_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_dep]
twins_dep_hq[,max_dep_normal := apply(.SD, 1, max), .SDcols = samples_normal_dep]

twins_vaf_hq[,mean_vaf_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf_hq[,median_vaf_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf_hq[,min_vaf_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf_hq[,max_vaf_normal := apply(.SD, 1, max), .SDcols = samples_normal_vaf]

# Now to the same for the low quality dataframe 
twins_mtr_lq = twins_lq[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep_lq = twins_lq[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf_lq = twins_lq[,c('mut_ID', samples_vaf), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf
twins_mtr_lq[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_mtr]
twins_mtr_lq[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_mtr]
twins_mtr_lq[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_mtr]
twins_mtr_lq[,max_mtr := apply(.SD, 1, max), .SDcols = samples_mtr]

twins_dep_lq[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_dep_lq[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_dep_lq[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_dep_lq[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]

twins_vaf_lq[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_vaf]
twins_vaf_lq[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_vaf]
twins_vaf_lq[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_vaf]
twins_vaf_lq[,max_vaf := apply(.SD, 1, max), .SDcols = samples_vaf]

# determine this specifically in normal samples (in case of ploidy changes in the tumour)
twins_mtr_lq[,mean_mtr_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr_lq[,median_mtr_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr_lq[,min_mtr_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr_lq[,max_mtr_normal := apply(.SD, 1, max), .SDcols = samples_normal_mtr]

twins_dep_lq[,mean_dep_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_dep]
twins_dep_lq[,median_dep_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_dep]
twins_dep_lq[,min_dep_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_dep]
twins_dep_lq[,max_dep_normal := apply(.SD, 1, max), .SDcols = samples_normal_dep]

twins_vaf_lq[,mean_vaf_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf_lq[,median_vaf_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf_lq[,min_vaf_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf_lq[,max_vaf_normal := apply(.SD, 1, max), .SDcols = samples_normal_vaf]

hist(twins_vaf_lq[,mean_vaf_normal])
hist(twins_vaf_hq[,mean_vaf_normal])

# perhaps we can test which are germline in hq and lq
twins_mtr_hq[,sum_mtr4 := rowSums(.SD>=4), .SDcols = samples_mtr]
twins_mtr_lq[,sum_mtr4 := rowSums(.SD>=4), .SDcols = samples_mtr]



hist(twins_mtr_hq[sum_mtr4!=22,sum_mtr4] %>% unlist(), xlab = 'Number of samples with variant read at MTR=4', main = 'Distribution of # samples with mutation')
hist(twins_mtr_lq[sum_mtr4!=22,sum_mtr4] %>% unlist(), xlab = 'Number of samples with variant read at MTR=4', main = 'Distribution of # samples with mutation')


