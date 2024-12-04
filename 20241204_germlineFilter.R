###################################################################################################################################
# SCRIPT 2b

# Script to show filtering of mutations thought to be on germline copy number-altered regions
# 2024-11-26
# Barbara Walkowiak bw18

# INPUT: 
# 1 pileup
# 2 existing filters 
# 3 List of mutations from previous filters 

# OUTPUT:
# 1 list of mutations to exclude based on presence in copy number germline variants 

# Script to identify regions of copy number alterations (deletions or duplications) in the germline
# Identify mutations in those regions where the VAF is 0.3 and so those can be explained as a germline copy number alteration 

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
twins_dt = fread('Data/twins_dt_20241204_1069.csv') # import high quality pileup with added filters

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Identify and exclude mutations that failed QC   
muts_failqc = twins_dt[(sum_req_filters >= 1 & f7_likelyGermline_bothTwins==0) | (sum_req_filters > 1 & f7_likelyGermline_bothTwins==1), mut_ID] %>% unlist()  
paste('Number of mutations that failed QC:', length(muts_failqc)) # 28,771 # so 8% fails QC, not great but not disastrous either?

# Informative sets of mutations 
# Mutations included int the final set of 1,069
muts_included = read.table('Data/mutations_include_20241204_1069.txt') %>% unlist()

# Mutations which were retained as non-germline out of the 1,069
muts_599 = read.table('Data/mutations_include_20241114_599.txt') %>% unlist()
muts_germline_cn = setdiff(muts_included, muts_599)

# Mutations which passed manual QC 
muts_dt = data.table(read.csv('Data/20241114_599muts_QCJbrowse.csv', header = T))
muts_255 = muts_dt[Jbrowse.quality == 'Y', mut_ID] %>% unlist()
paste('Number of mutations that passed QC:', length(muts_255)) # 255 
muts_failed_qc = setdiff(muts_599, muts_255)

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

col_germline = '#d3b28b'
col_excluded = '#a80505'

######################################################################################################
# SAMPLES

# Create lists of possible samples of interest
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
# AGGREGATED MTR, DEP and VAF (by sample category) 

# Aggregated VAF in tumour samples 
twins_dt_filtered[, dep_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_dt_filtered[, dep_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_dt_filtered[, dep_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_dt_filtered[, dep_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_dt_filtered[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]
twins_dt_filtered[, dep_all_tumour := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_dt_filtered[, dep_all := rowSums(.SD), .SDcols = samples_dep]

twins_dt_filtered[, mtr_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_dt_filtered[, mtr_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_dt_filtered[, mtr_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_dt_filtered[, mtr_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
twins_dt_filtered[, mtr_all_normal := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_dt_filtered[, mtr_all_tumour := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_dt_filtered[, mtr_all := rowSums(.SD), .SDcols = samples_mtr]

twins_dt_filtered[, vaf_all_normal_PD62341 := mtr_all_normal_PD62341 / dep_all_normal_PD62341]
twins_dt_filtered[, vaf_all_normal_PD63383 := mtr_all_normal_PD63383 / dep_all_normal_PD63383]
twins_dt_filtered[, vaf_all_tumour_PD62341 := mtr_all_tumour_PD62341 / dep_all_tumour_PD62341]
twins_dt_filtered[, vaf_all_tumour_PD63383 := mtr_all_tumour_PD63383 / dep_all_tumour_PD63383]
twins_dt_filtered[, vaf_all_normal := mtr_all_normal / dep_all_normal]
twins_dt_filtered[, vaf_all_tumour := mtr_all_tumour / dep_all_tumour]
twins_dt_filtered[, vaf_all := mtr_all / dep_all]

######################################################################################################
# Calculate distance to the next mutation 

# first, get chromosome and position columns
twins_dt_filtered[, Chrom := tstrsplit(mut_ID, '_', fixed=T, keep=1)]
twins_dt_filtered[, pos := tstrsplit(mut_ID, '_', fixed=T, keep=2)]
twins_dt_filtered[, pos := as.numeric(pos)]

# only for included mutations (1069)
twins_dt_included = twins_dt_filtered[mut_ID %in% muts_included]
setkey(twins_dt_included, Chrom, pos)
twins_dt_included[ , diff_up := pos - shift(pos), by = Chrom] # difference to the mutation before this one
twins_dt_included[ , diff_down := abs(pos - shift(pos, type = 'lead')), by = Chrom] # difference to the mutation after this one
twins_dt_included[, min_diff := pmin(diff_up, diff_down, na.rm=TRUE)] # nearest next mutation

######################################################################################################
# For each mutation, the aim is to plot coverage against distance to the next mutation 

twins_dt_included[, mut_cat := factor(fcase(
  mut_ID %in% muts_germline_cn, 'germline CN',
  !mut_ID %in% muts_germline_cn & mut_ID %in% muts_255, 'not germline, passed QC',
  !mut_ID %in% muts_germline_cn & !mut_ID %in% muts_255, 'not germline, failed QC'
))]

ggplot(twins_dt_included, aes(x = dep_all_normal, y = log10(min_diff), colour = mut_cat))+ # diff 0 if first mutation
  geom_point(alpha= 0.8)+
  theme_classic(base_size = 13)+
  labs(x = 'coverage (Aggregated normal samples)', y = 'log10(distance to nearest mutation)', colour = 'Mutation category')+
  ggtitle('1069 included mutations')+
  guides(color = guide_legend(override.aes = list(alpha = 1) ) )+
  scale_color_manual(values = c(col_germline, col_excluded, col_normal))
ggsave('Results/20241204_filteringGermline_covVslogDistance.pdf', height = 4, width = 4)
