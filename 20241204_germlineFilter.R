###################################################################################################################################
# SCRIPT 2

# Script to remove germline mutations in copy number-altered regions
# 2024-11-26
# Barbara Walkowiak bw18

# INPUT: 
# 1 PURPE dt (.cnv.segment) for each sample 

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
twins_dt = data.table(read.csv('Data/twins_dt_20241114_1134.csv')) # import high quality pileup with added filters

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Create a dataframe that only includes mutations retained post filtering  
muts_included = read.table('Data/mutations_include_20241114_1134.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts_included)) # 1134

muts_germline = read.table('Data/mutations_putativeGermline_20241114.txt') %>% unlist()
paste('Number of mutations identified as putative germline:', length(muts_germline)) # 333,031

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

col_background = '#c8c8c8'
col_germline = '#d3b28b'
col_other = '#25b6db'
col_normal_all = '#5265c9'
col_normal_all_removed = '#b9dff5'
col_clusters = '#bf068a'
col_coverage = '#bf7a06'
col_clusters_coverage = '#810cb9'
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
twins_dt[, dep_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_dt[, dep_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_dt[, dep_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_dt[, dep_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_dt[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]
twins_dt[, dep_all_tumour := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_dt[, dep_all := rowSums(.SD), .SDcols = samples_dep]

twins_dt[, mtr_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_dt[, mtr_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_dt[, mtr_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_dt[, mtr_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
twins_dt[, mtr_all_normal := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_dt[, mtr_all_tumour := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_dt[, mtr_all := rowSums(.SD), .SDcols = samples_mtr]

twins_dt[, vaf_all_normal_PD62341 := mtr_all_normal_PD62341 / dep_all_normal_PD62341]
twins_dt[, vaf_all_normal_PD63383 := mtr_all_normal_PD63383 / dep_all_normal_PD63383]
twins_dt[, vaf_all_tumour_PD62341 := mtr_all_tumour_PD62341 / dep_all_tumour_PD62341]
twins_dt[, vaf_all_tumour_PD63383 := mtr_all_tumour_PD63383 / dep_all_tumour_PD63383]
twins_dt[, vaf_all_normal := mtr_all_normal / dep_all_normal]
twins_dt[, vaf_all_tumour := mtr_all_tumour / dep_all_tumour]
twins_dt[, vaf_all := mtr_all / dep_all]
