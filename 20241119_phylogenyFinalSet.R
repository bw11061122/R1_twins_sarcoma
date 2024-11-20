###################################################################################################################################
# SCRIPT 4

# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-11-19
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# OUTPUT:
# 1 lists of mutations in specific categories of interest for the phylogeny
# 2 plots for each category of mutations of interest
# 3 Final reconstructed phylogeny 

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
library(RColorBrewer)
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

# Load QC-validated mutations (final set of mutations to be used)
muts_dt = data.table(read.csv('Data/20241114_599muts_QCJbrowse.csv', header = T))
muts = muts_dt[Jbrowse.quality == 'Y', mut_ID] %>% unlist()
paste('Number of mutations that passed QC:', length(muts)) # 258 

# Create a dataframe with mutations of interested retained
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Import dataframe with purity estimates
purity_dt = data.table(read.csv('Data/20241114_estimates_tumour_cont_27muts_median.csv'))

# Import list of driver genes (from Henry Lee-Six, 12/11/2024)
driver_genes_dt = data.table(read.csv('Data/HLS_fibromatoses_driver_list_with_fusions.csv', header=T))
driver_genes = driver_genes_dt[, gene] %>% unlist()

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
col_PD62341_spleen = "#148259"
col_PD63383_spleen = "#53088e"
col_PD63383_skin = '#d0bbe1'
col_bar = '#e87811'

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
twins_filtered_dt[, sum_all_mtr := rowSums(.SD>=4), .SDcols = c(samples_tumour_mtr, samples_normal_mtr)]
twins_filtered_dt[, sum_tumour_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, sum_normal_mtr := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_dt[, sum_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_dt[, sum_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_dt[, sum_tumour_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, sum_tumour_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, sum_normal_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, sum_normal_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# Add presence / absence based on VAF data
twins_filtered_dt[, sum_all_vaf := rowSums(.SD>=0.1), .SDcols = c(samples_tumour_vaf, samples_normal_vaf)]
twins_filtered_dt[, sum_tumour_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_filtered_dt[, sum_normal_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_vaf]
twins_filtered_dt[, sum_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD62341_vaf]
twins_filtered_dt[, sum_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD63383_vaf]
twins_filtered_dt[, sum_tumour_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_filtered_dt[, sum_tumour_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]
twins_filtered_dt[, sum_normal_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_filtered_dt[, sum_normal_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

twins_filtered_dt[, sum_all_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_all_mtr', 'sum_all_vaf')]
twins_filtered_dt[, sum_tumour_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_mtr', 'sum_tumour_vaf')]
twins_filtered_dt[, sum_normal_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_mtr', 'sum_normal_vaf')]
twins_filtered_dt[, sum_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_PD62341_mtr', 'sum_PD62341_vaf')]
twins_filtered_dt[, sum_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_PD63383_mtr', 'sum_PD63383_vaf')]
twins_filtered_dt[, sum_tumour_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_PD62341_mtr', 'sum_tumour_PD62341_vaf')]
twins_filtered_dt[, sum_tumour_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_PD63383_mtr', 'sum_tumour_PD63383_vaf')]
twins_filtered_dt[, sum_normal_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD62341_mtr', 'sum_normal_PD62341_vaf')]
twins_filtered_dt[, sum_normal_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD63383_mtr', 'sum_normal_PD63383_vaf')]

# Add a section for clean normal PD63383 samples (from script 2, evident that skin (PD63383bb) is contaminated)
samples_normal_PD63383_clean = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak")
samples_normal_PD63383_clean_dep = paste0(samples_normal_PD63383_clean, '_DEP')
samples_normal_PD63383_clean_mtr = paste0(samples_normal_PD63383_clean, '_MTR')
samples_normal_PD63383_clean_vaf = paste0(samples_normal_PD63383_clean, '_VAF')

twins_filtered_dt[, sum_normal_PD63383_clean_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_clean_mtr]
twins_filtered_dt[, sum_normal_PD63383_clean_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_clean_vaf]
twins_filtered_dt[, sum_normal_PD63383_clean_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD63383_clean_mtr', 'sum_normal_PD63383_clean_vaf')]

twins_filtered_dt[, agg_normal_PD63383_clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_clean_mtr]
twins_filtered_dt[, agg_normal_PD63383_clean_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_clean_dep]
twins_filtered_dt[, agg_normal_PD63383_clean_vaf := agg_normal_PD63383_clean_mtr / agg_normal_PD63383_clean_dep]

# Add a section for clean normal PD62341 samples (exclude spleen which we know got transfusion)
samples_normal_PD62341_clean = c("PD62341n", "PD62341h", "PD62341ad", "PD62341aa", "PD62341q")
samples_normal_PD62341_clean_dep = paste0(samples_normal_PD62341_clean, '_DEP')
samples_normal_PD62341_clean_mtr = paste0(samples_normal_PD62341_clean, '_MTR')
samples_normal_PD62341_clean_vaf = paste0(samples_normal_PD62341_clean, '_VAF')

twins_filtered_dt[, sum_normal_PD62341_clean_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_clean_mtr]
twins_filtered_dt[, sum_normal_PD62341_clean_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_clean_vaf]
twins_filtered_dt[, sum_normal_PD62341_clean_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD62341_clean_mtr', 'sum_normal_PD62341_clean_vaf')]

twins_filtered_dt[, agg_normal_PD62341_clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_clean_mtr]
twins_filtered_dt[, agg_normal_PD62341_clean_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_clean_dep]
twins_filtered_dt[, agg_normal_PD62341_clean_vaf := agg_normal_PD62341_clean_mtr / agg_normal_PD62341_clean_dep]

# Add a section to identify contaminated normal samples 
samples_normal_contaminated = c("PD62341aa", "PD62341h", "PD63383bb")
samples_normal_contaminated_dep = paste0(samples_normal_contaminated, '_DEP')
samples_normal_contaminated_mtr = paste0(samples_normal_contaminated, '_MTR')
samples_normal_contaminated_vaf = paste0(samples_normal_contaminated, '_VAF')

twins_filtered_dt[, sum_normal_contaminated_mtr := rowSums(.SD>=4), .SDcols = samples_normal_contaminated_mtr]
twins_filtered_dt[, sum_normal_contaminated_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_contaminated_vaf]
twins_filtered_dt[, sum_normal_contaminated_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_contaminated_mtr', 'sum_normal_contaminated_vaf')]

twins_filtered_dt[, agg_normal_contaminated_mtr := rowSums(.SD), .SDcols = samples_normal_contaminated_mtr]
twins_filtered_dt[, agg_normal_contaminated_dep := rowSums(.SD), .SDcols = samples_normal_contaminated_dep]
twins_filtered_dt[, agg_normal_contaminated_vaf := agg_normal_contaminated_mtr / agg_normal_contaminated_dep]

######################################################################################################
# Calculate confidence intervals for the values of VAF 

# CI VAF for each sample separately 
for (sample in samples_names){
  sample_dep = paste0(sample, '_DEP')
  sample_mtr = paste0(sample, '_MTR')
  twins_filtered_dt[, glue('{sample}_VAF_lowerCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, get(sample_dep), get(sample_mtr)/get(sample_dep))]
  twins_filtered_dt[, glue('{sample}_VAF_upperCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, get(sample_dep), get(sample_mtr)/get(sample_dep))]
}

# CI VAF for aggregate VAF values 
twins_filtered_dt[, vaf_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all, mtr_all/dep_all)]
twins_filtered_dt[, vaf_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all, mtr_all/dep_all)]

twins_filtered_dt[, vaf_all_normal_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all_normal, mtr_all_normal/dep_all_normal)]
twins_filtered_dt[, vaf_all_normal_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all_normal, mtr_all_normal/dep_all_normal)]
twins_filtered_dt[, vaf_all_normal_PD62341_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all_normal_PD62341, mtr_all_normal_PD62341/dep_all_normal_PD62341)]
twins_filtered_dt[, vaf_all_normal_PD62341_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all_normal_PD62341, mtr_all_normal_PD62341/dep_all_normal_PD62341)]
twins_filtered_dt[, vaf_all_normal_PD63383_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all_normal_PD63383, mtr_all_normal_PD63383/dep_all_normal_PD63383)]
twins_filtered_dt[, vaf_all_normal_PD63383_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all_normal_PD63383, mtr_all_normal_PD63383/dep_all_normal_PD63383)]

twins_filtered_dt[, vaf_all_tumour_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all_tumour, mtr_all_tumour/dep_all_tumour)]
twins_filtered_dt[, vaf_all_tumour_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all_tumour, mtr_all_tumour/dep_all_tumour)]
twins_filtered_dt[, vaf_all_tumour_PD62341_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all_tumour_PD62341, mtr_all_tumour_PD62341/dep_all_tumour_PD62341)]
twins_filtered_dt[, vaf_all_tumour_PD62341_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all_tumour_PD62341, mtr_all_tumour_PD62341/dep_all_tumour_PD62341)]
twins_filtered_dt[, vaf_all_tumour_PD63383_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, dep_all_tumour_PD63383, mtr_all_tumour_PD63383/dep_all_tumour_PD63383)]
twins_filtered_dt[, vaf_all_tumour_PD63383_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, dep_all_tumour_PD63383, mtr_all_tumour_PD63383/dep_all_tumour_PD63383)]

twins_filtered_dt[, agg_normal_contaminated_vaf_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, agg_normal_contaminated_dep, agg_normal_contaminated_mtr/agg_normal_contaminated_dep)]
twins_filtered_dt[, agg_normal_contaminated_vaf_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, agg_normal_contaminated_dep, agg_normal_contaminated_mtr/agg_normal_contaminated_dep)]
twins_filtered_dt[, agg_normal_PD62341_clean_vaf_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, agg_normal_PD62341_clean_dep, agg_normal_PD62341_clean_mtr/agg_normal_PD62341_clean_dep)]
twins_filtered_dt[, agg_normal_PD62341_clean_vaf_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, agg_normal_PD62341_clean_dep, agg_normal_PD62341_clean_mtr/agg_normal_PD62341_clean_dep)]
twins_filtered_dt[, agg_normal_PD63383_clean_vaf_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, agg_normal_PD63383_clean_dep, agg_normal_PD63383_clean_mtr/agg_normal_PD63383_clean_dep)]
twins_filtered_dt[, agg_normal_PD63383_clean_vaf_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, agg_normal_PD63383_clean_dep, agg_normal_PD63383_clean_mtr/agg_normal_PD63383_clean_dep)]

######################################################################################################
# Identifying mutations where VAF can be affected by PCR duplicates
# Calculate the VAF for forward and reverse strand separately 

# maybe first write column with mutant and wild type forward and reverse
forward_cols = grep("_FAZ|_FCZ|_FTZ|_FGZ", names(twins_filtered_dt), value = TRUE) 
reverse_cols = grep("_RAZ|_RCZ|_RTZ|_RGZ", names(twins_filtered_dt), value = TRUE) 

samples_forward = c(paste0(samples_names, c('_FAZ', '_FCZ', '_FGZ', '_FTZ')))
samples_reverse = c(paste0(samples_names, c('_RAZ', '_RCZ', '_RGZ', '_RTZ')))

for (sample in samples_names){
  
  twins_filtered_dt[, glue('forward_mut_{sample}') := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
    alt = Alt[row]
    cols_to_sum = paste0(sample, '_F', alt, 'Z')
    sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
  twins_filtered_dt[, glue('reverse_mut_{sample}') := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
    alt = Alt[row]
    cols_to_sum = paste0(sample, '_R', alt, 'Z')
    sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
  
  twins_filtered_dt[, glue('forward_wt_{sample}') := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
    ref = Ref[row]
    cols_to_sum = paste0(sample, '_F', ref, 'Z')
    sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
  twins_filtered_dt[, glue('reverse_wt_{sample}') := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
    ref = Ref[row]
    cols_to_sum = paste0(sample, '_R', ref, 'Z')
    sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
  
  twins_filtered_dt[, glue('forward_vaf_{sample}') := get(glue('forward_mut_{sample}')) / (get(glue('forward_mut_{sample}')) + get(glue('forward_wt_{sample}')))]
  twins_filtered_dt[, glue('reverse_vaf_{sample}') := get(glue('reverse_mut_{sample}')) / (get(glue('reverse_mut_{sample}')) + get(glue('reverse_wt_{sample}')))] 
  
  twins_filtered_dt[,  glue('forward_vaf_{sample}_lowerCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, get(glue('forward_mut_{sample}'))+get(glue('forward_wt_{sample}')), get(glue('forward_mut_{sample}'))/(get(glue('forward_mut_{sample}'))+get(glue('forward_wt_{sample}'))))]
  twins_filtered_dt[,  glue('forward_vaf_{sample}_upperCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, get(glue('forward_mut_{sample}'))+get(glue('forward_wt_{sample}')), get(glue('forward_mut_{sample}'))/(get(glue('forward_mut_{sample}'))+get(glue('forward_wt_{sample}'))))]

  twins_filtered_dt[,  glue('reverse_vaf_{sample}_lowerCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, get(glue('reverse_mut_{sample}'))+get(glue('reverse_wt_{sample}')), get(glue('reverse_mut_{sample}'))/(get(glue('reverse_mut_{sample}'))+get(glue('reverse_wt_{sample}'))))]
  twins_filtered_dt[,  glue('reverse_vaf_{sample}_upperCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, get(glue('reverse_mut_{sample}'))+get(glue('reverse_wt_{sample}')), get(glue('reverse_mut_{sample}'))/(get(glue('reverse_mut_{sample}'))+get(glue('reverse_wt_{sample}'))))]
}

# aggregate tumour 
samples_tumour_forward = c(paste0(samples_tumour, '_FAZ'), paste0(samples_tumour, '_FCZ'),
                           paste0(samples_tumour, '_FGZ'), paste0(samples_tumour, '_FTZ'))
samples_tumour_reverse = c(paste0(samples_tumour, '_RAZ'), paste0(samples_tumour, '_RCZ'),
                           paste0(samples_tumour, '_RGZ'), paste0(samples_tumour, '_RTZ'))

samples_normal_forward = c(paste0(samples_normal, '_FAZ'), paste0(samples_normal, '_FCZ'),
                           paste0(samples_normal, '_FGZ'), paste0(samples_normal, '_FTZ'))
samples_normal_reverse = c(paste0(samples_normal, '_RAZ'), paste0(samples_normal, '_RCZ'),
                           paste0(samples_normal, '_RGZ'), paste0(samples_normal, '_RTZ'))

samples_normal_PD62341_forward = c(paste0(samples_normal_PD62341, '_FAZ'), paste0(samples_normal_PD62341, '_FCZ'),
                                   paste0(samples_normal_PD62341, '_FGZ'), paste0(samples_normal_PD62341, '_FTZ'))
samples_normal_PD62341_reverse = c(paste0(samples_normal_PD62341, '_RAZ'), paste0(samples_normal_PD62341, '_RCZ'),
                                   paste0(samples_normal_PD62341, '_RGZ'), paste0(samples_normal_PD62341, '_RTZ'))

samples_normal_PD63383_forward = c(paste0(samples_normal_PD63383, '_FAZ'), paste0(samples_normal_PD63383, '_FCZ'),
                                   paste0(samples_normal_PD63383, '_FGZ'), paste0(samples_normal_PD63383, '_FTZ'))
samples_normal_PD63383_reverse = c(paste0(samples_normal_PD63383, '_RAZ'), paste0(samples_normal_PD63383, '_RCZ'),
                                   paste0(samples_normal_PD63383, '_RGZ'), paste0(samples_normal_PD63383, '_RTZ'))


twins_filtered_dt[, forward_mtr := rowSums(.SD), .SDcols = ]


######################################################################################################
# Analysis of QC-validated mutations 

# Plot the heatmap showing all mutations of interest 
mut_all_qc = as.matrix(twins_filtered_dt[, c(samples_vaf), with=FALSE])
rownames(mut_all_qc) = twins_filtered_dt[,1] %>% unlist()  
colnames(mut_all_qc) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame('Tumour cell fraction' = c(0, 0.1, 0.2, 0, 0.9, 0.3, 0.6, 0.5, 0.8, 0.8, 0, 0, 0, 0, 0, 0.9, 0.9, 0.3, 0.7, 0.4, 0, 0.5), # fraction of tumour cells
                            Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')),
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_all_qc)
annotation_colors = list('Tumour.cell.fraction' = colorRampPalette(c('#f2e7e7', 'darkred'))(11),
                         Status = c(normal=col_normal, tumour=col_tumour),
                         Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))
                        

# heatmap
pdf('Results/20241119_p4_heatmap_all_mutations_val258.pdf')
pheatmap(mut_all_qc,
         cellwidth=12, cellheight=0.4,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations, QC validated (258)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

pdf('Results/20241119_p4_heatmap_all_mutations_val258_large.pdf', height = 45)
pheatmap(mut_all_qc,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations, QC validated (259)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Trinucleotide context plots for mutations which are in the final 259 QC-validated set

# Function to identify the trinucleotide context of a given mutation, knowing its genomic cooordinates 
get_trinucs <- function(mybed, genome) {
  
  mybed$order <- 1:nrow(mybed) # order that rows are supplied in 
  gr <- GRanges(seqnames = mybed$Chrom, IRanges(start = mybed$Pos, width=1), ref=mybed$Ref, alt=mybed$Alt, order=mybed$order)
  # create a GRanges object with coordinates for the mutation
  
  if (all(substr(seqlevels(gr), 1, 3) != "chr")) {
    gr <- renameSeqlevels(gr, paste0("chr", seqlevels(gr)))
  }
  
  seqinfo(gr) <- seqinfo(genome)[seqlevels(gr)]
  gr <- sort(gr)
  bases <- c("A", "C", "G", "T")
  trinuc_levels <- paste0(rep(bases, each = 16), rep(rep(bases, each = 4), 4), rep(bases, 16))
  get_trinuc <- function(seqname) {
    pos <- start(gr[seqnames(gr) == seqname])
    view <- Views(genome[[seqname]], start = pos - 1, end = pos + 1)
    ans <- factor(as.character(view), levels = trinuc_levels, labels = 1:64)
    return(as.numeric(ans))
  }
  trinuc <- sapply(seqlevels(gr), get_trinuc)
  gr$trinuc <- factor(unlist(trinuc, use.names = FALSE), levels = 1:64, labels = trinuc_levels)
  remove(trinuc)
  gr$REF <- gr$ref
  gr$ALT <- gr$alt
  gr$context <- gr$trinuc
  torc <- which(gr$ref %in% c("A", "G"))
  gr$REF[torc] <- as.character(reverseComplement(DNAStringSet(gr$REF[torc])))
  gr$ALT[torc] <- as.character(reverseComplement(DNAStringSet(gr$ALT[torc])))
  gr$context[torc] <- as.character(reverseComplement(DNAStringSet(gr$context[torc])))
  gr$class <- paste(gr$REF, gr$ALT, "in", gr$context, sep = ".")
  class_levels <- paste(rep(c("C", "T"), each = 48), rep(c("A", "G", "T", "A", "C", "G"), each = 16), "in", paste0(rep(rep(bases, each = 4), 6), rep(c("C", "T"), each = 48), rep(bases, 24)), sep = ".")
  gr$class <- factor(gr$class, levels = class_levels)
  
  # match them up with one another using the original order that I put in - I think that they may have been reshuffled.
  grdf <- as.data.frame(gr)
  grdf <- grdf[with(grdf, order(order)),]
  return(grdf$class)
} 

mybed_pass = twins_dt[mut_ID %in% muts, c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_pass = get_trinucs(mybed_pass, BSgenome.Hsapiens.UCSC.hg38)
dt_pass = twins_dt[mut_ID %in% muts]
dt_pass$trins_pass=trins_pass

mut_sign_counts_pass = data.table(table(dt_pass[, trins_pass]))
setnames(mut_sign_counts_pass, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_pass[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_pass[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_pass[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_pass, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  ylim(c(0, 25))+
  labs(x = 'Context', y = 'Count', title = 'Final set of mutations (259)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", linewidth = 0.1)
ggsave('Results/20241119_p4_mut_trins_259muts_pass.pdf', width = 7.5, height = 3.5)

######################################################################################################
# Identify classes of mutations of interest 

# Mutations present in all normal samples
muts_all_normal = twins_filtered_dt[sum_normal_mtr_vaf>=9, mut_ID] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_all_normal, c('mut_ID', samples_vaf, 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf', 'vaf_all_tumour'), with=FALSE]
paste('Number of mutations present in all normal samples:', length(muts_all_normal)) # 5

# Mutations specific to PD62341 twin 
muts_PD62341_normal = twins_filtered_dt[sum_normal_PD63383_clean_mtr_vaf==0 & sum_normal_PD62341_mtr_vaf >= 4, mut_ID] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_PD62341_normal, c('mut_ID', 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf', 'vaf_all_tumour'), with=FALSE]
paste('Number of mutations specific to PD62341 twin:', length(muts_PD62341_normal)) # 7

# Mutations specific to PD63383 normal samples
muts_PD63383_normal = twins_filtered_dt[sum_normal_PD63383_mtr_vaf>=3 & sum_normal_PD62341_mtr_vaf <= 1, mut_ID] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_PD63383_normal, c('mut_ID', 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf', 'vaf_all_tumour'), with=FALSE]
paste('Number of mutations specific to PD63383 twin:', length(muts_PD63383_normal)) # 11

# Mutations specific to the tumour 
muts_tumour = twins_filtered_dt[sum_tumour_mtr_vaf >=3 & sum_normal_contaminated_mtr_vaf <= 3 & sum_normal_mtr_vaf <= 3, mut_ID] %>% unlist()
paste('Number of mutations present in tumour samples but not normal samples:', length(muts_tumour)) # 87

# Mutations present in 2 PD62341 tumour samples 
muts_tumour_PD62341 = twins_filtered_dt[sum_tumour_PD62341_mtr_vaf == 2 & sum_normal_contaminated_mtr_vaf <= 1 & sum_tumour_mtr_vaf == 2, mut_ID] %>% unlist()
paste('Number of mutations specific to PD62341 tumour:', length(muts_tumour_PD62341)) # 11

# Mutations specific to PD63383 tumour samples
muts_tumour_PD63383 = twins_filtered_dt[sum_tumour_PD63383_mtr_vaf == 2 & sum_normal_contaminated_mtr_vaf <= 1 & sum_tumour_mtr_vaf == 2, mut_ID] %>% unlist()
paste('Number of mutations specific to PD63383 tumour:', length(muts_tumour_PD63383)) # 24

# Mutations present in 1 sample only
muts_one_sample_normal = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_mtr_vaf==1, mut_ID] %>% unlist() # 22
muts_one_sample_tumour = twins_filtered_dt[sum_tumour_mtr_vaf==1 & sum_normal_mtr_vaf==0, mut_ID] %>% unlist() # 88 
muts_one_sample = c(muts_one_sample_normal, muts_one_sample_tumour) 
paste('Number of mutations present in only 1 sample:', length(muts_one_sample)) # 110

######################################################################################################
# Plot distribution of VAF for mutations of interest 

twins_vaf_melt = melt(twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE], id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_vaf_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v',
  'PD63383ae', 'PD63383ak', 'PD63383bb', 'PD63383t', 'PD63383u', 'PD63383w',
  'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u',
  'PD63383ap', 'PD63383aq'))]
twins_vaf_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_melt[, status := factor(fcase(sample %in% samples_normal, 'normal', sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, sample_type := factor(fcase(status == 'tumour', 'tumour', 
                                             status == 'normal' & twin == 'PD62341', 'PD62341 normal',
                                             status == 'normal' & twin == 'PD63383', 'PD63383 normal'))]

for (mut in muts_all_normal){
  
  ggplot(twins_vaf_melt[mut_ID==mut], aes(x = sample, y = value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'Sample', y = 'VAF', col = 'Sample')+
    ylim(c(0, 0.75))+
    ggtitle(glue('{mut}'))+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241114_p4_allNormal_dist_samples_{mut}.pdf'), height = 3, width = 6.5)
  
}

######################################################################################################
# create specific heatmaps for each class of mutation

# Mutations present in all normal samples
mut_all_normal = as.matrix(twins_filtered_dt[mut_ID %in% muts_all_normal, c(samples_vaf), with=FALSE])
rownames(mut_all_normal) = twins_filtered_dt[mut_ID %in% muts_all_normal,1] %>% unlist()  
colnames(mut_all_normal) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_allnormal4.pdf')
pheatmap(mut_all_normal,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations found in all normal samples", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Mutations specific to PD62341 twin 
mut_PD62341_normal = as.matrix(twins_filtered_dt[mut_ID %in% muts_PD62341_normal, c(samples_vaf), with=FALSE])
rownames(mut_PD62341_normal) = twins_filtered_dt[mut_ID %in% muts_PD62341_normal,1] %>% unlist()  
colnames(mut_PD62341_normal) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_allPD62341.pdf')
pheatmap(mut_PD62341_normal,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations specific to PD62341 twin", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Mutations specific to PD63383
mut_PD63383_normal = as.matrix(twins_filtered_dt[mut_ID %in% muts_PD63383_normal, c(samples_vaf), with=FALSE])
rownames(mut_PD63383_normal) = twins_filtered_dt[mut_ID %in% muts_PD63383_normal,1] %>% unlist()  
colnames(mut_PD63383_normal) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_allPD63383.pdf')
pheatmap(mut_PD63383_normal,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations specific to PD63383 twin", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Mutations specific to the tumour 
mut_tumour = as.matrix(twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383), c(samples_vaf), with=FALSE])
rownames(mut_tumour) = twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383),1] %>% unlist()  
colnames(mut_tumour) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_alltumour.pdf', height = 25)
pheatmap(mut_tumour,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations specific to the tumour", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

######################################################################################################
# Tumour tree

twins_filtered_tumour = twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383)]
twins_filtered_tumour[, mut_type := factor(fcase(
                        mut_ID %in% muts_tumour_PD62341, 'PD62341 subclone',
                        mut_ID %in% muts_tumour_PD63383, 'PD63383 subclone',
                        sum_tumour_vaf >= 8, 'likely MRCA',
                        sum_tumour_vaf < 8 & sum_tumour_vaf >= 3, 'likely subclonal'))]

ggplot(twins_filtered_tumour, aes(x = dep_all_tumour, y = vaf_all_tumour, col = mut_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage', y = 'VAF', col = 'Mutation category')+
  ylim(c(0, 0.6))+
  ggtitle(glue('Tumour mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_byGroup.pdf'), height = 3, width = 6.5)

ggplot(twins_filtered_tumour, aes(x = dep_all_tumour, y = vaf_all_tumour, col = factor(sum_tumour_vaf)))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage', y = 'VAF', col = 'Number of tumour samples')+
  ylim(c(0, 0.5))+
  ggtitle(glue('Tumour mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_byNrSamples.pdf'), height = 3, width = 6.5)

# What is the distribution of number of tumour samples with mutation?
table(twins_filtered_tumour[,sum_tumour_vaf])
# 27 and 35 in 10 or 9 samples
# conveniently you are using 27 mutations to estimate tumour purity
# would be good to check if these are the same ones 

# which tumour samples are mutations typically missing from?
colSums(twins_filtered_tumour[sum_tumour_vaf==9, c(samples_tumour_vaf), with=FALSE] < 0.1) # 33
# 32 missing from ak, 1 missing from b, 2 missing from ag 
# 11 are at VAF 0 at ak, so likely they really are not there (not just contamination)

colSums(twins_filtered_tumour[sum_tumour_vaf==8, c(samples_tumour_vaf), with=FALSE] < 0.1) # 4
# 4 missing from ak, 3 missing from ag, 1 missing from u

colSums(twins_filtered_tumour[sum_tumour_vaf==7, c(samples_tumour_vaf), with=FALSE] < 0.1) # 1
# 1 missing from ag, ak and u

colSums(twins_filtered_tumour[sum_tumour_vaf==6, c(samples_tumour_vaf), with=FALSE] < 0.1) # 6
# all missing in ak and ag
# 5 missing in u
# 4 missing in ae
# 3 missing in aj

colSums(twins_filtered_tumour[sum_tumour_vaf==5, c(samples_tumour_vaf), with=FALSE] > 0.1) # 7
# 0.4 in PD63383, 0.15 am, 0.1 ae, 0.16 b
# 0.4 in PD63383, 0.15 am, 0.15 ae, 0.1 ap
# 0.4 in PD63383, 0.15 am, 0.12 ae, 0.16 b
# 0.3 ae, 0.16 ag, 0.19 aj, 0.15 b, 0.12 u
# 0.3-0.4 in PD63383, 0.3 am, 0.2 b, 0.12 u
# 0.4 in PD63383, 0.1 am, 0.17 b, 0.1 ae
# 0.3 ae, 0.14 ag, 0.4 aj, 0.15 b, 0.1 u

colSums(twins_filtered_tumour[sum_tumour_vaf==4, c(samples_tumour_vaf), with=FALSE] > 0.1) # 5
# ae, aj, b, u (0.13 in b, 0.2 in the rest, 0.3 in ae), 0.7 in ag
# am, b, PD63383 (0.1 in am, 0.15 in b, 0.4 in PD63383)
# ae 0.4, ag 0.2, aj 0.2, 0.1 b, 0.09 in u
# am 0.2, 0.5 in PD63383, 0.2 in b, 0.07 in u
# am 0.14, 0.5 in PD63383, 0.13 in b, 0.08 in u

colSums(twins_filtered_tumour[sum_tumour_vaf==3, c(samples_tumour_vaf), with=FALSE] > 0.1) # 2
# ae, aj, ag (0.4 in ae, 0.2 in aj, 0.1 in ag)
# ap, aq, am (> 0.1 in am, ~0.5 in ap and aq)

ggplot(twins_filtered_tumour, aes(x = vaf_all_tumour, y = PD62341ak_VAF, col = factor(sum_tumour_vaf)))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF across the tumour', y = 'VAF PD62341ak', col = 'Number of tumour samples')+
  ylim(c(0, 0.5))+
  xlim(c(0, 0.5))+
  coord_equal(ratio=1)+
  ggtitle(glue('Tumour mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

twins_filtered_tumour[, agg_tumour_ak_mtr := rowSums(.SD), .SDcols = samples_tumour_mtr[samples_tumour_mtr!='PD62341ak_MTR']]
twins_filtered_tumour[, agg_tumour_ak_dep := rowSums(.SD), .SDcols = samples_tumour_dep[samples_tumour_mtr!='PD62341ak_DEP']]
twins_filtered_tumour[, agg_tumour_ak_vaf := agg_tumour_ak_mtr / agg_tumour_ak_dep]

ggplot(twins_filtered_tumour, aes(x = agg_tumour_ak_vaf, y = PD62341ak_VAF, col = factor(sum_tumour_vaf)))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF across the tumour\n(excluding PD62341ak)', y = 'VAF PD62341ak', col = 'Number of tumour samples')+
  ylim(c(0, 0.5))+
  xlim(c(0, 0.5))+
  coord_equal(ratio=1)+
  ggtitle(glue('Tumour mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# can I plot the VAF values for each sample for different nr of tumours
twins_tumour_melt = melt(twins_filtered_tumour[, c('mut_ID', samples_vaf, 'sum_tumour_vaf'), with=FALSE], id.vars = c('mut_ID', 'sum_tumour_vaf'))
twins_tumour_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_tumour_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'
))]
twins_tumour_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'
))]
twins_tumour_melt[, sample_type := paste(status, twin, sep = ', ')]

for (num in seq(2:11)){
  ggplot(twins_tumour_melt[status=='tumour' & sum_tumour_vaf==num], aes(x = sample, y = value))+
    geom_point(size=1.5, alpha = 0.8)+
    theme_classic(base_size = 14)+
    labs(x = 'Sample', y = 'VAF')+
    ylim(c(0, 0.7))+
    ggtitle(glue('Number of samples with mutation = {num}'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_tumourMuts_vafDist_{num}samples.pdf'), height = 3, width = 7)
}

# add sample tumour cell fraction to this
purity_dt[, sample := tstrsplit(sample, '_VAF', keep=1, fixed = TRUE)]
twins_tumour_melt = merge(twins_tumour_melt, purity_dt[, c('sample', 'tumour_cell_fraction'), with=FALSE])

for (num in seq(3:12)){
  ggplot(twins_tumour_melt[status=='tumour' & sum_tumour_vaf==num], aes(x = tumour_cell_fraction, y = value))+
    geom_point(size=1.5, alpha = 0.8)+
    theme_classic(base_size = 11)+
    labs(x = 'Tumour cell fraction', y = 'VAF')+
    ylim(c(0, 0.7))+
    ggtitle(glue('Mutations present in {num} samples'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_{num}samples.pdf'), height = 3, width = 3.5)
}

for (num in seq(3:12)){
  ggplot(twins_tumour_melt[status=='tumour' & sum_tumour_vaf==num], aes(x = tumour_cell_fraction, y = value, col = factor(sample)))+
    geom_point(size=1.5, alpha = 0.8)+
    theme_classic(base_size = 10)+
    labs(x = 'Tumour cell fraction', y = 'VAF', col = 'Sample')+
    ylim(c(0, 0.7))+
    ggtitle(glue('Mutations present in {num} samples'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))
  ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_{num}samples_colSample.pdf'), height = 3.5, width = 4.5)
}

for (num in seq(3:12)){
  ggplot(twins_tumour_melt[status=='tumour' & sum_tumour_vaf==num], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID)))+
    geom_point(size=1.5, alpha = 0.8)+
    guides(color = "none")+
    geom_line()+
    theme_classic(base_size = 10)+
    labs(x = 'Tumour cell fraction', y = 'VAF', col = 'Sample')+
    ylim(c(0, 0.7))+
    ggtitle(glue('Mutations present in {num} samples'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_{num}samples_colSample_lines.pdf'), height = 3.5, width = 4.5)
}

ggplot(twins_tumour_melt[status=='tumour'], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID)))+
  geom_point(size=1.5, alpha = 0.8)+
  guides(color = "none")+
  geom_line()+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Mutations present in {num} samples'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# it feels like if you treated this as a time course you could deconvolve this to clones
# okay first do heatmap and cluster rows 
mut_tumour = as.matrix(twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383), c(samples_vaf), with=FALSE])
rownames(mut_tumour) = twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383),1] %>% unlist()  
colnames(mut_tumour) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()

pdf('Results/20241119_p4_heatmap_muts_alltumour.pdf', height = 25)
pheatmap(mut_tumour,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations specific to the tumour", 
         legend = T, 
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# okay maybe reshape and do clustering on this
twins_tumour_data = twins_filtered_tumour[, c('mut_ID', samples_vaf), with=FALSE]

for (k in seq(from = 2, to = 10)){
  
  kmeans_clusters = kmeans(twins_tumour_data[, c(samples_vaf), with=FALSE], centers = k)
  twins_tumour_data[, cluster := kmeans_clusters$cluster]
  
  twins_tumour_melt = melt(twins_tumour_data[, c('mut_ID', samples_vaf, 'cluster'), with=FALSE], id.vars = c('mut_ID', 'cluster'))
  twins_tumour_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
  twins_tumour_melt[, status := factor(fcase(
    sample %in% samples_normal, 'normal',
    sample %in% samples_tumour, 'tumour'))]
  twins_tumour_melt = merge(twins_tumour_melt, purity_dt[, c('sample', 'tumour_cell_fraction'), with=FALSE])
  
  ggplot(twins_tumour_melt[status=='tumour'], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = factor(cluster)))+
    geom_point(size=1.5, alpha = 0.8)+
    geom_line()+
    theme_classic(base_size = 10)+
    labs(x = 'Tumour cell fraction', y = 'VAF', col = 'Cluster')+
    ylim(c(0, 0.7))+
    ggtitle(glue('Tumour samples, k-means clustering with {k} clusters'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans_{k}.pdf'), height = 3.5, width = 4.5)
}

# plot clusters separately to see if they are different enough 
kmeans_clusters = kmeans(twins_tumour_data[, c(samples_vaf), with=FALSE], centers = 2)
twins_tumour_data[, cluster := kmeans_clusters$cluster]
  
twins_tumour_melt = melt(twins_tumour_data[, c('mut_ID', samples_vaf, 'cluster'), with=FALSE], id.vars = c('mut_ID', 'cluster'))
twins_tumour_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_tumour_melt[, status := factor(fcase(
    sample %in% samples_normal, 'normal',
    sample %in% samples_tumour, 'tumour'))]
twins_tumour_melt = merge(twins_tumour_melt, purity_dt[, c('sample', 'tumour_cell_fraction'), with=FALSE])
  
ggplot(twins_tumour_melt[status=='tumour'&cluster==1], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID)))+
    geom_point(size=1.5, alpha = 0.8)+
    geom_line()+
    theme_classic(base_size = 10)+
    labs(x = 'Tumour cell fraction', y = 'VAF')+
    ylim(c(0, 0.7))+
    ggtitle(glue('Tumour samples, cluster 1 (kmeans 2)'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans2_cluster1.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==2], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID)))+
  geom_point(size=1.5, alpha = 0.8)+
  geom_line()+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 2 (kmeans 2)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans2_cluster2.pdf'), height = 3.5, width = 4.5)

# do this for 3 clusters
kmeans_clusters = kmeans(twins_tumour_data[, c(samples_vaf), with=FALSE], centers = 3)
twins_tumour_data[, cluster := kmeans_clusters$cluster]

twins_tumour_melt = melt(twins_tumour_data[, c('mut_ID', samples_vaf, 'cluster'), with=FALSE], id.vars = c('mut_ID', 'cluster'))
twins_tumour_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_tumour_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_tumour_melt = merge(twins_tumour_melt, purity_dt[, c('sample', 'tumour_cell_fraction'), with=FALSE])

ggplot(twins_tumour_melt[status=='tumour'&cluster==1], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID)))+
  geom_point(size=1.5, alpha = 0.8)+
  geom_line()+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 1 (kmeans 3)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans3_cluster1.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==2], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID)))+
  geom_point(size=1.5, alpha = 0.8)+
  geom_line()+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 2 (kmeans 3)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans3_cluster2.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==3], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID)))+
  geom_point(size=1.5, alpha = 0.8)+
  geom_line()+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 3 (kmeans 3)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans3_cluster3.pdf'), height = 3.5, width = 4.5)

# do this for 4 clusters
set.seed(10)
kmeans_clusters = kmeans(twins_tumour_data[, c(samples_vaf), with=FALSE], centers = 4)
twins_tumour_data[, cluster := kmeans_clusters$cluster]

twins_tumour_melt = melt(twins_tumour_data[, c('mut_ID', samples_vaf, 'cluster'), with=FALSE], id.vars = c('mut_ID', 'cluster'))
twins_tumour_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_tumour_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_tumour_melt = merge(twins_tumour_melt, purity_dt[, c('sample', 'tumour_cell_fraction'), with=FALSE])

ggplot(twins_tumour_melt[status=='tumour'&cluster==1], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 1 (kmeans 4)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans4_cluster1.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==2], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 2 (kmeans 4)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans4_cluster2.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==3], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 3 (kmeans 4)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans4_cluster3.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==4], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ylim(c(0, 0.7))+
  ggtitle(glue('Tumour samples, cluster 4 (kmeans 4)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans4_cluster4.pdf'), height = 3.5, width = 4.5)

# okay try 5 clusters but no more
# do this for 4 clusters
set.seed(10)
kmeans_clusters = kmeans(twins_tumour_data[, c(samples_vaf), with=FALSE], centers = 5)
twins_tumour_data[, cluster := kmeans_clusters$cluster]

twins_tumour_melt = melt(twins_tumour_data[, c('mut_ID', samples_vaf, 'cluster'), with=FALSE], id.vars = c('mut_ID', 'cluster'))
twins_tumour_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_tumour_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_tumour_melt = merge(twins_tumour_melt, purity_dt[, c('sample', 'tumour_cell_fraction'), with=FALSE])

ggplot(twins_tumour_melt[status=='tumour'&cluster==1], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ggtitle(glue('Tumour samples, cluster 1 (kmeans 5)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans5_cluster1.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==2], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ggtitle(glue('Tumour samples, cluster 2 (kmeans 5)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans5_cluster2.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==3], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ggtitle(glue('Tumour samples, cluster 3 (kmeans 5)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans5_cluster3.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==4], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ggtitle(glue('Tumour samples, cluster 4 (kmeans 5)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans5_cluster4.pdf'), height = 3.5, width = 4.5)

ggplot(twins_tumour_melt[status=='tumour'&cluster==5], aes(x = tumour_cell_fraction, y = value, group = factor(mut_ID), col = sample))+
  geom_line(col = 'grey')+
  geom_point(size=2.5, alpha = 1)+
  theme_classic(base_size = 10)+
  labs(x = 'Tumour cell fraction', y = 'VAF')+
  ggtitle(glue('Tumour samples, cluster 5 (kmeans 5)'))+
  scale_color_manual(values = c('#09ddbd', '#167288', '#a89a49', '#d48c84', '#dd0d0d', '#8cdaec', '#dd8709', '#9bddb1', '#d2ace0', '#6f0b9e'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_tumourMuts_tumourCellFraction_vs_Vaf_clustering_kmeans5_cluster4.pdf'), height = 3.5, width = 4.5)








