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
paste('Number of mutations that passed QC:', length(muts)) # 256 

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
  twins_filtered_dt[, glue('{sample}_vafLowerCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, get(sample_dep), get(sample_mtr)/get(sample_dep))]
  twins_filtered_dt[, glue('{sample}_vafUpperCI') := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, get(sample_dep), get(sample_mtr)/get(sample_dep))]
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

twins_filtered_dt[, agg_normal_contaminated_vafLowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, agg_normal_contaminated_dep, agg_normal_contaminated_mtr/agg_normal_contaminated_dep)]
twins_filtered_dt[, agg_normal_contaminated_vafUpperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, agg_normal_contaminated_dep, agg_normal_contaminated_mtr/agg_normal_contaminated_dep)]
twins_filtered_dt[, agg_normal_PD62341_clean_vafLowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, agg_normal_PD62341_clean_dep, agg_normal_PD62341_clean_mtr/agg_normal_PD62341_clean_dep)]
twins_filtered_dt[, agg_normal_PD62341_clean_vafUpperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, agg_normal_PD62341_clean_dep, agg_normal_PD62341_clean_mtr/agg_normal_PD62341_clean_dep)]
twins_filtered_dt[, agg_normal_PD63383_clean_vafLowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, agg_normal_PD63383_clean_dep, agg_normal_PD63383_clean_mtr/agg_normal_PD63383_clean_dep)]
twins_filtered_dt[, agg_normal_PD63383_clean_vafUpperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, agg_normal_PD63383_clean_dep, agg_normal_PD63383_clean_mtr/agg_normal_PD63383_clean_dep)]

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
  
  twins_filtered_dt[, glue('forward_wt_{sample}') := sapply(1:.N, function(row) { # identify forward strands which are wt 
    ref = Ref[row]
    cols_to_sum = paste0(sample, '_F', ref, 'Z')
    sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
  twins_filtered_dt[, glue('reverse_wt_{sample}') := sapply(1:.N, function(row) { # identify reverse strands which are wt  
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

# Aggregate values for all samples
twins_filtered_dt[, forward_mut_all := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_names, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_all := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_names, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_all := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_names, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_all := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_names, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_all := forward_mut_all / (forward_mut_all + forward_wt_all)]
twins_filtered_dt[, reverse_vaf_all := reverse_mut_all / (reverse_mut_all + reverse_wt_all)] 
twins_filtered_dt[,  forward_vaf_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_all + forward_wt_all, forward_vaf_all)]
twins_filtered_dt[,  forward_vaf_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_all + forward_wt_all, forward_vaf_all)]
twins_filtered_dt[,  reverse_vaf_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_all + reverse_wt_all, reverse_vaf_all)]
twins_filtered_dt[,  reverse_vaf_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_all + reverse_wt_all, reverse_vaf_all)]

# Aggregate normal values 
twins_filtered_dt[, forward_mut_normal_all := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_normal_all := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_normal_all := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_normal_all := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_normal_all := forward_mut_normal_all / (forward_mut_normal_all + forward_wt_normal_all)]
twins_filtered_dt[, reverse_vaf_normal_all := reverse_mut_normal_all / (reverse_mut_normal_all + reverse_wt_normal_all)] 
twins_filtered_dt[,  forward_vaf_normal_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_normal_all + forward_wt_normal_all, forward_vaf_normal_all)]
twins_filtered_dt[,  forward_vaf_normal_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_normal_all + forward_wt_normal_all, forward_vaf_normal_all)]
twins_filtered_dt[,  reverse_vaf_normal_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_normal_all + reverse_wt_normal_all, reverse_vaf_normal_all)]
twins_filtered_dt[,  reverse_vaf_normal_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_normal_all + reverse_wt_normal_all, reverse_vaf_normal_all)]

# Aggregate PD62341 normal
twins_filtered_dt[, forward_mut_normal_PD62341 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD62341, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_normal_PD62341 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD62341, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_normal_PD62341 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD62341, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_normal_PD62341 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD62341, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_normal_PD62341 := forward_mut_normal_PD62341 / (forward_mut_normal_PD62341 + forward_wt_normal_PD62341)]
twins_filtered_dt[, reverse_vaf_normal_PD62341 := reverse_mut_normal_PD62341 / (reverse_mut_normal_PD62341 + reverse_wt_normal_PD62341)] 
twins_filtered_dt[,  forward_vaf_normal_PD62341_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_normal_PD62341 + forward_wt_normal_PD62341, forward_vaf_normal_PD62341)]
twins_filtered_dt[,  forward_vaf_normal_PD62341_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_normal_PD62341 + forward_wt_normal_PD62341, forward_vaf_normal_PD62341)]
twins_filtered_dt[,  reverse_vaf_normal_PD62341_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_normal_PD62341 + reverse_wt_normal_PD62341, reverse_vaf_normal_PD62341)]
twins_filtered_dt[,  reverse_vaf_normal_PD62341_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_normal_PD62341 + reverse_wt_normal_PD62341, reverse_vaf_normal_PD62341)]

# Aggregate PD63383 normal
twins_filtered_dt[, forward_mut_normal_PD63383 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD63383, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_normal_PD63383 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD63383, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_normal_PD63383 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD63383, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_normal_PD63383 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD63383, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_normal_PD63383 := forward_mut_normal_PD63383 / (forward_mut_normal_PD63383 + forward_wt_normal_PD63383)]
twins_filtered_dt[, reverse_vaf_normal_PD63383 := reverse_mut_normal_PD63383 / (reverse_mut_normal_PD63383 + reverse_wt_normal_PD63383)] 
twins_filtered_dt[,  forward_vaf_normal_PD63383_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_normal_PD63383 + forward_wt_normal_PD63383, forward_vaf_normal_PD63383)]
twins_filtered_dt[,  forward_vaf_normal_PD63383_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_normal_PD63383 + forward_wt_normal_PD63383, forward_vaf_normal_PD63383)]
twins_filtered_dt[,  reverse_vaf_normal_PD63383_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_normal_PD63383 + reverse_wt_normal_PD63383, reverse_vaf_normal_PD63383)]
twins_filtered_dt[,  reverse_vaf_normal_PD63383_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_normal_PD63383 + reverse_wt_normal_PD63383, reverse_vaf_normal_PD63383)]

# Aggregate tumour values  
twins_filtered_dt[, forward_mut_tumour_all := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_tumour, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_tumour_all := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_tumour, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_tumour_all := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_tumour, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_tumour_all := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_tumour, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_tumour_all := forward_mut_tumour_all / (forward_mut_tumour_all + forward_wt_tumour_all)]
twins_filtered_dt[, reverse_vaf_tumour_all := reverse_mut_tumour_all / (reverse_mut_tumour_all + reverse_wt_tumour_all)] 
twins_filtered_dt[,  forward_vaf_tumour_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_tumour_all + forward_wt_tumour_all, forward_vaf_tumour_all)]
twins_filtered_dt[,  forward_vaf_tumour_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_tumour_all + forward_wt_tumour_all, forward_vaf_tumour_all)]
twins_filtered_dt[,  reverse_vaf_tumour_all_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_tumour_all + reverse_wt_tumour_all, reverse_vaf_tumour_all)]
twins_filtered_dt[,  reverse_vaf_tumour_all_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_tumour_all + reverse_wt_tumour_all, reverse_vaf_tumour_all)]

# Aggregate PD62341 tumour 
twins_filtered_dt[, forward_mut_tumour_PD62341 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_tumour_PD62341, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_tumour_PD62341 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_tumour_PD62341, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_tumour_PD62341 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_tumour_PD62341, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_tumour_PD62341 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_tumour_PD62341, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_tumour_PD62341 := forward_mut_tumour_PD62341 / (forward_mut_tumour_PD62341 + forward_wt_tumour_PD62341)]
twins_filtered_dt[, reverse_vaf_tumour_PD62341 := reverse_mut_tumour_PD62341 / (reverse_mut_tumour_PD62341 + reverse_wt_tumour_PD62341)] 
twins_filtered_dt[,  forward_vaf_tumour_PD62341_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_tumour_PD62341 + forward_wt_tumour_PD62341, forward_vaf_tumour_PD62341)]
twins_filtered_dt[,  forward_vaf_tumour_PD62341_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_tumour_PD62341 + forward_wt_tumour_PD62341, forward_vaf_tumour_PD62341)]
twins_filtered_dt[,  reverse_vaf_tumour_PD62341_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_tumour_PD62341 + reverse_wt_tumour_PD62341, reverse_vaf_tumour_PD62341)]
twins_filtered_dt[,  reverse_vaf_tumour_PD62341_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_tumour_PD62341 + reverse_wt_tumour_PD62341, reverse_vaf_tumour_PD62341)]

# Aggregate PD63383 tumour 
twins_filtered_dt[, forward_mut_tumour_PD63383 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_tumour_PD63383, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_tumour_PD63383 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_tumour_PD63383, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_tumour_PD63383 := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_tumour_PD63383, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_tumour_PD63383 := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_tumour_PD63383, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_tumour_PD63383 := forward_mut_tumour_PD63383 / (forward_mut_tumour_PD63383 + forward_wt_tumour_PD63383)]
twins_filtered_dt[, reverse_vaf_tumour_PD63383 := reverse_mut_tumour_PD63383 / (reverse_mut_tumour_PD63383 + reverse_wt_tumour_PD63383)] 
twins_filtered_dt[,  forward_vaf_tumour_PD63383_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_tumour_PD63383 + forward_wt_tumour_PD63383, forward_vaf_tumour_PD63383)]
twins_filtered_dt[,  forward_vaf_tumour_PD63383_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_tumour_PD63383 + forward_wt_tumour_PD63383, forward_vaf_tumour_PD63383)]
twins_filtered_dt[,  reverse_vaf_tumour_PD63383_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_tumour_PD63383 + reverse_wt_tumour_PD63383, reverse_vaf_tumour_PD63383)]
twins_filtered_dt[,  reverse_vaf_tumour_PD63383_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_tumour_PD63383 + reverse_wt_tumour_PD63383, reverse_vaf_tumour_PD63383)]

# In addition, do this for clean samples (possibly most useful for comparisons)

# clean PD62341 (excluding PD62341v - spleen contaminated due to twin-twin transfusion)
twins_filtered_dt[, forward_mut_normal_PD62341_clean := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD62341_clean, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_normal_PD62341_clean := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD62341_clean, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_normal_PD62341_clean := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD62341_clean, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_normal_PD62341_clean := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD62341_clean, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_normal_PD62341_clean := forward_mut_normal_PD62341_clean / (forward_mut_normal_PD62341_clean + forward_wt_normal_PD62341_clean)]
twins_filtered_dt[, reverse_vaf_normal_PD62341_clean := reverse_mut_normal_PD62341_clean / (reverse_mut_normal_PD62341_clean + reverse_wt_normal_PD62341_clean)] 
twins_filtered_dt[,  forward_vaf_normal_PD62341_clean_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_normal_PD62341_clean + forward_wt_normal_PD62341_clean, forward_vaf_normal_PD62341_clean)]
twins_filtered_dt[,  forward_vaf_normal_PD62341_clean_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_normal_PD62341_clean + forward_wt_normal_PD62341_clean, forward_vaf_normal_PD62341_clean)]
twins_filtered_dt[,  reverse_vaf_normal_PD62341_clean_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_normal_PD62341_clean + reverse_wt_normal_PD62341_clean, reverse_vaf_normal_PD62341_clean)]
twins_filtered_dt[,  reverse_vaf_normal_PD62341_clean_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_normal_PD62341_clean + reverse_wt_normal_PD62341_clean, reverse_vaf_normal_PD62341_clean)]

# clean PD63383 (excluding PD63383bb - skin contaminated by PD62341-derived tumour)
twins_filtered_dt[, forward_mut_normal_PD63383_clean := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD63383_clean, '_F', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_mut_normal_PD63383_clean := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = paste0(samples_normal_PD63383_clean, '_R', alt, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]
twins_filtered_dt[, forward_wt_normal_PD63383_clean := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD63383_clean, '_F', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]  
twins_filtered_dt[, reverse_wt_normal_PD63383_clean := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  ref = Ref[row]
  cols_to_sum = paste0(samples_normal_PD63383_clean, '_R', ref, 'Z')
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)})]

twins_filtered_dt[, forward_vaf_normal_PD63383_clean := forward_mut_normal_PD63383_clean / (forward_mut_normal_PD63383_clean + forward_wt_normal_PD63383_clean)]
twins_filtered_dt[, reverse_vaf_normal_PD63383_clean := reverse_mut_normal_PD63383_clean / (reverse_mut_normal_PD63383_clean + reverse_wt_normal_PD63383_clean)] 
twins_filtered_dt[,  forward_vaf_normal_PD63383_clean_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, forward_mut_normal_PD63383_clean + forward_wt_normal_PD63383_clean, forward_vaf_normal_PD63383_clean)]
twins_filtered_dt[,  forward_vaf_normal_PD63383_clean_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, forward_mut_normal_PD63383_clean + forward_wt_normal_PD63383_clean, forward_vaf_normal_PD63383_clean)]
twins_filtered_dt[,  reverse_vaf_normal_PD63383_clean_lowerCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.05, reverse_mut_normal_PD63383_clean + reverse_wt_normal_PD63383_clean, reverse_vaf_normal_PD63383_clean)]
twins_filtered_dt[,  reverse_vaf_normal_PD63383_clean_upperCI := mapply(function(a, s, p) qbinom(a, s, p) / s, 0.95, reverse_mut_normal_PD63383_clean + reverse_wt_normal_PD63383_clean, reverse_vaf_normal_PD63383_clean)]

######################################################################################################
# Compare VAF of forward and reverse strands of mutations with suspected VAF issues

# Compare VAF of forward and reverse strands across all mutations
pdf('Results/20241119_p4_hist_log2_forward_reverse_vaf.pdf')
hist(log2(twins_filtered_dt[,forward_vaf_all] / twins_filtered_dt[,reverse_vaf_all]),
     xlab = 'log2(forward VAF/reverse VAF)', main = 'VAF across all samples (22)')
dev.off()

sum(abs(log2(twins_filtered_dt[,forward_vaf_all] / twins_filtered_dt[,reverse_vaf_all]))>1) # 42
sum(abs(log2(twins_filtered_dt[,forward_vaf_all] / twins_filtered_dt[,reverse_vaf_all]))>1.5) # 21

twins_filtered_dt[, log2_forward_reverse := log2(forward_vaf_all / reverse_vaf_all)]
twins_filtered_dt[abs(log2_forward_reverse)>1, mut_ID] %>% unlist()

sum(twins_filtered_dt[abs(log2_forward_reverse)>1, mut_ID] %>% unlist() %in% muts_all_normal) # 0
sum(twins_filtered_dt[abs(log2_forward_reverse)>1.5, mut_ID] %>% unlist() %in% muts_all_normal) # 1
sum(twins_filtered_dt[abs(log2_forward_reverse)>1, mut_ID] %>% unlist() %in% muts_PD62341_normal) # 0
sum(twins_filtered_dt[abs(log2_forward_reverse)>1, mut_ID] %>% unlist() %in% muts_PD63383_normal) # 0
sum(twins_filtered_dt[abs(log2_forward_reverse)>1.5, mut_ID] %>% unlist() %in% muts_PD63383_normal) # 1

sum(twins_filtered_dt[abs(log2_forward_reverse)>1, mut_ID] %>% unlist() %in% muts_tumour) # 0
sum(twins_filtered_dt[abs(log2_forward_reverse)>1, mut_ID] %>% unlist() %in% muts_tumour_PD62341) # 3
sum(twins_filtered_dt[abs(log2_forward_reverse)>1.5, mut_ID] %>% unlist() %in% muts_tumour_PD62341) # 1
sum(twins_filtered_dt[abs(log2_forward_reverse)>1, mut_ID] %>% unlist() %in% muts_tumour_PD63383) # 2
sum(twins_filtered_dt[abs(log2_forward_reverse)>1.5, mut_ID] %>% unlist() %in% muts_tumour_PD63383) # 0

# Show distribution of VAF values on a histogram
pdf('Results/20241119_p4_hist_reverse_forward_vaf_chr1_388.pdf', height = 8)
par(mfrow = c(2,1))
hist(twins_filtered_dt[mut_ID == 'chr1_38827952_C_A', c(paste0('forward_vaf_', samples_names)), with=FALSE] %>% unlist(), 
     breaks=10, xlab = 'VAF on forward strand (aggregate samples)', main = 'chr1_38827952_C_A, forward', xlim = c(0, 1))
hist(twins_filtered_dt[mut_ID == 'chr1_38827952_C_A', c(paste0('reverse_vaf_', samples_names)), with=FALSE] %>% unlist(), 
     breaks=10, xlab = 'VAF on reverse strand (aggregate samples)', main = 'chr1_38827952_C_A, reverse', xlim = c(0, 1))
dev.off()

# Plot VAF distribution for each sample individually (VAF from forward + reverse strands)
# Plot distribution for each sample individually 
twins_vaf_melt = data.table(melt(twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE], id.vars = 'mut_ID'))
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep =1)]
twins_vaf_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v',
  'PD63383ak', 'PD63383bb', 'PD63383t', 'PD63383ae', 'PD63383u', 'PD63383w',
  'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u',
  'PD63383ap', 'PD63383aq'))]
twins_vaf_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_melt[, status := factor(fcase(sample %in% samples_normal, 'normal', sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, sample_type := factor(fcase(status == 'tumour', 'tumour', 
                                             status == 'normal' & twin == 'PD62341', 'PD62341 normal',
                                             status == 'normal' & twin == 'PD63383', 'PD63383 normal'))]
twins_vaf_melt[, tissue := factor(fcase(
  sample %in% c('PD62341ad', 'PD63383ak'), 'cerebellum',
  sample %in% c('PD62341aa', 'PD63383bb'), 'skin',
  sample %in% c('PD62341v', 'PD63383w'), 'spleen',
  sample %in% c('PD62341h', 'PD63383t'), 'liver',
  sample %in% c('PD62341q', 'PD63383u'), 'pancreas',
  sample %in% c('PD62341n', 'PD63383ae'), 'heart', 
  sample %in% samples_tumour, 'tumour'
))]

# Plot distribution of  VAF across samples
ggplot(twins_vaf_melt[mut_ID=='chr1_38827952_C_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'F+R VAF', col = 'Sample')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_VAF_chr1_388_agg.pdf'), height = 3, width = 6.5)

# Plots by tissue 
ggplot(twins_vaf_melt[mut_ID=='chr1_38827952_C_A'&status=='normal'], aes(x = tissue, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'F+R VAF', col = 'Twin')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_VAF_chr1_388_tissue_agg.pdf'), height = 3, width = 6.5)

# Plot distribution for each sample individually (VAF from forward strand)
twins_vaf_forward_melt = data.table(melt(twins_filtered_dt[, c('mut_ID', paste0('forward_vaf_', samples_names)), with=FALSE], id.vars = 'mut_ID'))
twins_vaf_forward_melt[, sample := tstrsplit(variable, 'forward_vaf_', fixed = TRUE, keep =2)]
twins_vaf_forward_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v',
  'PD63383ak', 'PD63383bb', 'PD63383t', 'PD63383ae', 'PD63383u', 'PD63383w',
  'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u',
  'PD63383ap', 'PD63383aq'))]
twins_vaf_forward_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_forward_melt[, status := factor(fcase(sample %in% samples_normal, 'normal', sample %in% samples_tumour, 'tumour'))]
twins_vaf_forward_melt[, sample_type := factor(fcase(status == 'tumour', 'tumour', 
                                                     status == 'normal' & twin == 'PD62341', 'PD62341 normal',
                                                     status == 'normal' & twin == 'PD63383', 'PD63383 normal'))]
twins_vaf_forward_melt[, tissue := factor(fcase(
  sample %in% c('PD62341ad', 'PD63383ak'), 'cerebellum',
  sample %in% c('PD62341aa', 'PD63383bb'), 'skin',
  sample %in% c('PD62341v', 'PD63383w'), 'spleen',
  sample %in% c('PD62341h', 'PD63383t'), 'liver',
  sample %in% c('PD62341q', 'PD63383u'), 'pancreas',
  sample %in% c('PD62341n', 'PD63383ae'), 'heart', 
  sample %in% samples_tumour, 'tumour'
))]

# Plot distribution of forward VAF across samples
ggplot(twins_vaf_forward_melt[mut_ID=='chr1_38827952_C_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'forward VAF', col = 'Sample')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_forwardVAF_chr1_388.pdf'), height = 3, width = 6.5)

# Plots by tissue 
ggplot(twins_vaf_forward_melt[mut_ID=='chr1_38827952_C_A'&status=='normal'], aes(x = tissue, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'forward VAF', col = 'Twin')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_forwardVAF_chr1_388_tissue.pdf'), height = 3, width = 6.5)


# Plot distribution for each sample individually (VAF from reverse strand)
twins_vaf_reverse_melt = data.table(melt(twins_filtered_dt[, c('mut_ID', paste0('reverse_vaf_', samples_names)), with=FALSE], id.vars = 'mut_ID'))
twins_vaf_reverse_melt[, sample := tstrsplit(variable, 'reverse_vaf_', fixed = TRUE, keep =2)]
twins_vaf_reverse_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v',
  'PD63383ak', 'PD63383bb', 'PD63383t', 'PD63383ae', 'PD63383u', 'PD63383w',
  'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u',
  'PD63383ap', 'PD63383aq'))]
twins_vaf_reverse_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_reverse_melt[, status := factor(fcase(sample %in% samples_normal, 'normal', sample %in% samples_tumour, 'tumour'))]
twins_vaf_reverse_melt[, sample_type := factor(fcase(status == 'tumour', 'tumour', 
                                                     status == 'normal' & twin == 'PD62341', 'PD62341 normal',
                                                     status == 'normal' & twin == 'PD63383', 'PD63383 normal'))]
twins_vaf_reverse_melt[, tissue := factor(fcase(
  sample %in% c('PD62341ad', 'PD63383ak'), 'cerebellum',
  sample %in% c('PD62341aa', 'PD63383bb'), 'skin',
  sample %in% c('PD62341v', 'PD63383w'), 'spleen',
  sample %in% c('PD62341h', 'PD63383t'), 'liver',
  sample %in% c('PD62341q', 'PD63383u'), 'pancreas',
  sample %in% c('PD62341n', 'PD63383ae'), 'heart', 
  sample %in% samples_tumour, 'tumour'
))]

# Plot distribution of reverse VAF across samples
ggplot(twins_vaf_reverse_melt[mut_ID=='chr1_38827952_C_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'reverse VAF', col = 'Sample')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_reverseVAF_chr1_388.pdf'), height = 3, width = 6.5)

# Plots by tissue 
ggplot(twins_vaf_reverse_melt[mut_ID=='chr1_38827952_C_A'&status=='normal'], aes(x = tissue, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'reverse VAF', col = 'Twin')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241119_p4_reverseVAF_chr1_388_tissue.pdf'), height = 3, width = 6.5)
