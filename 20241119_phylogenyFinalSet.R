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
twins_dt = fread('Data/pileup_merged_20241016.tsv') # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Load QC-validated mutations (final set of mutations to be used)
muts_dt = data.table(read.csv('Data/20241114_599muts_QCJbrowse.csv', header = T))
muts = muts_dt[Jbrowse.quality == 'Y', mut_ID] %>% unlist()
paste('Number of mutations that passed QC:', length(muts)) # 255 

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

col_tumour_PD62341.2 = "#e56161"
col_tumour_PD63383.2 = "#bd0f71"

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
pdf('Results/20241119_p4_heatmap_all_mutations_val255.pdf')
pheatmap(mut_all_qc,
         cellwidth=10, cellheight=1.5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations, QC validated (255)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

pdf('Results/20241119_p4_heatmap_all_mutations_val255_large.pdf', height = 45)
pheatmap(mut_all_qc,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations, QC validated (255)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Trinucleotide context plots for mutations which are in the final 255 QC-validated set

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
  labs(x = 'Context', y = 'Count', title = 'Final set of mutations (255)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", linewidth = 0.1)
ggsave('Results/20241119_p4_mut_trins_255muts_pass.pdf', width = 8, height = 2.5)

# show mutations which were excluded by removing germline mutations on copy number change-regions and artefacts
muts1134_dt = data.table(read.table('Data/mutations_include_20241114_1134.txt'))
muts1134 = muts1134_dt[,V1] %>% unlist()
muts599_dt = data.table(read.table('Data/mutations_include_20241114_599.txt'))
muts599 = muts599_dt[,V1] %>% unlist()

mybed = twins_dt[mut_ID %in% setdiff(muts1134, muts599), c('Chrom', 'Pos', 'Ref', 'Alt')]
trins = get_trinucs(mybed, BSgenome.Hsapiens.UCSC.hg38)
dt = twins_dt[mut_ID %in% setdiff(muts1134, muts599)]
dt$trins=trins

mut_sign_counts = data.table(table(dt[, trins]))
setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'n = 535')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", linewidth = 0.1)
ggsave('Results/20241119_p4_mut_trins_excluded_535_germline.pdf', width = 8, height = 2.5)

mybed = twins_dt[mut_ID %in% setdiff(muts599, muts), c('Chrom', 'Pos', 'Ref', 'Alt')]
trins = get_trinucs(mybed, BSgenome.Hsapiens.UCSC.hg38)
dt = twins_dt[mut_ID %in% setdiff(muts599, muts)]
dt$trins=trins

mut_sign_counts = data.table(table(dt[, trins]))
setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'n = 344')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", linewidth = 0.1)
ggsave('Results/20241119_p4_mut_trins_excluded_344_qc.pdf', width = 8, height = 2.5)

######################################################################################################
# Identify classes of mutations of interest 

# Mutations present in all normal samples
muts_all_normal = twins_filtered_dt[sum_normal_mtr_vaf>=9, mut_ID] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_all_normal, c('mut_ID', samples_vaf, 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf', 'vaf_all_tumour'), with=FALSE]
paste('Number of mutations present in all normal samples:', length(muts_all_normal)) # 1

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

# Have all mutations been assigned to exactly one group?
print(length(c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal,
               muts_one_sample_normal, muts_one_sample_tumour, muts_tumour_PD62341, muts_tumour_PD63383, muts_tumour)))

# Manually assign remaining mutations 
# "chr3_195014637_A_G" "chr4_15905566_C_T"  "chr5_44907911_G_A"  "chr8_137568368_T_A"
muts_tumour = c(muts_tumour, 'chr5_44907911_G_A', 'chr8_137568368_T_A') # not assigned to tumour because also present in PD62341q (VAF = 0.11-0.12); however absent from n, v, ad therefore likely contamination, especially that VAF in some tumour samples is ~0.5
muts_PD63383_normal = c(muts_PD63383_normal, 'chr4_15905566_C_T') # present in some PD63383 samples and PD62341v
# chr4_15905566_C_T this one could be shared, I need to double check with Henry 

# Create and save a dataframe with all mut IDs used to reconstruct phylogeny (255) and assignment to the group
muts_classes = data.table(muts)
setnames(muts_classes, 'muts', 'mut_ID')
muts_classes[, mut_class := factor(fcase(
  mut_ID %in% muts_all_normal, 'shared_early_embryonic',
  mut_ID %in% muts_PD62341_normal, 'PD62341_specific_embryonic',
  mut_ID %in% muts_PD63383_normal, 'PD63383_specific_embryonic',
  mut_ID %in% muts_tumour, 'tumour',
  mut_ID %in% muts_tumour_PD62341, 'tumour_PD62341_only',
  mut_ID %in% muts_tumour_PD63383, 'tumour_PD63383_only',
  mut_ID %in% muts_one_sample_normal, '1_normal_sample',
  mut_ID %in% muts_one_sample_tumour, '1_tumour_sample',
  !mut_ID %in% c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal, muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383, muts_one_sample_normal, muts_one_sample_tumour), 'unassigned'
))]

write.table(muts_classes, 'Data/20241205_muts_classes_255.csv', sep = ',', quote=F, col.names=T, row.names=F)

######################################################################################################
# Analysis of early embryonic mutations 

# Heatmap to visualize each class of mutations of interest 

# Mutations present in all normal samples (only show normal samples)
mut_ee = as.matrix(twins_filtered_dt[mut_ID %in% c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal), c(samples_normal_vaf), with=FALSE])
rownames(mut_ee) = twins_filtered_dt[mut_ID %in% c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal),1] %>% unlist()  
colnames(mut_ee) = tstrsplit(colnames(mut_ee), '_VAF', fixed=TRUE, keep=1) %>% unlist()
col_annotation_normal = data.frame('TCF' = c(0, 0.1, 0.2, 0, 0, 0, 0, 0, 0, 0.3, 0.4, 0), # fraction of tumour cells
                                   Twin = c(rep('PD62341', 4), rep('PD63383', 6), rep('PD62341', 2)))
rownames(col_annotation_normal) = colnames(mut_ee)
annotation_colors_normal = list('TCF' = colorRampPalette(c('#f2e7e7', 'darkred'))(3),
                                Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))
pdf('Results/20241119_p4_heatmap_muts_earlyDev.pdf')
pheatmap(mut_ee,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation_normal,
         annotation_colors = annotation_colors_normal,
         main="Early embryonic mutations", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# The same heatmap, but without rownames (mutation IDs)
pdf('Results/20241119_p4_heatmap_muts_earlyDev_noRowNames.pdf')
pheatmap(mut_ee,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation_normal,
         annotation_colors = annotation_colors_normal,
         main="Early embryonic mutations", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()


# Plot VAF values for each mutation 
twins_vaf_melt = data.table::melt(twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE], id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
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
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, tissue := factor(tissue, levels = c('heart', 'spleen', 'cerebellum', 'pancreas', 'liver', 'skin'))]

for (mut in c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal)){
  
  # VAF of the mutation in each sample (tumour + normal)
  ggplot(twins_vaf_melt[mut_ID==mut], aes(x = sample, y = value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'Sample', y = 'VAF', col = 'Sample')+
    ylim(c(0, 0.8))+
    ggtitle(glue('{mut}'))+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_allNormal_dist_samples_{mut}.pdf'), height = 3, width = 6.5)
  
  # VAF of the mutation across normal samples
  ggplot(twins_vaf_melt[mut_ID==mut&status=='normal'], aes(x = tissue, y = value, col = twin))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'Tissue', y = 'VAF', col = 'Twin')+
    ggtitle(glue('{mut}'))+
    scale_color_manual(values = c(col_PD62341, col_PD63383))+
    ylim(c(0, 0.8))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_allNormal_dist_tissuse_{mut}.pdf'), height = 3, width = 6.5)
}

# Plot across tissues clean samples (not contaminated by cells from the other twin: exclude PD63383bb / skin and PD62341v / spleen)
twins_vaf_clean_melt = data.table::melt(twins_filtered_dt[, c('mut_ID', paste0(c('PD62341ad', 'PD62341h', 'PD62341n', 'PD62341q',
  'PD63383ak', 'PD63383t', 'PD63383ae', 'PD63383u'), '_VAF')), with=FALSE], id.vars = 'mut_ID')
twins_vaf_clean_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_vaf_clean_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341h', 'PD62341n', 'PD62341q',
  'PD63383ak', 'PD63383t', 'PD63383ae', 'PD63383u'))]
twins_vaf_clean_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_clean_melt[, tissue := factor(fcase(
  sample %in% c('PD62341ad', 'PD63383ak'), 'cerebellum',
  sample %in% c('PD62341h', 'PD63383t'), 'liver',
  sample %in% c('PD62341q', 'PD63383u'), 'pancreas',
  sample %in% c('PD62341n', 'PD63383ae'), 'heart', 
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_clean_melt[, tissue := factor(tissue, levels = c('heart', 'cerebellum', 'pancreas', 'liver'))]

for (mut in c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal)){
  
  # plot VAF of the mutation in each sample (tumour + clean normal)
  ggplot(twins_vaf_clean_melt[mut_ID==mut], aes(x = sample, y = value, col = twin))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'Sample', y = 'VAF', col = 'Sample')+
    ylim(c(0, 0.8))+
    ggtitle(glue('{mut}'))+
    scale_color_manual(values = c(col_PD62341, col_PD63383))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_allNormal_dist_samples_clean_{mut}.pdf'), height = 3, width = 6.5)
  
  # compare VAF across tissues from both twins (only clean normal samples)
  ggplot(twins_vaf_clean_melt[mut_ID==mut], aes(x = tissue, y = value, col = twin))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'Tissue', y = 'VAF', col = 'Twin')+
    ggtitle(glue('{mut}'))+
    scale_color_manual(values = c(col_PD62341, col_PD63383))+
    ylim(c(0, 0.8))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241119_p4_allNormal_dist_tissuse_clean_{mut}.pdf'), height = 3, width = 6.5)
  
}

# Show VAF across aggregated samples
agg_vaf_dt_PD62341 = twins_filtered_dt[, c('mut_ID', 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD62341_clean_vafLowerCI', 'agg_normal_PD62341_clean_vafUpperCI'), with=FALSE]
agg_vaf_dt_PD62341[,twin := 'PD62341']
setnames(agg_vaf_dt_PD62341, c('agg_normal_PD62341_clean_vaf', 'agg_normal_PD62341_clean_vafLowerCI','agg_normal_PD62341_clean_vafUpperCI'), c('vaf', 'lowerCI', 'upperCI'))
agg_vaf_dt_PD63383 = twins_filtered_dt[, c('mut_ID', 'agg_normal_PD63383_clean_vaf', 'agg_normal_PD63383_clean_vafLowerCI', 'agg_normal_PD63383_clean_vafUpperCI'), with=FALSE]
agg_vaf_dt_PD63383[,twin := 'PD63383']
setnames(agg_vaf_dt_PD63383, c('agg_normal_PD63383_clean_vaf', 'agg_normal_PD63383_clean_vafLowerCI','agg_normal_PD63383_clean_vafUpperCI'), c('vaf', 'lowerCI', 'upperCI'))
agg_vaf_dt = rbind(agg_vaf_dt_PD62341, agg_vaf_dt_PD63383)
ggplot(agg_vaf_dt[mut_ID %in% muts_all_normal], aes(x = twin, y = vaf, color = twin))+
  geom_point(size = 1.5) + 
  theme_classic(base_size = 14)+
  labs(x = 'Twin', y = 'VAF', col = 'Twin')+
  ggtitle(glue('chr1_38827952_C_A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  ylim(c(0, 0.8))+
  theme(legend.position = 'None')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, width = 0),
                position = position_dodge(width = 0.90))
ggsave('Results/20241119_aggvaf_withCI_earlyMuts.pdf', height = 3, width = 3)

ggplot(agg_vaf_dt[mut_ID %in% muts_PD62341_normal], aes(x = twin, y = vaf, color = twin, group = twin))+
  geom_point(size = 2) + 
  facet_wrap(~mut_ID, nrow=2) + 
  theme_classic(base_size = 12)+
  labs(x = 'Twin', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  ylim(c(0, 0.8))+
  theme(legend.position = 'None')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, width = 0),
                position = position_dodge(width = 0.90))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341.pdf', height = 4, width = 8)

ggplot(agg_vaf_dt[mut_ID %in% muts_PD63383_normal], aes(x = twin, y = vaf, color = twin, group = twin))+
  geom_point(size = 2) + 
  facet_wrap(~mut_ID, nrow=2) + 
  theme_classic(base_size = 12)+
  labs(x = 'Twin', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  ylim(c(0, 0.8))+
  theme(legend.position = 'None')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, width = 0),
                position = position_dodge(width = 0.90))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383.pdf', height = 4, width = 10)

# Aggregate, mutation on y axis
agg_vaf_dt[, mut_ID := factor(mut_ID, levels = {
  agg_vaf_dt[, .(mean_col = mean(vaf)), by = mut_ID][order(mean_col), mut_ID]})]

ggplot(agg_vaf_dt[mut_ID %in% muts_PD62341_normal], aes(x = vaf, y = mut_ID, color = twin, group = twin))+
  geom_point(size = 2) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  xlim(c(0, 0.8))+
  geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, width = 0))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_aggCI.pdf', height = 2.5, width = 5)

ggplot(agg_vaf_dt[mut_ID %in% muts_PD63383_normal], aes(x = vaf, y = mut_ID, color = twin, group = twin))+
  geom_point(size = 2) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  xlim(c(0, 0.8))+
  geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, width = 0))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_aggCI.pdf', height = 3.5, width = 5)

ggplot(agg_vaf_dt[mut_ID %in% muts_PD63383_normal], aes(x = vaf, y = mut_ID, color = twin, group = twin))+
  geom_point(size = 2) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  xlim(c(0, 0.4))+
  geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, width = 0))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_aggCI_xlim2.pdf', height = 3.5, width = 5)

# For mutations enriched in one twin, do a scatterplot
twins_vaf_PD62341_melt = melt(twins_filtered_dt[, c('mut_ID', samples_normal_PD62341_clean_vaf), with=FALSE], id.vars = 'mut_ID')
twins_vaf_PD63383_melt = melt(twins_filtered_dt[, c('mut_ID', samples_normal_PD63383_clean_vaf), with=FALSE], id.vars = 'mut_ID')
setnames(twins_vaf_PD62341_melt, 'value', 'vaf_PD62341')
setnames(twins_vaf_PD63383_melt, 'value', 'vaf_PD63383')
twins_vafs_merged = cbind(twins_vaf_PD62341_melt[,c('mut_ID', 'vaf_PD62341'), with=F], 
                          twins_vaf_PD63383_melt[,'vaf_PD63383', with=F])

ggplot(twins_vafs_merged[mut_ID %in% muts_PD62341_normal], aes(x = vaf_PD62341, y = vaf_PD63383))+
  geom_point(size = 2.5) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF in PD62341', y = 'VAF in PD63383', col = 'Twin')+
  ggtitle(glue('PD62341'))+
  ylim(c(0, 0.7))+
  xlim(c(0, 0.7))+
  coord_equal(ratio = 1)+
  theme(legend.position = 'None')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_scatter.pdf', height = 3, width = 3)

ggplot(twins_vafs_merged[mut_ID %in% muts_PD63383_normal], aes(x = vaf_PD63383, y = vaf_PD62341))+
  geom_point(size = 2.5) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF in PD63383', y = 'VAF in PD62341', col = 'Twin')+
  ggtitle(glue('PD63383'))+
  ylim(c(0, 0.4))+
  xlim(c(0, 0.4))+
  coord_equal(ratio = 1)+
  theme(legend.position = 'None')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_scatter.pdf', height = 3, width = 3)

# Alternative way of displaying mutations
twins_vaf_melt = data.table::melt(twins_filtered_dt[, c('mut_ID', samples_normal_PD62341_clean_vaf, samples_normal_PD63383_clean_vaf), with=FALSE], id.vars = 'mut_ID')
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
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, mut_ID := factor(mut_ID, levels = {
  twins_vaf_melt[, .(mean_col = mean(value)), by = mut_ID][order(mean_col), mut_ID]
})]

# Mutation on the y axis 
ggplot(twins_vaf_melt[mut_ID %in% muts_PD62341_normal&status=='normal'], aes(x = value, y = mut_ID, color = twin, group = twin))+
  geom_point(size = 2.5, alpha = 0.6) +
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  xlim(c(0, 0.8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_dotplot.pdf', height = 3, width = 6)

ggplot(twins_vaf_melt[mut_ID %in% muts_PD63383_normal&status=='normal'], aes(x = value, y = mut_ID, color = twin, group = twin))+
  geom_point(size = 2.5, alpha = 0.6) +
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  xlim(c(0, 0.8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_dotplot.pdf', height = 3.5, width = 6)

ggplot(twins_vaf_melt[mut_ID %in% muts_PD63383_normal&status=='normal'], aes(x = value, y = mut_ID, color = twin, group = twin))+
  geom_point(size = 2.5, alpha = 0.6) +
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  xlim(c(0, 0.4))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_dotplot_xlim2.pdf', height = 3.5, width = 6)

# Mutation on the x axis
ggplot(twins_vaf_melt[mut_ID %in% muts_PD62341_normal&status=='normal'], aes(x = mut_ID, y = value, color = twin, group = twin))+
  geom_point(size = 2.5, alpha = 0.6) +
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.8))+
  geom_hline(yintercept = 0.4, colour="black", size = 0.2)+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 5))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_dotplot_xaxis.pdf', height = 3, width = 6.5)

ggplot(twins_vaf_melt[mut_ID %in% muts_PD63383_normal&status=='normal'], aes(x = mut_ID, y = value, color = twin, group = twin))+
  geom_point(size = 2.5, alpha = 0.6) +
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.4))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 5))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_dotplot_xaxis.pdf', height = 3, width = 6.5)

# Plots to compare VAFs across tissues (can this tell us anything about the contribution to different tissues later on?)
ggplot(twins_vaf_melt[mut_ID %in% muts_PD63383_normal&status=='normal'&twin=='PD63383'], aes(x = tissue, y = value, color = mut_ID, group = mut_ID))+
  geom_point(size = 2.5, alpha = 0.6) +
  geom_line()+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'VAF', col = 'Mutation')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.4))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_acrossTissues.pdf', height = 3.5, width = 6.5)

twins_vaf_melt[, tissue_twin := paste0(twin, '_', tissue)]
ggplot(twins_vaf_melt[mut_ID %in% muts_PD63383_normal&status=='normal'], aes(x = tissue_twin, y = value, color = mut_ID, group = mut_ID))+
  geom_point(size = 2.5, alpha = 0.6) +
  geom_line()+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'VAF', col = 'Mutation')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.4))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_acrossTissues_bothTwins.pdf', height = 3.6, width = 6.5)

ggplot(twins_vaf_melt[mut_ID %in% muts_PD62341_normal&status=='normal'&twin=='PD62341'], aes(x = tissue, y = value, color = mut_ID, group = mut_ID))+
  geom_point(size = 2.5, alpha = 0.6) +
  geom_line()+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'VAF', col = 'Mutation')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_acrossTissues.pdf', height = 3.5, width = 6.5)

ggplot(twins_vaf_melt[mut_ID %in% c(muts_PD63383_normal, 'chr1_38827952_C_A') &status=='normal'&twin=='PD63383'], aes(x = tissue, y = value, color = mut_ID, group = mut_ID))+
  geom_point(size = 2.5, alpha = 0.6) +
  geom_line()+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'VAF', col = 'Mutation')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_acrossTissues_pluschr1.pdf', height = 3.5, width = 6.5)

ggplot(twins_vaf_melt[mut_ID %in% c(muts_PD62341_normal, 'chr1_38827952_C_A') &status=='normal'&twin=='PD62341'], aes(x = tissue, y = value, color = mut_ID, group = mut_ID))+
  geom_point(size = 2.5, alpha = 0.6) +
  geom_line()+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'VAF', col = 'Mutation')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_acrossTissues_pluschr1.pdf', height = 3.5, width = 6.5)

ggplot(twins_vaf_melt[mut_ID %in% c(muts_PD63383_normal) &status=='normal'&twin=='PD62341'], aes(x = tissue, y = value, color = mut_ID, group = mut_ID))+
  geom_point(size = 2.5, alpha = 0.6) +
  geom_line()+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'VAF', col = 'Mutation')+
  ggtitle(glue('PD62341 samples'))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.3))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_acrossTissues_inPD62341.pdf', height = 3.5, width = 6.5)

for (mut in c(muts_PD63383_normal)){
  ggplot(twins_vaf_melt[mut_ID ==mut &status=='normal'&twin=='PD63383'], aes(x = tissue, y = value, group=mut_ID))+
    geom_point(size = 2.5, alpha = 0.6) +
    geom_line()+
    theme_classic(base_size = 14)+
    labs(x = 'Tissue', y = 'VAF')+
    ggtitle(glue('{mut}'))+
    guides(color = guide_legend(override.aes = list(alpha = 1)))+
    ylim(c(0, 0.4))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
  ggsave(glue('Results/20241119_aggvaf_withCI_earlyMutsPD63383_acrossTissues_{mut}.pdf'), height = 3.5, width = 6.5)
  
}

ggplot(twins_vaf_melt[mut_ID =='chr1_38827952_C_A' &status=='normal'&twin=='PD63383'], aes(x = tissue, y = value, group=mut_ID))+
  geom_point(size = 2.5, alpha = 0.6) +
  geom_line()+
  theme_classic(base_size = 14)+
  labs(x = 'Tissue', y = 'VAF')+
  ggtitle(glue('chr1_38827952_C_A'))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ylim(c(0, 0.8))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
ggsave(glue('Results/20241119_aggvaf_withCI_earlyMutsPD63383_acrossTissues_chr1_388.pdf'), height = 3.5, width = 6.5)


for (mut in c(muts_PD62341_normal, 'chr1_38827952_C_A')){
  ggplot(twins_vaf_melt[mut_ID ==mut &status=='normal'&twin=='PD62341'], aes(x = tissue, y = value, group=mut_ID))+
    geom_point(size = 2.5, alpha = 0.6) +
    geom_line()+
    theme_classic(base_size = 14)+
    labs(x = 'Tissue', y = 'VAF')+
    ggtitle(glue('{mut}'))+
    guides(color = guide_legend(override.aes = list(alpha = 1)))+
    ylim(c(0, 0.8))+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size = 8))
  ggsave(glue('Results/20241119_aggvaf_withCI_earlyMutsPD62341_acrossTissues_{mut}.pdf'), height = 3.5, width = 6.5)
  
}

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

# Are PD63383-specific mutations really absent from PD62341?
for (mut in muts_PD63383_normal){
  ggplot(agg_vaf_dt[mut_ID == mut], aes(x = twin, y = vaf, color = twin))+
    geom_point(size = 2.5) + 
    theme_classic(base_size = 12)+
    labs(x = 'Twin', y = 'VAF', col = 'Twin')+
    ggtitle(glue('{mut}'))+
    scale_color_manual(values = c(col_PD62341, col_PD63383))+
    ylim(c(0, 0.3))+
    theme(legend.position = 'None')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    geom_errorbar(aes(ymin = lowerCI, ymax = upperCI, width = 0),
                  position = position_dodge(width = 0.90))
  ggsave(glue('Results/20241119_aggvaf_withCI_earlyMutsPD63383_{mut}.pdf'), height = 3, width = 3)
}

######################################################################################################
# Tumour evolution: origin 

# Q1: which twin did the tumour originate from?

# Heatmap: early embryonic mutations across normal and tumour samples
# keep order the same as on the previous heatmap 
mut_order = c('chr1_38827952_C_A', 
              'chr3_77633967_C_T', 'chr4_74625500_G_T', 'chr13_50815806_A_G', 'chr7_149688370_C_T',
              'chr1_103587565_A_C', 'chr21_40193588_G_A', 'chr7_73831920_C_T', 'chr4_75704880_G_A',
              'chr6_165179306_G_A', 'chr11_34011887_C_T', 'chrX_115066661_C_T',
              'chr3_62055057_C_G', 'chr3_62055077_G_C', 'chr16_5479739_C_T', 'chr3_50106043_C_T', 
              'chr14_105458006_C_A', 'chr17_33422229_C_A', 'chr2_95662131_G_A')
twins_filtered_dt_ee = twins_filtered_dt[mut_ID %in% c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal), c('mut_ID', samples_vaf), with=F]
twins_filtered_dt_ee = twins_filtered_dt_ee[, mut_ID := factor(mut_ID, levels = c(mut_order))]
twins_filtered_dt_ee = twins_filtered_dt_ee[order(mut_ID)]

# keep order of columns the same as in previous heatmaps
setcolorder(twins_filtered_dt_ee, c('mut_ID', paste0(c('PD63383aq', 'PD63383ap', 'PD62341ae', 'PD62341am', 'PD62341ap',
                                                       'PD62341b', 'PD62341aj', 'PD62341u', 'PD62341ak', 'PD62341ag',
                                                       'PD62341n', 'PD62341h', 'PD62341q', 'PD62341aa', 'PD62341ad', 'PD62341v',
                                                       'PD63383bb', 'PD63383t', 'PD63383ae', 'PD63383u', 'PD63383w', 'PD63383ak'), '_VAF')))

mut_eet = as.matrix(twins_filtered_dt_ee[mut_ID %in% c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal), 
                                         paste0(c('PD63383aq', 'PD63383ap', 'PD62341ae', 'PD62341am', 'PD62341ap','PD62341b', 'PD62341aj', 'PD62341u', 'PD62341ak', 'PD62341ag',
                                                  'PD62341n', 'PD62341h', 'PD62341q', 'PD62341aa', 'PD62341ad', 'PD62341v','PD63383bb', 'PD63383t', 'PD63383ae', 'PD63383u', 'PD63383w', 'PD63383ak'), '_VAF'), with=FALSE])
rownames(mut_eet) = twins_filtered_dt_ee[mut_ID %in% c(muts_all_normal, muts_PD62341_normal, muts_PD63383_normal),1] %>% unlist()  
colnames(mut_eet) = tstrsplit(colnames(mut_eet), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation_tumour = data.frame('TCF' = c(0.9, 0.9, 0.9, 0.8, 0.8, 0.7, 0.6, 0.5, 0.5, 0.3, 
                                              0, 0.5, 0.1, 0.2, 0, 0, 0.3, 0, 0, 0, 0, 0), # fraction of tumour cells
                            Status = c(rep('tumour', 10), rep('normal', 12)),
                            Twin = c(rep('PD63383', 2),  rep('PD62341', 14), rep('PD63383', 6)))
rownames(col_annotation_tumour) = colnames(mut_eet)
annotation_colors_tumour = list('TCF' = colorRampPalette(c('#f2e7e7', 'darkred'))(11),
                         Status = c(normal=col_normal, tumour=col_tumour),
                         Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

pdf('Results/20241119_p4_heatmap_muts_earlyDev_withTumour.pdf')
pheatmap(mut_eet,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation_tumour,
         annotation_colors = annotation_colors_tumour,
         main="Early embryonic mutations", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = F, cluster_cols = F, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Show VAF across aggregated samples, adding tumour samples 
agg_vaf_dt_PD62341 = twins_filtered_dt[, c('mut_ID', 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD62341_clean_vafLowerCI', 'agg_normal_PD62341_clean_vafUpperCI'), with=FALSE]
agg_vaf_dt_PD62341[,twin := 'PD62341']
agg_vaf_dt_PD62341[,status := 'normal']
setnames(agg_vaf_dt_PD62341, c('agg_normal_PD62341_clean_vaf', 'agg_normal_PD62341_clean_vafLowerCI','agg_normal_PD62341_clean_vafUpperCI'), c('vaf', 'lowerCI', 'upperCI'))
agg_vaf_dt_PD63383 = twins_filtered_dt[, c('mut_ID', 'agg_normal_PD63383_clean_vaf', 'agg_normal_PD63383_clean_vafLowerCI', 'agg_normal_PD63383_clean_vafUpperCI'), with=FALSE]
agg_vaf_dt_PD63383[,twin := 'PD63383']
agg_vaf_dt_PD63383[,status := 'normal']
setnames(agg_vaf_dt_PD63383, c('agg_normal_PD63383_clean_vaf', 'agg_normal_PD63383_clean_vafLowerCI','agg_normal_PD63383_clean_vafUpperCI'), c('vaf', 'lowerCI', 'upperCI'))
agg_vaf_tumour_PD62341 = twins_filtered_dt[, c('mut_ID', 'vaf_all_tumour_PD62341', 'vaf_all_tumour_PD62341_lowerCI', 'vaf_all_tumour_PD62341_upperCI'), with=FALSE]
agg_vaf_tumour_PD62341[,twin := 'PD62341']
agg_vaf_tumour_PD62341[,status := 'tumour']
setnames(agg_vaf_tumour_PD62341, c('vaf_all_tumour_PD62341', 'vaf_all_tumour_PD62341_lowerCI', 'vaf_all_tumour_PD62341_upperCI'), c('vaf', 'lowerCI', 'upperCI'))
agg_vaf_tumour_PD63383 = twins_filtered_dt[, c('mut_ID', 'vaf_all_tumour_PD63383', 'vaf_all_tumour_PD63383_lowerCI', 'vaf_all_tumour_PD63383_upperCI'), with=FALSE]
agg_vaf_tumour_PD63383[,twin := 'PD63383']
agg_vaf_tumour_PD63383[,status := 'tumour']
setnames(agg_vaf_tumour_PD63383, c('vaf_all_tumour_PD63383', 'vaf_all_tumour_PD63383_lowerCI', 'vaf_all_tumour_PD63383_upperCI'), c('vaf', 'lowerCI', 'upperCI'))
agg_vaf_tumour = twins_filtered_dt[, c('mut_ID', 'vaf_all_tumour', 'vaf_all_tumour_lowerCI', 'vaf_all_tumour_upperCI'), with=FALSE]
agg_vaf_tumour[,twin := 'both']
agg_vaf_tumour[,status := 'tumour']
setnames(agg_vaf_tumour, c('vaf_all_tumour', 'vaf_all_tumour_lowerCI', 'vaf_all_tumour_upperCI'), c('vaf', 'lowerCI', 'upperCI'))

agg_vaf_dt_tumour = rbind(agg_vaf_dt_PD62341, agg_vaf_dt_PD63383, agg_vaf_tumour_PD62341, agg_vaf_tumour_PD63383)
agg_vaf_dt_tumour[, sample := paste0(twin, ', ', status)]
agg_vaf_dt_tumour[, sample_type := factor(fcase(
  status =='normal', paste0(twin, ', ', status),
  status == 'tumour', 'tumour'))]
agg_vaf_dt_tumour[, mut_ID := factor(mut_ID, levels = {
  agg_vaf_dt[, .(mean_col = mean(vaf)), by = mut_ID][order(mean_col), mut_ID]})]

ggplot(agg_vaf_dt_tumour[mut_ID %in% muts_PD62341_normal], aes(x = vaf, y = mut_ID, color = sample, group = sample))+
  geom_point(size = 2) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  scale_color_manual(values = c(col_PD62341, col_tumour_PD62341, col_PD63383, col_tumour_PD63383))+
  xlim(c(0, 0.8))+
  geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, width = 0))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_aggCI_xlim2.pdf', height = 2.5, width = 5)

ggplot(agg_vaf_dt_tumour[mut_ID %in% muts_PD63383_normal], aes(x = vaf, y = mut_ID, color = sample, group = sample))+
  geom_point(size = 2) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_tumour_PD62341, col_PD63383, col_tumour_PD63383))+
  xlim(c(0, 0.8))+
  geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, width = 0))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_aggCI_xlim2.pdf', height = 3.5, width = 5)

agg_vaf_dt_tumour2 = rbind(agg_vaf_dt_PD62341, agg_vaf_dt_PD63383, agg_vaf_tumour)
agg_vaf_dt_tumour2[, sample_type := factor(fcase(
  status =='normal', paste0(twin, ', ', status),
  status == 'tumour', 'tumour'))]
agg_vaf_dt_tumour2[, mut_ID := factor(mut_ID, levels = {
  agg_vaf_dt[, .(mean_col = mean(vaf)), by = mut_ID][order(mean_col), mut_ID]})]

ggplot(agg_vaf_dt_tumour2[mut_ID %in% muts_PD62341_normal], aes(x = vaf, y = mut_ID, color = sample_type, group = sample_type))+
  geom_point(size = 2) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample')+
  ggtitle(glue('Mutations enriched in PD62341'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  xlim(c(0, 0.8))+
  geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, width = 0))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD62341_aggCI_xlim2_sampleType.pdf', height = 2.5, width = 5)

ggplot(agg_vaf_dt_tumour2[mut_ID %in% muts_PD63383_normal], aes(x = vaf, y = mut_ID, color = sample_type, group = sample_type))+
  geom_point(size = 2) + 
  theme_classic(base_size = 12)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample')+
  ggtitle(glue('Mutations enriched in PD63383'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  xlim(c(0, 0.8))+
  geom_errorbar(aes(xmin = lowerCI, xmax = upperCI, width = 0))
ggsave('Results/20241119_aggvaf_withCI_earlyMutsPD63383_aggCI_xlim2_sampleType.pdf', height = 3.5, width = 5)

######################################################################################################
# Tumour evolution: subclonal structure 

# First, plot a heatmap to show tumour-specific mutations 
mut_tumour = as.matrix(twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383), c(samples_vaf), with=FALSE])
rownames(mut_tumour) = twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383),1] %>% unlist()  
colnames(mut_tumour) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_alltumour.pdf')
pheatmap(mut_tumour,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation_tumour,
         annotation_colors = annotation_colors_tumour,
         main="Mutations specific to the tumour", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Plot coverage vs VAF to identify possible subclones 
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

# Mutations present in one sample only 
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

# Can PCR duplicates specifically on the mutant strand in chr1_388 explain the inflated VAF?
qbinom(0.05, 174, 104/174)/174
# this is agg from 5 samples
# say each sample has 2 duplicate reads 
qbinom(0.05, 164, 94/164)/164
# okay this does actually get us to a reasonable number 
# but you would require > 2 duplicates per sample and your VAF is still really high  

qbinom(0.05, 199, 56/199)/199
qbinom(0.05, 189, 46/189)/189
# this would be 0.15 so still too high 
# again, you would have to have > 2 duplicates per sample for this to be viable




