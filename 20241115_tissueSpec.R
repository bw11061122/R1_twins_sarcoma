###################################################################################################################################
# SUPPLEMENTARY SCRIPT 

# Script to analyse whether there is any tissue specificity
# November 2024 
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# Cleaned up script to identify mutations of interest 
# Script to look at tumour-specific mutations and tumour evolution 

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
library(viridis)
library(grid)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

###################################################################################################################################
# INPUT FILES 

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = fread('Data/pileup_merged_20241016.tsv') # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Out/F1/F1_mutations_final_20241208_255.txt') 
muts = muts$V1 %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 255

# Create dataframe with mutations that passed current filters 
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Import dataframe with purity estimates
purity_dt = fread('Out/F2/F2_estimates_tumour_cont_26muts_median_20241208.csv')

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
twins_filtered_dt[, sum_tumour_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, sum_normal_mtr := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_dt[, sum_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_dt[, sum_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_dt[, sum_tumour_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, sum_tumour_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, sum_normal_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, sum_normal_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# Add presence / absence based on VAF data
twins_filtered_dt[, sum_tumour_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_filtered_dt[, sum_normal_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_vaf]
twins_filtered_dt[, sum_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD62341_vaf]
twins_filtered_dt[, sum_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD63383_vaf]
twins_filtered_dt[, sum_tumour_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_filtered_dt[, sum_tumour_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]
twins_filtered_dt[, sum_normal_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_filtered_dt[, sum_normal_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

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

######################################################################################################
# ARE THERE MUTATIONS PRESENT IN A SPECIFIC TISSUE OF BOTH TWINS?

# search through samples which are from the same tissue 
samples_spleen = c('PD62341v_VAF', 'PD63383w_VAF') 
samples_skin = c('PD62341aa_VAF', 'PD63383bb_VAF') 
samples_cerebellum = c('PD62341ad_VAF', 'PD63383ak_VAF') 
samples_liver = c('PD62341h_VAF', 'PD63383t_VAF') 
samples_heart = c('PD62341n_VAF', 'PD63383ae_VAF') 
samples_pancreas = c('PD62341q_VAF', 'PD63383u_VAF') 

twins_filtered_dt[, mean_vaf_spleen := apply(.SD, 1, mean), .SDcols = samples_spleen]
twins_filtered_dt[, mean_vaf_skin := apply(.SD, 1, mean), .SDcols = samples_skin]
twins_filtered_dt[, mean_vaf_cerebellum := apply(.SD, 1, mean), .SDcols = samples_cerebellum]
twins_filtered_dt[, mean_vaf_liver := apply(.SD, 1, mean), .SDcols = samples_liver]
twins_filtered_dt[, mean_vaf_heart := apply(.SD, 1, mean), .SDcols = samples_heart]
twins_filtered_dt[, mean_vaf_pancreas := apply(.SD, 1, mean), .SDcols = samples_pancreas]
twins_filtered_dt[, mean_vaf_normal := apply(.SD, 1, mean), .SDcols = samples_normal_vaf]

# you can also specify sample order and plot them as a heatmap - that will be easy to see 
mut_tissues_nt = twins_filtered_dt[, c('mut_ID', samples_normal_vaf), with=FALSE] # only include mutations that passed your filter 
setcolorder(mut_tissues_nt, c('mut_ID', 'PD62341ad_VAF', 'PD63383ak_VAF', 
                              'PD62341aa_VAF', 'PD63383bb_VAF',
                              'PD62341q_VAF', 'PD63383u_VAF',
                              'PD62341h_VAF', 'PD63383t_VAF',
                              'PD62341n_VAF', 'PD63383ae_VAF',
                              'PD62341v_VAF', 'PD63383w_VAF'))
mut_tissues_nt_mat = as.matrix(mut_tissues_nt[, c(samples_normal_vaf), with=FALSE])
rownames(mut_tissues_nt_mat) = mut_tissues_nt[,1] %>% unlist()  
colnames(mut_tissues_nt_mat) = tstrsplit(colnames(mut_tissues_nt_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Twin = rep(c('PD62341', 'PD63383'), 6), 
                            Tissue = rep(c('cerebellum', 'skin', 'pancreas', 'liver', 'heart', 'spleen'), each = 2))
rownames(col_annotation) = colnames(mut_tissues_nt_mat)
annotation_colors = list(Twin = c(PD62341=col_PD62341, PD63383=col_PD63383), 
                         Tissue = c(cerebellum = '#0957cc', skin = '#8accfc', 
                                    pancreas = '#bc7f0c', liver = '#e4bf28', 
                                    heart = '#d1140d', spleen = '#f1732a'))

# heatmap (plot, don't save - not as useful)
pheatmap(mut_tissues_nt_mat,
         cellwidth=40, cellheight=1,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="255 mutations", 
         legend = T, 
         treeheight_row = 0, 
         cluster_rows = T, cluster_cols = F, 
         show_rownames = F, show_colnames = T,
         fontsize=20, cexCol=2) 

# search for mutations which seem to be enriched in a given tissue over the others 
twins_filtered_dt[, ratio_spleen_normal := mean_vaf_spleen / mean_vaf_normal]
twins_filtered_dt[, ratio_skin_normal := mean_vaf_skin / mean_vaf_normal]
twins_filtered_dt[, ratio_cerebellum_normal := mean_vaf_cerebellum / mean_vaf_normal]
twins_filtered_dt[, ratio_liver_normal := mean_vaf_liver / mean_vaf_normal]
twins_filtered_dt[, ratio_heart_normal := mean_vaf_heart / mean_vaf_normal]
twins_filtered_dt[, ratio_pancreas_normal := mean_vaf_pancreas / mean_vaf_normal]

twins_filtered_dt[, ratio_spleen_normal_log := log2(ratio_spleen_normal)]
twins_filtered_dt[, ratio_skin_normal_log := log2(ratio_skin_normal)]
twins_filtered_dt[, ratio_cerebellum_normal_log := log2(ratio_cerebellum_normal)]
twins_filtered_dt[, ratio_liver_normal_log := log2(ratio_liver_normal)]
twins_filtered_dt[, ratio_heart_normal_log := log2(ratio_heart_normal)]
twins_filtered_dt[, ratio_pancreas_normal_log := log2(ratio_pancreas_normal)]

# can I identify any instances where the log is higher?
# spleen
twins_filtered_dt[ratio_spleen_normal_log > 1 & PD62341v_VAF > 0.1 & PD63383w_VAF > 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
twins_filtered_dt[ratio_spleen_normal_log < -1 & ratio_skin_normal != 0 & PD62341v_VAF < 0.1 & PD63383w_VAF < 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
# empty

# skin
twins_filtered_dt[ratio_skin_normal_log > 1 & PD62341aa_VAF > 0.1 & PD63383bb_VAF > 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
twins_filtered_dt[ratio_skin_normal_log < -1 & ratio_skin_normal != 0 & PD62341aa_VAF < 0.1 & PD63383bb_VAF < 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
# empty 

# cerebellum
twins_filtered_dt[ratio_cerebellum_normal_log > 1 & PD62341ad_VAF > 0.1 & PD63383ak_VAF > 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
twins_filtered_dt[ratio_cerebellum_normal_log < -1 & ratio_cerebellum_normal != 0 & PD62341ad_VAF < 0.1 & PD63383ak_VAF < 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
# empty

# liver 
twins_filtered_dt[ratio_liver_normal_log > 1 & PD62341h_VAF > 0.1 & PD63383t_VAF > 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
twins_filtered_dt[ratio_liver_normal_log < -1 & ratio_liver_normal != 0 & PD62341h_VAF < 0.1 & PD63383t_VAF < 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
# empty

# heart 
twins_filtered_dt[ratio_heart_normal_log > 1 & PD62341n_VAF > 0.1 & PD63383ae_VAF > 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
twins_filtered_dt[ratio_heart_normal_log < -1 & ratio_heart_normal != 0 & PD62341n_VAF < 0.1 & PD63383ae_VAF < 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
# empty 

# pancreas
twins_filtered_dt[ratio_pancreas_normal_log > 1 & PD62341q_VAF > 0.1 & PD63383u_VAF > 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
twins_filtered_dt[ratio_pancreas_normal_log < -1 & ratio_pancreas_normal != 0 & PD62341q_VAF < 0.1 & PD63383u_VAF < 0.1 & sum_normal_PD62341_vaf > 1 & sum_normal_PD63383_vaf > 1]
# empty 