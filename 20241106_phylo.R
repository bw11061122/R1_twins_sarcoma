###################################################################################################################################
# SCRIPT 4

# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-11-06
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
library(beeswarm)
library(viridis)
library(grid)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

###################################################################################################################################
# INPUT FILES 

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241114_599.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 599
twins_filtered_dt = twins_dt[mut_ID %in% muts]

muts_qc = read.table('Data/mutations_passedQuality_20241114.txt') %>% unlist()

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

######################################################################################################
# Identify mutations for which presence is not clear (MTR >= 4 but VAF < 0.1 or VAF > 0.1 but MTR < 4)

twins_filtered_dt[sum_all_mtr_vaf == sum_all_vaf & sum_all_mtr_vaf == sum_all_mtr] # 331 
paste('Number of mutations with consistent VAF and MTR:', dim(twins_filtered_dt[sum_all_mtr_vaf == sum_all_vaf & sum_all_mtr_vaf == sum_all_mtr])[1])

twins_filtered_dt[, diff_mtr_vaf := sum_all_vaf - sum_all_mtr]
table(twins_filtered_dt[, diff_mtr_vaf])

pdf('Results/20241114_p4_genotyping_diffVafMtr.pdf')
hist(twins_filtered_dt[, diff_mtr_vaf], breaks = 20, xlab = 'Difference in the # samples with mutation (VAF - MTR)',
     main = 'Genotyping: difference in # samples with mutation')
dev.off()

# are mutations with large difference more often assigned or not assigned?
muts_diff0 = twins_filtered_dt[diff_mtr_vaf == 0, mut_ID] %>% unlist()
muts_diff5 = twins_filtered_dt[abs(diff_mtr_vaf) >= 5, mut_ID] %>% unlist()

sum(muts_diff0 %in% muts_assigned) # 324 (out of 331 so not bad)
sum(muts_diff5 %in% muts_assigned) # 22 (very rare)
sum(muts_diff5 %in% muts_unassigned) # 53 (1/3 unassigned mutations!)

# Create a dataframe with 3 levels: present (VAF AND MTR), maybe (VAF OR MTR but not both), absent (neither VAF nor MTR)
twins_filtered_dt_binary = data.table(twins_filtered_dt[, c('mut_ID', samples_vaf, samples_mtr), with=FALSE])

for (sample in samples_names){
  sample_vaf = paste0(sample, '_VAF')
  sample_mtr = paste0(sample, '_MTR')
  twins_filtered_dt_binary[, glue('{sample}_binary') := fcase(
    get(sample_vaf) >= 0.1 & get(sample_mtr) >= 4, 1,
    get(sample_vaf) < 0.1 & get(sample_mtr) < 4, 0,
    get(sample_vaf) >= 0.1 & get(sample_mtr) < 4 | get(sample_vaf) < 0.1 & get(sample_mtr) >= 4, 0.5
  )]
}

samples_binary = paste0(samples_names, '_binary')
twins_filtered_dt_binary = twins_filtered_dt_binary[, c('mut_ID', samples_binary), with=FALSE]
mut_all_binary = as.matrix(twins_filtered_dt_binary[, c(samples_binary), with=FALSE])
rownames(mut_all_binary) = twins_filtered_dt_binary[,1] %>% unlist()  
colnames(mut_all_binary) = tstrsplit(colnames(mut_all_binary), '_binary', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)),
                            Purity = c(0, 0.1, 0.2, 0, 0.9, 0.3, 0.5, 0.6, 0.8, 0.8, 0, 0, 0, 0, 0, 0.9, 0.9, 0.3, 0.7, 0.4, 0, 0.5)) # fraction of tumour cells 
rownames(col_annotation) = colnames(mut_all_binary)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), 
                         Twin = c(PD62341=col_PD62341, PD63383=col_PD63383),
                         Purity = colorRampPalette(c('#f2e7e7', 'darkred'))(11))

# heatmap
pdf('Results/20241114_p4_heatmap_all_mutations599_binary.pdf')
pheatmap(mut_all_binary,
         cellwidth=10, cellheight=0.4,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations (599)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off() # actually looks worse than with actual VAF data 

######################################################################################################
# Identify groups of mutations of interest 

# MUTATIONS OF INTEREST:
# Present in all samples
# Present in only a single sample
# Present in all normal samples 
# Tumour specific  
# Absent from PD63383 normal samples (PD62341 subclonal)

# Create lists of mutations 
muts = twins_filtered_dt[, mut_ID] %>% unlist()
muts_all_samples = twins_filtered_dt[sum_tumour_mtr_vaf==10 & sum_normal_mtr_vaf==12, mut_ID] %>% unlist()

muts_all_normal = twins_filtered_dt[sum_normal_mtr_vaf==12, mut_ID] %>% unlist()
muts_all_normal_nt = setdiff(muts_all_normal, muts_all_samples)

muts_all_tumour = twins_filtered_dt[sum_tumour_mtr_vaf==10, mut_ID] %>% unlist()
muts_all_tumour_nt = setdiff(muts_all_tumour, muts_all_samples)

muts_only_normal = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_mtr_vaf>1, mut_ID] %>% unlist() # exclude mutations present in a single sample
muts_only_normal_PD62341 = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_PD62341_mtr_vaf>1 & sum_normal_PD63383_mtr_vaf==0, mut_ID] %>% unlist() # exclude mutations present in a single sample
muts_only_normal_PD63383 = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_PD62341_mtr_vaf==0 & sum_normal_PD63383_mtr_vaf>1, mut_ID] %>% unlist() # exclude mutations present in a single sample

muts_only_tumour = twins_filtered_dt[sum_normal_mtr_vaf==0 & sum_tumour_mtr_vaf>1, mut_ID] %>% unlist()
muts_only_tumour_cont = twins_filtered_dt[sum_normal_mtr_vaf<=3 & sum_normal_mtr_vaf>=1 & sum_tumour_mtr_vaf>5, mut_ID] %>% unlist()
muts_tumour_PD63383 = twins_filtered_dt[sum_normal_PD63383_mtr_vaf<=1 & sum_normal_PD62341_mtr_vaf==0 & sum_tumour_PD63383_mtr_vaf==2, mut_ID] %>% unlist()
muts_only_tumour_PD63383 = setdiff(muts_tumour_PD63383, c(muts_only_tumour, muts_only_tumour_cont)) # 48 
muts_tumour_spec = c(muts_only_tumour, muts_only_tumour_cont, muts_only_tumour_PD63383) # 140

muts_one_sample_normal = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_mtr_vaf==1, mut_ID] %>% unlist()
muts_one_sample_tumour = twins_filtered_dt[sum_tumour_mtr_vaf==1 & sum_normal_mtr_vaf==0, mut_ID] %>% unlist()
muts_one_sample = c(muts_one_sample_normal, muts_one_sample_tumour)

muts_absent_PD63383_all = twins_filtered_dt[sum_normal_PD63383_clean_mtr_vaf==0, mut_ID] %>% unlist() # 282
muts_absent_PD63383 = setdiff(muts_absent_PD63383_all, c(muts_only_normal, muts_all_tumour_nt, muts_tumour_spec, muts_one_sample)) # 13

# Check number of mutations in each category 
paste('Mutations present in all samples:', length(muts_all_samples)) # 78

paste('Number of muts present in all normal samples:', length(muts_all_normal)) # 97 (632 - 535: makes sense)
paste('Number of muts present in all tumour samples:', length(muts_all_tumour)) # 144
paste('Number of muts present in all normal samples but not all tumour samples:', length(muts_all_normal_nt)) # 19 (97 - 78 = 19)
paste('Number of muts present in all tumour samples but not all normal samples:', length(muts_all_tumour_nt)) # 66 (144 - 78 = 66)

paste('Number of muts present in only normal samples:', length(muts_only_normal)) # 20
paste('Number of muts present in only PD62341 normal samples:', length(muts_only_normal_PD62341)) # 1
paste('Number of muts present in only PD63383 normal samples:', length(muts_only_normal_PD63383)) # 7

paste('Number of muts present in only tumour samples:', length(muts_only_tumour)) # 35
paste('Number of muts likely in only tumour samples:', length(muts_only_tumour_cont)) # 81
paste('Number of muts present in PD63383 tumour only:', length(muts_only_tumour_PD63383)) # 24 
paste('Number of muts likely to be tumour-specific:', length(c(muts_tumour_spec))) # 140

paste('Number of muts present in a single sample (normal):', length(muts_one_sample_normal)) # 38 
paste('Number of muts present in a single sample (tumour):', length(muts_one_sample_tumour)) # 99 
paste('Number of muts present in only a single sample:', length(c(muts_one_sample))) # 137

paste('Number of muts absent from PD63383 normal (not yet identified):', length(c(muts_absent_PD63383))) # 10

# Create list with all mutations that have been assigned to a group (groups)
muts_assigned = c(muts_all_samples, muts_all_normal_nt, muts_all_tumour_nt,
                  muts_only_normal, muts_tumour_spec, muts_one_sample, muts_absent_PD63383) %>% unique()
muts_unassigned = setdiff(muts, muts_assigned)

paste('Number of mutations assigned to a group:', length(muts_assigned)) # 449
paste('Number of mutations unassigned to a group:', length(muts_unassigned)) # 150 

######################################################################################################
# Heatmap to show all classes of mutations
mut_all_dt = twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE]
mut_all_mat = as.matrix(mut_all_dt[, c(samples_vaf), with=FALSE])
rownames(mut_all_mat) = mut_all_dt[,1] %>% unlist()  
colnames(mut_all_mat) = tstrsplit(colnames(mut_all_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)),
                            Purity = c(0, 0.1, 0.2, 0, 0.9, 0.3, 0.5, 0.6, 0.8, 0.8, 0, 0, 0, 0, 0, 0.9, 0.9, 0.3, 0.7, 0.4, 0, 0.5)) # fraction of tumour cells 
rownames(col_annotation) = colnames(mut_all_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), 
                         Twin = c(PD62341=col_PD62341, PD63383=col_PD63383),
                         Purity = colorRampPalette(c('#f2e7e7', 'darkred'))(11))

# heatmap
pdf('Results/20241114_p4_heatmap_all_mutations599.pdf')
pheatmap(mut_all_mat,
         cellwidth=10, cellheight=0.4,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations (599)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# save larger format to be able to see mutation names
pdf('Results/20241114_p4_heatmap_all_mutations599_large.pdf', height = 90)
pheatmap(mut_all_mat,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations (599)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# Heatmap: mutations assigned to a group
mut_assigned_dt = twins_filtered_dt[mut_ID %in% muts_assigned, c('mut_ID', samples_vaf), with=FALSE]
mut_assigned_mat = as.matrix(mut_assigned_dt[, c(samples_vaf), with=FALSE])
rownames(mut_assigned_mat) = mut_assigned_dt[,1] %>% unlist()  
colnames(mut_assigned_mat) = tstrsplit(colnames(mut_assigned_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

pdf('Results/20241114_p4_heatmap_assigned_mutations449.pdf')
pheatmap(mut_assigned_mat,
         cellwidth=10, cellheight=0.4,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Assigned mutations (449)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# Heatmap: mutations I didn't assign to any group
mut_unassigned_dt = twins_filtered_dt[!mut_ID %in% muts_assigned, c('mut_ID', samples_vaf), with=FALSE]
mut_unassigned_mat = as.matrix(mut_unassigned_dt[, c(samples_vaf), with=FALSE])
rownames(mut_unassigned_mat) = mut_unassigned_dt[,1] %>% unlist()  
colnames(mut_unassigned_mat) = tstrsplit(colnames(mut_unassigned_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

pdf('Results/20241114_p4_heatmap_unassigned_mutations150.pdf')
pheatmap(mut_unassigned_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Unassigned mutations (150)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

pdf('Results/20241114_p4_heatmap_unassigned_mutations150_large.pdf', height = 30)
pheatmap(mut_unassigned_mat,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Unassigned mutations (150)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Plot tri-nucleotide context of assigned and unassigned mutations

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

# Plot the trinucleotide context for unassigned mutations
mybed_fin = twins_dt[mut_ID %in% muts,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_fin = get_trinucs(mybed_fin, BSgenome.Hsapiens.UCSC.hg38)
dt_fin = twins_dt[mut_ID %in% muts]
dt_fin$trins_fin=trins_fin

# plot the distribution of different mutations across different contexts 
mut_sign_counts_fin = data.table(table(dt_fin[, trins_fin]))
setnames(mut_sign_counts_fin, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_fin[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_fin[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_fin[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_fin, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'All mutations (599)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_all_muts599.pdf', width = 7.5, height = 3.5)

# Plot the trinucleotide distribution for assigned mutations
mybed_asn = twins_dt[mut_ID %in% muts_assigned,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_asn = get_trinucs(mybed_asn, BSgenome.Hsapiens.UCSC.hg38)
dt_asn = twins_dt[mut_ID %in% muts_assigned]
dt_asn$trins_asn=trins_asn

# plot the distribution of different mutations across different contexts 
mut_sign_counts_asn = data.table(table(dt_asn[, trins_asn]))
setnames(mut_sign_counts_asn, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_asn[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_asn[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_asn[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_asn, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Assigned mutations (449)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_assigned_muts.pdf', width = 7.5, height = 3.5)

# Plot the trinucleotide context for unassigned mutations
mybed_usn = twins_dt[mut_ID %in% muts_unassigned,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_usn = get_trinucs(mybed_usn, BSgenome.Hsapiens.UCSC.hg38)
dt_usn = twins_dt[mut_ID %in% muts_unassigned]
dt_usn$trins_usn=trins_usn

# plot the distribution of different mutations across different contexts 
mut_sign_counts_usn = data.table(table(dt_usn[, trins_usn]))
setnames(mut_sign_counts_usn, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_usn[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_usn[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_usn[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_usn, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Unassigned mutations (150)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_unassigned_muts.pdf', width = 7.5, height = 3.5)

# Trinucleotide context for different groups
# Present in all normal samples
# Plot the trinucleotide context for unassigned mutations
mybed1 = twins_dt[mut_ID %in% muts_all_normal,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins1 = get_trinucs(mybed1, BSgenome.Hsapiens.UCSC.hg38)
dt1 = twins_dt[mut_ID %in% muts_all_normal]
dt1$trins1=trins1

# plot the distribution of different mutations across different contexts 
mut_sign_counts1 = data.table(table(dt1[, trins1]))
setnames(mut_sign_counts1, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts1[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts1[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts1[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts1, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations in all normal samples (97)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_allNormal_muts.pdf', width = 7.5, height = 3.5)

# Likely tumour specific 
mybed2 = twins_dt[mut_ID %in% muts_tumour_spec,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins2 = get_trinucs(mybed2, BSgenome.Hsapiens.UCSC.hg38)
dt2 = twins_dt[mut_ID %in% muts_tumour_spec]
dt2$trins2=trins2

# plot the distribution of different mutations across different contexts 
mut_sign_counts2 = data.table(table(dt2[, trins2]))
setnames(mut_sign_counts2, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts2[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts2[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts2[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts2, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations likely specific to the tumour (140)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_tumourSpec_muts.pdf', width = 7.5, height = 3.5)

# PD63383 tumour specific 
# Plot the trinucleotide context for unassigned mutations
mybed3 = twins_dt[mut_ID %in% muts_tumour_PD63383,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins3 = get_trinucs(mybed3, BSgenome.Hsapiens.UCSC.hg38)
dt3 = twins_dt[mut_ID %in% muts_tumour_PD63383]
dt3$trins3=trins3

# plot the distribution of different mutations across different contexts 
mut_sign_counts3 = data.table(table(dt3[, trins3]))
setnames(mut_sign_counts3, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts3[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts3[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts3[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts3, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations specific to PD63383 tumour (24)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_tumourPD63383_muts.pdf', width = 7.5, height = 3.5)

# PD63383 tumour specific 
# Plot the trinucleotide context for unassigned mutations
mybed4 = twins_dt[mut_ID %in% muts_one_sample,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins4 = get_trinucs(mybed4, BSgenome.Hsapiens.UCSC.hg38)
dt4 = twins_dt[mut_ID %in% muts_one_sample]
dt4$trins4=trins4

# plot the distribution of different mutations across different contexts 
mut_sign_counts4 = data.table(table(dt4[, trins4]))
setnames(mut_sign_counts4, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts4[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts4[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts4[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts4, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Mutations present in a single sample (137)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p4_mut_trins_singleSample_muts.pdf', width = 7.5, height = 3.5)

######################################################################################################
# CLASSES OF MUTATIONS 1: MUTATIONS PRESENT IN ALL NORMAL SAMPLES 
# Inspect mutations that are present in all samples
# Why were these not called as germline?

# For all mutations, VAF is clearly < 0.5
pdf('Results/20241114_p4_hist_vaf_97muts_in_all_normal_samples.pdf', height = 4, width = 4)
hist(twins_filtered_dt[mut_ID %in% muts_all_normal, vaf_all_normal] %>% unlist(), 
     xlab = 'VAF (agg normal samples)', main = '97 muts: all normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_filtered_dt[mut_ID %in% muts_all_normal, vaf_all_normal]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241114_p4_hist_vaf_97muts_in_PD62341_normal_samples.pdf', height = 4, width = 4)
hist(twins_filtered_dt[mut_ID %in% muts_all_normal, vaf_all_normal_PD62341] %>% unlist(), 
     xlab = 'VAF (agg PD62341 normal samples)', main = '97 muts: PD62341 normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_filtered_dt[mut_ID %in% muts_all_normal, vaf_all_normal_PD62341]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241114_p4_hist_vaf_97muts_in_PD63383_normal_samples.pdf', height = 4, width = 4)
hist(twins_filtered_dt[mut_ID %in% muts_all_normal, vaf_all_normal_PD63383] %>% unlist(), 
     xlab = 'VAF (agg PD63383 normal samples)', main = '97 muts: PD63383 normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_filtered_dt[mut_ID %in% muts_all_normal, vaf_all_normal_PD63383]), col = 'purple', lwd = 2.5)
dev.off()

# plot aggregated VAF vs coverage
ggplot(twins_filtered_dt[mut_ID %in% muts_all_normal], aes(x = dep_all_normal, y = vaf_all_normal))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'VAF (agg normal)', y = 'Coverage (agg normal)')+
  ggtitle(glue('97 mutations (all normal samples)'))+
  xlim(c(0, 600))+
  ylim(c(0, 0.5))
ggsave(glue('Results/20241114_p4_vaf_vs_dep_allNormal.pdf'), heigh = 4, width = 4)

# plot VAF in PD62341 vs PD63383 (are any of those mutations enriched?)
ggplot(twins_filtered_dt[mut_ID %in% muts_all_normal], aes(x = vaf_all_normal_PD62341, y = vaf_all_normal_PD63383))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'VAF (agg PD62341 normal)', y = 'VAF (agg PD63383 normal)')+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  ggtitle(glue('97 mutations (all normal samples)'))+
  coord_equal(ratio=1)
ggsave(glue('Results/20241114_p4_vaf_normal_PD62341_vs_PD63383.pdf'), heigh = 4, width = 4)

# identify enriched mutations 

# log fold change 
twins_filtered_dt[, ratio_PD62341_PD63383_normal := vaf_all_normal_PD62341 / vaf_all_normal_PD63383]
twins_filtered_dt[, ratio_PD62341_PD63383_normal_log := log2(ratio_PD62341_PD63383_normal)]

pdf('Results/20241114_p4_hist_vaf_97muts_in_all_normal_samples_ratio_PD62341toPD63383.pdf', heigh = 4, width = 4)
hist(twins_filtered_dt[mut_ID %in% muts_all_normal, ratio_PD62341_PD63383_normal] %>% unlist(), 
     xlab = 'VAF ratio PD62341 / PD63383', main = '97 mutations', breaks=30)
dev.off()

pdf('Results/20241114_p4_hist_vaf_97muts_in_all_normal_samples_ratio_PD62341toPD63383_log.pdf', heigh = 4, width = 4)
hist(twins_filtered_dt[mut_ID %in% muts_normal_all, ratio_PD62341_PD63383_normal_log] %>% unlist(), 
     xlab = 'log2(VAF ratio PD62341 / PD63383)', main = '97 mutations', breaks=30)
abline(v = -0.3, col = 'red', lwd = 2.5)
abline(v = 0.3, col = 'red', lwd = 2.5)
dev.off()

twins_filtered_dt[, mut_cat := factor(fcase(
  vaf_all_normal_PD62341 > 0.45 & vaf_all_normal_PD63383 > 0.45, 'putative germline', 
  ratio_PD62341_PD63383_normal_log < -0.3, 'PD63383 enriched',
  ratio_PD62341_PD63383_normal_log > 0.3, 'PD62341 enriched',
  ratio_PD62341_PD63383_normal_log < 0.3 & ratio_PD62341_PD63383_normal_log > -0.3, 'no difference'))]

ggplot(twins_filtered_dt[mut_ID %in% muts_normal_all], aes(x = vaf_all_normal_PD62341, y = vaf_all_normal_PD63383, col = mut_cat))+
  geom_point(size=1.5, alpha = 1)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF agg PD62341 normal', y = 'VAF agg PD63383 normal', col = 'Mutation category')+
  scale_color_manual(values = c('grey', col_PD62341, col_PD63383, 'darkred'))+
  ggtitle(glue('97 mutations (all normal samples)'))+
  xlim(c(0, 0.6))+
  ylim(c(0, 0.6))+
  coord_equal(ratio=1)
ggsave(glue('Results/20241114_p4_vaf_all_normal_PD62341_PD63383_col_cat.pdf'), width = 6, height = 4)

# inspect those mutations in detail 
twins_filtered_dt[mut_ID %in% muts_all_normal & mut_cat == 'putative germline', mut_ID] # 1
# "chr12_131428827_C_A" # poor mapping 
twins_filtered_dt[mut_ID %in% muts_normal_all & mut_cat == 'PD62341 enriched', mut_ID] # 3
# "chr10_13474524_G_A" # poor quality (deletions)
# "chr11_4219150_T_G"  # poor mapping 
# "chr13_38383628_T_C" # looks quite okay
# "chr18_56506391_G_C" # poor quality (split reads)
# "chr19_50081651_T_C" # poor mapping 
# "chr1_16864604_G_A" # poor mapping 
# "chr2_87964592_G_T" # poor mapping 
# "chr6_81469373_G_A"  # split reads 
# "chr8_123201027_T_C" # looks real
# "chr8_232785_G_C" # poor mapping of (-) strands
twins_filtered_dt[mut_ID %in% muts_normal_all & mut_cat == 'PD63383 enriched', mut_ID] # 1
# "chr10_79789115_C_T" # poor mapping 
# "chr1_38827952_C_A"  # looks real 
# "chr22_34554733_G_A" # poor mapping 
# "chr7_5942730_C_T"  # poor mapping 
# "chr9_135524032_G_A" # poor mapping 
# "chr9_63875526_T_G"  # poor quality (indels)

# plot the distribution of VAFs across samples for the good looking mutations
twins_vaf_melt = melt(twins_filtered_vaf[, c('mut_ID', samples_vaf), with=FALSE], id.vars = 'mut_ID')
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

ggplot(twins_vaf_melt[mut_ID=='chr8_123201027_T_C'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ggtitle(glue('Mutation: chr8:123201027, T>C'))+
  ylim(c(0, 0.75))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p4_allNormal_dist_samples_mut_chr8_PD62341_enriched.pdf'), height = 3, width = 6.5)

ggplot(twins_vaf_melt[mut_ID=='chr1_38827952_C_A'], aes(x = sample, y = value, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Sample')+
  ylim(c(0, 0.75))+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241114_p4_allNormal_dist_samples_mut_chr1_PD63383_enriched.pdf'), height = 3, width = 6.5)

# plot distribution of mutations present in all normal samples across the genome 
Chrom = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
          "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
          "chr20", "chr21", "chr22", "chrX")
# Source: https://www.ncbi.nlm.nih.gov/grc/human/data (Hg38), accessed 14/11/2024
Chrom_length = as.numeric(c(248956422, 242193529, 198295559, 190214555,
                            181538259, 170805979, 159345973, 145138636,
                            138394717, 133797422, 135086622, 133275309,
                            114364328, 107043718, 101991189, 90338345,
                            83257441, 80373285, 58617616, 64444167,
                            46709983, 50818468, 156040895))
chr_lengths_dt = data.table(cbind(Chrom, Chrom_length))

# plot against all mutations identified in the set 
twins_dt = merge(twins_dt, chr_lengths_dt, by = 'Chrom')
twins_dt[, Chrom := factor(Chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                             "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                             "chr20", "chr21", "chr22", "chrX"))]
twins_dt[, pos := tstrsplit(mut_ID, '_', fixed=T, keep=2)]
twins_dt[, pos := as.numeric(pos)]

twins_dt[, cat_all_normal := factor(fcase(
  mut_ID %in% muts_all_normal, 'present in all normal samples (97)',
  !mut_ID %in% muts_all_normal, 'other'
))]

twins_dt[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]

for (chr in Chrom){
  dt = twins_dt[mut_ID %in% muts_qc & Chrom == chr] # exclude mutations removed due to poor quality 
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt %>% arrange(cat_all_normal), aes(x = pos, y = dep_all_normal, col=cat_all_normal, alpha=cat_all_normal, order=cat_all_normal))+
    geom_point(size=1.5)+
    scale_color_manual(values = c('grey', col_normal))+
    scale_alpha_discrete(range = c(0.1, 0.9))+
    guides(alpha = "none")+
    theme_classic(base_size = 12)+
    labs(x = 'Genomic position', y = 'Total coverage (normal samples)', col='Mutation category')+
    ggtitle(glue('{chr}, {nr} mutations'))+
    guides(color = guide_legend(override.aes = list(alpha = 1) ) )+
    xlim(c(0, length))
  ggsave(glue('Results/20241114_p4_dist_allNormal97_{chr}.pdf'), height=3, width=6.5)
}

######################################################################################################
# MUTATION CLASS 2: Mutations present in all normal, but not all tumour samples

muts_all_normal_nt # 19
# "chr10_79512492_C_G" # poor quality 
# "chr13_38383628_T_C" # likely present in all samples: looks okay (some mis-mapping of minus strands and insertions) # 3 reads in PD63383aq
# "chr16_20575297_A_G" # poor quality 
# "chr16_88009831_A_G" # likely present in all samples: some insertions (VAF 0.09 in PD62341ap - has lots of insertion and only 3 reads with the mutation)
# "chr18_56506391_G_C" # poor quality: insertions + split reads 
# "chr19_40853394_C_A" # poor quality mapping 
# "chr19_50081651_T_C" # poor mapping of (-) strand reads, PD62341 has 5 reads (3 on better-mapping reads)
# "chr1_103717717_G_A" # poor quality mapping 
# "chr1_13269066_T_C"  # poor quality mapping
# "chr1_21423048_C_A"  # poor quality mapping of (-) strand reads
# "chr1_38827952_C_A"  # loss of chr1 copy in PD62341am, PD63383ap, PD63383aq
# "chr1_47072883_A_T"  # poor quality mapping 
# "chr2_1440752_T_C"   # likely presence in all samples: looks okay - 5 reads in PD62341ak and VAF 0.09 
# "chr3_190409217_C_T" # likely presence in all samples: looks okay - 0.064 in PD62341u, but high-quality reads
# "chr6_94069453_A_T"  # poor quality (all 'mutations' likely due to insertions)
# "chr9_133062780_G_C" # poor mapping  
# "chr9_77762599_A_T"  # poor mapping 
# "chrX_140138202_C_T" # poor mapping of (-) strand reads (explains lower VAF in PD62341u)
# "chrX_66502763_G_A" # likely presence in all samples: looks good (3 reads in am and ap, likely due to lower coverage in those samples)

# Main conclusion: all of those are either artefacts or present in all samples 

######################################################################################################
# MUTATION CLASS 3: Mutations present in only normal samples 

muts_only_normal # 20
# "chr11_34011887_C_T" # looks okay
# "chr13_24392029_C_G" # poor quality 
# "chr13_50815806_A_G" # looks okay
# "chr17_20476029_C_T" # poor mapping 
# "chr18_56262973_G_A" # poor quality (mapping reverse strands + indels)
# "chr1_103587565_A_C" # looks real (3 PD63383 samples)
# "chr21_28681933_A_C" # poor quality 
# "chr21_40193588_G_A" # looks okay
# "chr3_165901319_C_A" # to consider - mapping is quite okay but maps to several places on blat
# "chr3_77633967_C_T" # looks okay
# "chr4_15905566_C_T" # looks okay (present at low levels, but mapping is good quality, I do not note any insertions, no double-mapping on blat) 
# "chr4_74625500_G_T" # looks okay 
# "chr4_75704880_G_A" # looks okay
# "chr6_159851462_T_C" # poor quality (deletions) 
# "chr6_165179306_G_A" # looks okay
# "chr7_110349267_G_A" # plausible (+ strand reads map with poor quality)
# "chr7_143846161_G_A" # poor mapping 
# "chr7_149688370_C_T" # looks okay
# "chr7_73831920_C_T"  # looks okay
# "chrX_115066661_C_T" # looks okay 

muts_only_normal_val = c("chr11_34011887_C_T", "chr13_50815806_A_G", "chr1_103587565_A_C", "chr21_40193588_G_A", 
                         "chr3_165901319_C_A", "chr3_77633967_C_T", "chr4_15905566_C_T", "chr4_74625500_G_T", 
                         "chr4_75704880_G_A",  "chr6_165179306_G_A", "chr7_149688370_C_T", "chr7_73831920_C_T", "chrX_115066661_C_T")

twins_vaf_sub = twins_filtered_vaf[mut_ID %in% muts_only_normal_val, c('mut_ID', samples_vaf), with=FALSE]
twins_vaf_melt = melt(twins_vaf_sub, id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_', keep = 1, fixed = TRUE)]
twins_vaf_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_melt[, sample_type := paste(status, twin, sep = ', ')]
twins_vaf_melt[, sample_type1 := factor(fcase(
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

# calculate the median for 5 clean PD63383 normal samples and arrange mutations by this value
twins_vaf_sub[, median_PD63383_clean := apply(.SD, 1, function(x) median(x)), .SDcols = samples_normal_PD63383_clean_vaf]
twins_vaf_sub[, median_PD63383 := apply(.SD, 1, function(x) median(x)), .SDcols = samples_normal_PD63383_vaf]
twins_vaf_sub[, median_PD62341_clean := apply(.SD, 1, function(x) median(x)), .SDcols = c('PD62341aa_VAF', 'PD62341q_VAF', 'PD62341h_VAF', 'PD62341n_VAF', 'PD62341ad_VAF')]

twins_vaf_sub[, mean_PD63383_clean := apply(.SD, 1, function(x) mean(x)), .SDcols = samples_normal_PD63383_clean_vaf]
twins_vaf_sub[, mean_PD63383 := apply(.SD, 1, function(x) mean(x)), .SDcols = samples_normal_PD63383_vaf]
twins_vaf_sub[, mean_PD62341_clean := apply(.SD, 1, function(x) mean(x)), .SDcols = c('PD62341aa_VAF', 'PD62341q_VAF', 'PD62341h_VAF', 'PD62341n_VAF', 'PD62341ad_VAF')]

ggplot(twins_vaf_sub, aes(x = median_PD63383, y = median_PD62341_clean))+
  geom_point(size=2.5)+
  theme_classic(base_size = 12)+
  labs(x = 'Median in PD63383', y = 'Median in PD62341 (excluded PD62341v)')+
  ggtitle(glue('Mutations only in normal samples (13)'))+
  xlim(c(0, 0.25))+
  ylim(c(0, 0.25))+
  coord_equal(ratio=1)
ggsave('Results/20241114_p4_muts_onlyNormal13_medians.pdf', height = 4, width = 4)

ggplot(twins_vaf_sub, aes(x = mean_PD63383, y = mean_PD62341_clean))+
  geom_point(size=2.5)+
  theme_classic(base_size = 12)+
  labs(x = 'Mean in PD63383', y = 'Mean in PD62341 (excluded PD62341v)')+
  ggtitle(glue('Mutations only in normal samples (13)'))+
  xlim(c(0, 0.25))+
  ylim(c(0, 0.25))+
  coord_equal(ratio=1)
ggsave('Results/20241114_p4_muts_onlyNormal13_means.pdf', height = 4, width = 4)

twins_vaf_melt = merge(twins_vaf_melt, twins_vaf_sub[, c('mut_ID', 'median_PD63383_clean'), with=FALSE], by = 'mut_ID')
twins_vaf_melt[, mut_ID := factor(mut_ID, levels=unique(mut_ID[order(median_PD63383_clean, mut_ID)]), ordered = TRUE)]

ggplot(twins_vaf_melt, aes(x = value, y = mut_ID, col = sample_type1))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations present only in normal samples (13)'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))
ggsave('Results/20241114_p4_muts_onlyNormal13.pdf', height = 5, width = 7.5)

twins_vaf_melt[, sample_type2 := factor(fcase(
  sample == 'PD62341v', 'PD62341, normal: spleen',
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]
ggplot(twins_vaf_melt, aes(x = value, y = mut_ID, col = sample_type2))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations present only in normal samples (13)'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD62341_spleen, col_PD63383, col_tumour))
ggsave('Results/20241114_p4_muts_onlyNormal13_spleen.pdf', height = 5, width = 7.5)

twins_vaf_melt[, sample_type3 := factor(fcase(
  sample == 'PD62341v', 'PD62341, normal: spleen',
  sample == 'PD63383bb', 'PD63383, normal: skin',
  sample == 'PD63383w', 'PD63383, normal: spleen',
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

ggplot(twins_vaf_melt %>% arrange(sample_type3), aes(x = value, y = mut_ID, col = sample_type3))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations present only in normal samples (13)'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD62341_spleen, col_PD63383, col_PD63383_skin, col_PD63383_spleen, col_tumour))
ggsave('Results/20241114_p4_muts_onlyNormal13_spleen_skin.pdf', height = 5, width = 7.5)

######################################################################################################
# MUTATION CLASS 4: Mutations absent from clean PD63383 normal samples

muts_absent_PD63383
# "chr10_46846776_C_T"  # poor quality
# "chr12_114403258_G_A" # looks okay (4 PD62341 tumour, 1 PD62341 normal - likely tumour specific subclonal)
# "chr13_18565005_G_A"  # poor quality 
# "chr17_1199637_C_G"   # poor quality (split reads)
# "chr18_1379467_G_T"  # looks okay (4 PD62341 tumour, 1 PD62341 normal - likely tumour specific subclonal)
# "chr1_12834696_C_T"  # poor quality ((-) strand maps very poorly)  
# "chr1_148178303_A_G" # poor mapping 
# "chr20_54999342_C_A" # looks okay (seems lower in the tumour too?)   
# "chr5_71157247_C_T"  # poor mapping   
# "chr6_57048152_G_A"  # poor mapping + low coverage 
# "chr6_8310684_G_T"   # looks okay (5 PD62341 tumour, 1 PD62341 normal - likely tumour specific subclonal)
# "chr8_137568368_T_A" # looks okay (9 tumour, 3 PD62341 normal q, aa, h - possibly tumour specific)
# "chrX_135867449_G_A" # very poor mapping

muts_absent_PD63383_val = c("chr12_114403258_G_A", "chr18_1379467_G_T", "chr20_54999342_C_A",  "chr6_8310684_G_T", "chr8_137568368_T_A")

twins_vaf_melt = melt(twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE])
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


for (mut in muts_absent_PD63383_val){
  ggplot(twins_vaf_melt[mut_ID==mut], aes(x = sample, y = value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'Sample', y = 'VAF', col = 'Sample')+
    ggtitle(glue('Mutation: {mut}'))+
    ylim(c(0, 0.6))+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241114_p4_PD63383absent_{mut}.pdf'), height = 3, width = 6.5)
}

######################################################################################################
# CLASSES OF MUTATIONS 5: MUTATIONS PRESENT IN ALL TUMOUR BUT ONLY SOME NORMAL SAMPLES

# lists of relevant mutations 
muts_tumour_all_nt
muts_tumour_all_nt_adj

# how many normal samples are these usually present in?
pdf('Results/20241109_p7_hist_nr_normal_66muts_tumour_nt.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt, sum_normal],
     xlab = 'Number of normal samples', main = '66 mutations (all tumour, not all normal)')
dev.off()

# heatmap to have a look at how this looks like generally
mut_tumour_nt = twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt, c('mut_ID', samples_vaf), with=FALSE]
mut_tumour_nt_mat = as.matrix(mut_tumour_nt[, c(samples_vaf), with=FALSE])
rownames(mut_tumour_nt_mat) = mut_tumour_nt[,1] %>% unlist()  
colnames(mut_tumour_nt_mat) = tstrsplit(colnames(mut_tumour_nt_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_tumour_nt_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_all_tumour_some_normal.pdf')
pheatmap(mut_tumour_nt_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="66 muts: all tumour, but not all normal", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

pdf('Results/20241109_p6_heatmap_all_tumour_some_normal_rownames.pdf', height = 80)
pheatmap(mut_tumour_nt_mat,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="66 muts: all tumour, but not all normal", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# heatmap for adjusted vaf values
mut_tumour_nt_adj = twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj, c('mut_ID', samples_normal_adj, samples_tumour_adj), with=FALSE]
mut_tumour_nt_mat_adj = as.matrix(mut_tumour_nt_adj[, c(samples_normal_adj, samples_tumour_adj), with=FALSE])
rownames(mut_tumour_nt_mat_adj) = mut_tumour_nt_adj[,1] %>% unlist()  
colnames(mut_tumour_nt_mat_adj) = tstrsplit(colnames(mut_tumour_nt_mat_adj), '_adj', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_tumour_nt_mat_adj)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_all_tumour_some_normal_vaf_adj.pdf')
pheatmap(mut_tumour_nt_mat_adj,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="66 muts: all tumour, but not all normal, adj VAF", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# are there any mutations present in PD63383 normal samples?
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 0] # 11 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 1] # 21
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 2] # 1

# are there any mutations that are only shared with normal PD63383?
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD62341 == 0] # 0

# get out mutations absent from PD63383 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1] # 33 

# which normal PD62341 samples are these shared with?
pdf('Results/20241109_p7_hist_nr_normal_66muts_tumour_nt_PD62341.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1, sum_normal_PD62341],
     xlab = 'Number of normal PD62341 samples with mutation', main = '66 mutations (all tumour, not all normal)') # 33 
dev.off()

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 2] # 8 
# likely all contamination
# chr10_100754461_C_T aa, h
# chr14_104500332_C_T aa, h
# chr15_23462705_C_A aa, h
# chr17_40061856_G_C aa, h
# chr4_147908994_T_C aa, h
# chr6_95827754_C_A aa, h
# chrX_66719643_C_G aa, h
# chrX_68803487_C_T aa, h
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 3] # 1 
# chr5_44907911_G_A q, aa, h (adjusted VAF in q = 0.086, so maybe not in the end)
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 4] # 2 
# chr3_62055057_C_G q, aa, h, n
# chr3_62055077_G_C q, aa, h, n 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 5] # 3 
# chr16_5479739_C_T q, aa, ad, h, n
# chr2_95662131_G_A q, aa, ad, h, n
# chr3_50106043_C_T q, aa, ad, h, n
# absent from v - but it is likely that MTR is lower due to twin-twin contamination 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 <= 1 & sum_normal_PD62341 == 6] # 2
# chr14_105458006_C_A # all PD62341
# chr17_33422229_C_A # all PD62341

twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 1] # 13
# only sharing with aa or h
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 4] # 2
# chr3_62055057_C_G q, aa, h, n 
# chr3_62055077_G_C q, aa, h, n 
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 5] # 3
# chr16_5479739_C_T q, aa, ad, h, n
# chr2_95662131_G_A q, aa, ad, h, n
# chr3_50106043_C_T q, aa, ad, h, n
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 6] # 2
# chr14_105458006_C_A # all PD62341
# chr17_33422229_C_A # all PD62341

twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 1] # 13
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 4] # 2
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 5] # 3
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj <= 1 & sum_normal_PD62341_adj == 6] # 2

twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 1] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 2] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 3] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 4] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt_adj & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 5] # 0
twins_filtered_vaf_adj[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383_adj == 2 & sum_normal_PD62341_adj == 6] # 1
# chr15_49480646_T_A (bb and t - why t?)

# plot VAF 

######################################################################################################
# MUTATION CLASSES 6: MUTATIONS IN ONLY TUMOUR SAMPLES

ggplot(twins_vaf_adj_melt[mut_ID %in% muts_tumour_only], aes(x = value, y = mut_ID, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only tumour samples): 31'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour_PD62341, col_tumour_PD63383))
ggsave('Results/20241109_p4_muts_only_tumour_samples_vaf_adj.pdf', height = 6.5, width = 7.5)

twins_vaf_melt = melt(twins_filtered_vaf[, c('mut_ID', samples_vaf), with=FALSE], id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_vaf_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_melt[, sample_type := paste(status, twin, sep = ', ')]
twins_vaf_melt[, sample_type1 := factor(fcase(
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

ggplot(twins_vaf_melt[mut_ID %in% muts_tumour_only], aes(x = value, y = mut_ID, col = sample_type))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only tumour samples): 31'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour_PD62341, col_tumour_PD63383))
ggsave('Results/20241109_p4_muts_only_tumour_samples_vaf.pdf', height = 6.5, width = 7.5)

######################################################################################################
# MUTATION CLASSES 6B: MUTATIONS ONLY IN PD63383 TUMOUR

# what about mutations present in PD63383 only?
twins_filtered_vaf[sum_tumour_PD63383 > 0 & sum_tumour_PD62341 == 0] 
twins_filtered_vaf_adj[sum_tumour_PD63383 > 0 & sum_tumour_PD62341 == 0] 
# everything looks great unless stated otherwise 
# the same set of mutations comes up in adj and non-adj VAF
# "chr10_100308403_C_T" 
# "chr10_31333522_G_A"  
# "chr11_80115755_C_T"  
# "chr12_43132316_G_A" 
# "chr12_80762031_G_A"  
# "chr19_43348626_C_A" # not as good but probably fine  
# "chr1_160088990_G_A"  
# "chr1_51714096_G_C"  
# "chr22_25281865_C_T" 
# "chr22_43290224_C_T" # double check
# "chr2_72939205_C_T"   
# "chr2_82147172_T_C"   
# "chr3_137508691_C_A" 
# "chr4_179587218_G_A"  
# "chr5_157248612_A_G" # present in only 1 sample (VAF = 0.09, but 4 reads)
# "chr5_28014472_C_T"   
# "chr5_54017866_G_T" 
# "chr6_86549064_C_T"  
# "chr6_92179408_A_G"   
# "chr7_13677323_A_G"   
# "chr7_148446413_G_A"  
# "chr7_49732059_G_A"  
# "chr9_98275090_C_T"   
# "chrX_117418608_G_A"  
# "chrX_119915655_G_A"  
# "chrX_36779694_C_T"  
# "chrX_48453603_G_A"   
# "chrX_70352624_G_A"   
# "chrX_72671786_A_T" # looks plausible but really maps to EVERYWHERE on Blat    
# "chrX_77681162_G_A"

muts_PD63383_tumour_only = twins_filtered_vaf_adj[sum_tumour_PD63383 > 0 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist() 

# in PD63383 only but not tumour only 
# these are due to persisting contamination in the skin (even after adjusting for contamination)
setdiff(muts_PD63383_tumour_only, muts_tumour_only) # 18

# in tumour only but not only in PD63383
# are these mutations only present in PD62341 tumour?
setdiff(muts_tumour_only, muts_PD63383_tumour_only) # 19
twins_filtered_vaf[mut_ID %in% setdiff(muts_tumour_only, muts_PD63383_tumour_only), 
         c('mut_ID', 'sum_tumour', 'sum_tumour_PD62341')]
# chr8_21444725_G_A # present in both PD63383 tumours and 2 PD62341 tumour samples
# chr7_152562130_G_A # present in both PD63383 tumours and 5 PD62341 tumour samples

twins_filtered_vaf[mut_ID %in% setdiff(muts_tumour_only, muts_PD63383_tumour_only), c('mut_ID', samples_tumour_PD62341_vaf), with=FALSE]
# note that this is with un-adjusted values and some samples we estimate to be of very low purity
# "chr10_72139813_C_T" ak, ap // ak, ap
# "chr11_43029745_A_T" ae, aj, b // ae, aj, b, u
# "chr12_19555916_C_T" ae, aj // ae, aj, u
# "chr13_28953051_G_T" ae, ag, aj // ae, ag, aj, b, u
# "chr13_62110424_G_A" am, ap // am, ap
# "chr15_68707110_G_A" ak, ap // ak, ap
# "chr16_17805832_C_T" ae, ag, aj, b, u // ae, ag, aj, b, u
# "chr1_14644486_C_A" ae, aj // ae, aj, b
# "chr1_38311094_G_T" ak, ap // ak, ap
# "chr2_231591536_G_T" ae, am // ae, am, u
# "chr3_189240074_C_A" am // ak, am, ap
# "chr3_95470911_T_A" ak // ak, ap
# "chr5_58182780_T_A" ak, ap // ak, ap
# "chr7_152562130_G_A" ae, aj, am, ap, b // ae, aj, am, ak, ap, b
# "chr7_153784080_C_G" ae, ag, am, b // ae, ag, ak, am, ap, b
# "chr8_21444725_G_A" am, b // am, b, u
# "chr8_91252357_A_C" ak // ak, ap
# "chrX_64561555_C_T" ak, ap // ag, ak, ap
# "chrX_9804641_G_A" ak, ap // ak, ap

# can I get a heatmap of this?
PD62341_subclonal = twins_filtered_vaf_adj[mut_ID %in% setdiff(muts_tumour_only, muts_PD63383_tumour_only), c('mut_ID', samples_tumour_vaf), with=FALSE]
PD62341_subclonal_mat = as.matrix(PD62341_subclonal[, c(samples_tumour_vaf), with=FALSE])
rownames(PD62341_subclonal_mat) = PD62341_subclonal[,1] %>% unlist()  
colnames(PD62341_subclonal_mat) = tstrsplit(colnames(PD62341_subclonal_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Twin = c(rep('PD62341', 6), rep('PD63383', 2), rep('PD62341', 2)))
rownames(col_annotation) = colnames(PD62341_subclonal_mat)
annotation_colors = list(Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# annotation data frame 
pdf('Results/20241109_p6_heatmap_PD62341_subclonal.pdf')
pheatmap(PD62341_subclonal_mat,
         cellwidth=10, cellheight=5,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="19 PD62341 subclonal mutations", 
         legend = T,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Analysis of tumour evolution: coverage vs VAF plot

twins_agg_vaf[, mut_class := factor(fcase(
  mut_ID %in% muts_all, 'all samples',
  mut_ID %in% muts_normal_all_nt, 'all normal samples, not all tumour',
  mut_ID %in% muts_tumour_all_nt, 'all tumour samples, not all normal',
  mut_ID %in% muts_normal_only, 'normal only',
  mut_ID %in% muts_tumour_only, 'tumour only',
  !mut_ID %in% c(muts_tumour_only, muts_normal_only, muts_all, 
                 muts_normal_all_nt, muts_tumour_all_nt), 'other')) ]
twins_agg_vaf[, mut_class := factor(mut_class, levels = c(
  'all samples', 'all normal samples, not all tumour', 'all tumour samples, not all normal',
  'normal only', 'tumour only', 'other'))]

myColors <- c('lightblue', 'purple', '#07a94b', '#f3782b',  'darkred', 'grey')
names(myColors) <- levels(twins_agg_vaf[, mut_class])

# Plot coverage vs VAF 
ggplot(twins_agg_vaf %>% arrange(mut_class), aes(x = dep_tumour_all, y = vaf_tumour_all, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  scale_color_manual(values = myColors)+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002.pdf'), height = 5, width = 6)

ggplot(twins_agg_vaf %>% arrange(mut_class), aes(x = dep_normal_all, y = vaf_normal_all, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all normal samples)', y = 'VAF (all normal samples)', col = 'Mutation class')+
  ggtitle(glue('Coverage vs VAF across normal samples'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  scale_color_manual(values = myColors)
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002_normal.pdf'), height = 5, width = 6)

# color by type of mutation
twins_agg_vaf[, c('Ref', 'Alt') := tstrsplit(mut_ID, '_', fixed=T, keep=3:4)]
twins_agg_vaf[, mut_type := paste0(Ref, '>', Alt)]
twins_agg_vaf[, mut_type := as.factor(fcase(
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

ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = mut_type))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_mut_type.pdf'))

# color by chromosome 
twins_agg_vaf[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = Chrom))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_chrom.pdf'))

twins_agg_vaf[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'))]

ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = loss))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = 'Chromosome')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_loss.pdf'))

# add number of tumour samples the mutation is present in 
twins_agg_vaf = merge(twins_agg_vaf, twins_filtered_vaf[, c('mut_ID', 'sum_tumour'), with=FALSE], by = 'mut_ID')
twins_agg_vaf[, sum_tumour := as.factor(sum_tumour)]
ggplot(twins_agg_vaf, aes(x = dep_tumour_all, y = vaf_tumour_all, col = sum_tumour))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (all tumour samples)', y = 'VAF (all tumour samples)', col = '# tumour samples\nmut detected in')+
  ggtitle(glue('Coverage vs VAF across the tumour'))+
  geom_hline(yintercept = 0.35, col = 'black')+
  annotate('text', 525, 0.365, label = 'Estimated clonal VAF', col = 'black')
ggsave(glue('Results/20241109_cov_vs_vaf_sum_tumour.pdf'))

######################################################################################################
# Plot VAF and coverage across the genome 
twins_agg_vaf[, pos := tstrsplit(mut_ID, '_', fixed=T, keep=2)]
twins_agg_vaf[, pos := as.numeric(pos)]
twins_agg_vaf[, Chrom_nr := tstrsplit(Chrom, 'chr', fixed=T, keep=2)]
twins_agg_vaf[, Chrom_nr := factor(Chrom_nr, levels = 
                                        c('1', '2', '3', '4', '5',
                                          '6', '7', '8', '9', '10',
                                          '11', '12', '13', '14', '15',
                                          '16', '17', '18', '19', '20',
                                          '21', '22', 'X'))]

ggplot(twins_agg_vaf, aes(x = pos, y = vaf_tumour_all, col=mut_class))+
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
ggsave(glue('Results/20241109_vaf_across_the_genome_1002.pdf'))

######################################################################################################
# Plot distribution of mutations present in all tumour samples across the genome 

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
  ggsave(glue('Results/20241109_dist_across_genome_641_tumour_{chr}.pdf'), height=3, width=4.5)
}

twins_agg[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all_nt, 'all normal, but not all tumour',
  mut_ID %in% muts_tumour_all_nt, 'all tumour, but not all normal samples',
  mut_ID %in% muts_all, 'all samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  mut_ID %in% muts_normal_only, 'only normal samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_all, muts_normal_only, muts_tumour_all_only, muts_tumour_only), 'other'))]

# make sure colors are always the same
myColors <- c('darkred', 'purple', '#07a94b', '#f3782b',  '#2fb1f3', 'grey')
names(myColors) <- levels(twins_agg[, mut_class])

for (chr in Chrom){
  dt = twins_agg[Chrom == chr]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = tumour_vaf, col = mut_class))+
    geom_point(size=2.5)+
    theme_classic(base_size = 10)+
    scale_color_manual(values = myColors)+
    labs(x = 'Genomic position', y = 'VAF (agg tumour)', col = 'Mutation category')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241109_dist_across_genome_1002_tumour_{chr}_cat.pdf'), height=3, width=5)
}

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
pdf('Results/20241109_p3_hist_vaf_dist_630muts.pdf')
hist(twins_tumour_agg[mut_ID %in% muts_mapped_tumour, tumour_vaf], breaks = 20,
     xlab = 'VAF (total tumour)', main = '630 mutations (mapped all tumour samples)', xlim = c(0, 1))
abline(v = median(twins_tumour_agg[mut_ID %in% muts_mapped_tumour, tumour_vaf]), col = 'purple')
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
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_630mapped.pdf'))

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
ggsave(glue('Results/20241109_vaf_tumour_vs_PD62341_630mapped.pdf'))

ggplot(twins_agg_filt, aes(x=tumour_vaf, y=PD63383_vaf))+
  geom_point()+
  theme_classic(base_size=15)+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  labs(x = 'VAF (total tumour)', y = glue('VAF in PD63383'), col = 'Mutation category')+
  ggtitle(glue('630 mutations mapped in the tumour, PD63383'))
ggsave(glue('Results/20241109_vaf_tumour_vs_PD63383_630mapped.pdf'))

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
# Looking at tumour-clonal mutations

twins_filtered_vaf[sum_tumour > 1 & sum_normal == 0] # 28
twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj == 0] # 52

twins_filtered_vaf[sum_tumour > 1 & sum_normal >= 1 & sum_normal <= 3] # 117
twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj >= 1 & sum_normal_adj <= 3] # 106

# which normal samples are these mutations present in?
colSums(twins_filtered_vaf[sum_tumour > 1 & sum_normal >= 1 & sum_normal <= 3, c(samples_normal_vaf), with=FALSE] >= 0.1)
# PD62341aa = 27
# PD62341h = 75
# PD63383bb = 86

colSums(twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj >= 1 & sum_normal_adj <= 3, c(samples_normal_adj), with=FALSE] >= 0.1)
# PD62341aa = 12
# PD62341h = 52
# PD63383bb = 54

# identify mutations which are present in normal samples other than contaminated ones
contaminated_samples = c('PD62341aa', 'PD62341h', 'PD63383bb')
contaminated_samples_vaf = paste0(contaminated_samples, '_VAF')
contaminated_samples_adj = paste0(contaminated_samples, '_adj')

twins_filtered_vaf[, sum_normal_cont := rowSums(.SD >= 0.1), .SDcols = contaminated_samples_vaf]
twins_filtered_vaf[sum_tumour > 1 & sum_normal <= 3 & sum_normal == sum_normal_cont] # 124

twins_filtered_vaf_adj[, sum_normal_cont_adj := rowSums(.SD >= 0.1), .SDcols = contaminated_samples_adj]
twins_filtered_vaf_adj[sum_tumour_adj > 1 & sum_normal_adj <= 3 & sum_normal_adj == sum_normal_cont_adj] # 136

# check the other mutations (present in normal samples but not contaminated ones)
twins_filtered_vaf[sum_tumour > 1 & sum_normal <= 3 & sum_normal != sum_normal_cont, mut_ID] 
# "chr10_49977949_A_G" # poor mapping 
# "chr13_18565005_G_A" # poor mapping 
# "chr15_30811818_A_G" # poor mapping 
# "chr15_56691722_C_T" # looks real q, bb
# "chr18_71485446_G_A" # looks real q, bb, h
# "chr19_43132229_T_C" # looks real aa, PD63383u ???
# "chr19_7878006_A_G" # not real  
# "chr1_83136026_A_T" # poor mapping 
# "chr20_44114996_C_T" # looks real n, h (n after adjustment)
# "chr2_57814739_A_C" # looks real q, h (q after adjustment)
# "chr5_70784531_C_A" # poor mapping 
# "chr7_120677593_C_T" # looks real n, h, bb (n, h after adjustment)
# "chr7_139659050_G_A" # looks real q, h
# "chr8_131571989_A_T" # looks real q, h, bb (q after adjustment)
# "chr8_46678612_T_C" # looks real q, h, n (q, n after adjustment)
# "chr9_100061581_C_T" # double check, not convinced PD63383ak, w
# "chr9_100061582_G_A" # double check, not convinced PD63383ak, w
# "chr9_12824689_G_A" # looks good q, bb, h
# "chr9_41808224_G_T" # poor mapping  
# "chrX_124531709_C_T" # looks good q, h, bb

muts_tumour_spec = twins_filtered_vaf[sum_tumour > 1 & sum_normal <= 3, mut_ID] %>% unlist() # 144 (NB this is just VAF-based)

# how many tumour samples are these mutations usually in?
pdf('Results/20241109_p4_hist_nr_tumour_samples_144muts.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_spec, sum_tumour], 
     xlab = 'Number of tumour samples with mutation', main = '144 mutations (min 2 tumour, max 3 normal)')
dev.off()

pdf('Results/20241109_p4_hist_vaf_144muts.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_tumour_spec, vaf_tumour_all], breaks = 20,
     xlab = 'VAF (agg tumour)', main = '144 mutations (min 2 tumour, max 3 normal)')
dev.off()

# heatmap so it is easier to think about it 
mut_tumour_144 = twins_agg_vaf[mut_ID %in% muts_tumour_spec, c('mut_ID', samples_vaf), with=FALSE]
mut_tumour_144_mat = as.matrix(mut_tumour_144[, c(samples_vaf), with=FALSE])
rownames(mut_tumour_144_mat) = mut_tumour_144[,1] %>% unlist()  
colnames(mut_tumour_144_mat) = tstrsplit(colnames(mut_tumour_144_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_tumour_144_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_tumour_144muts.pdf')
pheatmap(mut_tumour_144_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="144 mutations: likely tumour samples only", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# more investigation of PD63383 mutations
twins_vaf_sub = twins_agg_vaf[, c('mut_ID', samples_vaf), with=FALSE]
twins_vaf_melt = melt(twins_vaf_sub, id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_vaf_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'
))]
twins_vaf_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'
))]
twins_vaf_melt[, sample_type := paste(status, twin, sep = ', ')]
                 
ggplot(twins_vaf_melt[mut_ID %in% muts_PD63383_tumour_only], aes(x = value, y = mut_ID, col = sample_type))+
  geom_point()+
  theme_classic(base_size=14)+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour, col_tumour_PD63383))+
  labs(x = 'VAF', y = glue('Mutation'), col = 'Sample type')+
  ggtitle(glue('PD63383 tumour-restricted mutations'))
ggsave(glue('Results/20241109_vaf_PD63383restircted_tumour.pdf'), height = 8, width = 6)

# investigation of other classes of mutations
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2 & sum_tumour_PD63383 == 2] # 28

# mutations present in 2 tumour samples not in PD63383 tumour 
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2 & sum_tumour_PD63383 != 2] # 13
# chr10_72139813_C_T # ak, ap looks okay
# chr11_43029745_A_T # ae, b
# chr12_19555916_C_T # ae, aj
# chr13_18565005_G_A # PD63383aq, PD62341ap
# chr13_62110424_G_A # am, ap
# chr15_30811818_A_G # ag, ap
# chr15_68707110_G_A # ak, ap
# chr1_14644486_C_A # ae, aj
# chr1_38311094_G_T # ak, ap
# chr2_231591536_G_T # ae, aj
# chr5_58182780_T_A # ak, ap
# chrX_64561555_C_T # ak, ap
# chrX_9804641_G_A # ak, ap

twins_filtered_vaf[sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0] # 30
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0] # 28
setdiff(twins_filtered_vaf[sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist(), twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour_PD63383 >= 1 & sum_tumour_PD62341 == 0, mut_ID] %>% unlist())
# "chr22_43290224_C_T" - not buying it 
# "chr5_157248612_A_G"

twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 10] # 27
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 9] # 38
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 8] # 8
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 7] # 1
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 6] # 8
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 5] # 8
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 4] # 10
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 3] # 3
twins_filtered_vaf[mut_ID %in% muts_tumour & sum_tumour == 2] # 41

# add a column for names of samples which the mutation is missing in:
twins_filtered_vaf[, missing_tumour_sample := apply(twins_filtered_vaf[, ..samples_tumour_vaf], 1, function(row) {
  cols_vaf01 = samples_tumour_vaf[which(row < 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

twins_filtered_vaf[, present_tumour_sample := apply(twins_filtered_vaf[, ..samples_tumour_vaf], 1, function(row) {
  cols_vaf01 = samples_tumour_vaf[which(row >= 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 9, missing_tumour_sample])
# 1 missing in b, 2 in ag, 35 in ak

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 8, missing_tumour_sample])
# 1 missing in ag and u, 2 missing in ak and u, 5 missing in ak and ag

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 7, missing_tumour_sample])
# missing in ak, ag, u

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 6, missing_tumour_sample])
# 1 missing in ak, ag, ae, aj
# 3 missing in ak, ag, u, ae
# 2 missing in ak, ag, u, aj
# 1 missing in aj, ap, b, PD63383ap # chr1 - loss? (mapping isn't great so not 100% convinced)
# 1 missing in am, ap, b, PD63383ap # chr19 - loss?

# mutations in 5 tumour samples or less
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 5]
# chr10_49977949_A_G # not very good 
# chr11_26060345_G_T # looks okay
# chr12_112652964_C_T # looks okay
# chr13_103318073_T_A # looks okay
# chr16_17805832_C_T # looks okay
# chr17_70958446_C_T # looks okay
# chr4_135461656_G_A # looks okay
# chr6_8310684_G_T # looks okay
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 4]
# chr12_114403258_G_A # looks okay
# chr12_82832044_C_A # looks okay
# chr18_1379467_G_T # looks okay
# chr19_7878006_A_G # looks horrible
# chr5_70784531_C_A # terrible mapping 
# chr6_168731130_C_T # looks okay
# chr7_153784080_C_G # only on poorly mapped reads
# chr8_21444725_G_A # looks okay
# chr9_100061581_C_T # looks okay but a bit funny 
# chr9_41808224_G_A # poor mapping 
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 3]
# chr13_28953051_G_T # looks okay
# chr8_141860196_G_A # looks okay
# chr9_100061582_G_A # looks okay but a bit funny 
twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2] # checked these and generally look fine

table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 5, present_tumour_sample])
table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 4, present_tumour_sample])
table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 3, present_tumour_sample])
table(twins_filtered_vaf[mut_ID %in% muts_tumour_spec & sum_tumour == 2, present_tumour_sample])

######################################################################################################
# what about mutations present in all tumour but not all normal samples?
muts_tumour_all_nt[muts_tumour_all_nt %in% muts_tumour_spec] # 24
# makes sense as the other two are on chr3 G>C and C>G

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 0] # 11
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 == 1] # 21

# likely shared with all normal samples 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 >= 5] # 31

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt & sum_normal_PD63383 %in% c(2, 3, 4)]
# chr15_49480646_T_A # looks okay (probably PD62341 specific and contamination)
# chr17_18402053_G_C # poor mapping 
# chrX_798031_G_A # looks okay 

######################################################################################################

twins_filtered_vaf[sum_normal == 1 & sum_tumour == 0] # 1
# chr19_54746517_T_C Ql not terrible, 4 reads in PD63383t and generally low everywhere, wouldn't buy it 

twins_filtered_vaf[sum_normal == 2 & sum_tumour == 0] # 4
# chr17_20476029_C_T # poor mapping 
# chr4_15905566_C_T # looks fine (but then, PD62341v and PD63383ae?)
# chr6_159851462_T_C # not real 
# chr7_110349267_G_A # mapping to several places (Blat)

twins_filtered_vaf[sum_normal == 0 & sum_tumour == 1] # 4
# chr3_189240074_C_A am looks okay
# chr3_95470911_T_A ak looks okay
# chr5_157248612_A_G  PD63383aq looks okay
# chr8_91252357_A_C ak looks okay

######################################################################################################
# check my mutation classes again

# create lists of mutations: presence defined based on MTR AND VAF 
muts = twins_filtered_vaf[, mut_ID] %>% unlist()
muts_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist())) 
muts_normal = Reduce(intersect, list(twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist())) 
muts_tumour = Reduce(intersect, list(twins_filtered_vaf[sum_tumour>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist()))
muts_normal_all = Reduce(intersect, list(twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()))
muts_tumour_all = Reduce(intersect, list(twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()))
muts_tumour_only = Reduce(intersect, list(twins_filtered_vaf[sum_normal==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==0, mut_ID] %>% unlist()))
muts_tumour_all_only = Reduce(intersect, list(muts_tumour_all, muts_tumour_only)) %>% unlist()
muts_normal_only = Reduce(intersect, list(twins_filtered_vaf[sum_tumour==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==0, mut_ID] %>% unlist()))
muts_normal_all_only = Reduce(intersect, list(muts_normal_all, muts_normal_only)) %>% unlist()
muts_normal_all_nt = setdiff(muts_normal_all, muts_all)
muts_tumour_all_nt = setdiff(muts_tumour_all, muts_all)

paste('Mutations present in all samples:', length(muts_all)) # 575
paste('Number of muts present in min 1 normal sample:', length(muts_normal)) # 963
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour)) # 981
paste('Number of muts present in all normal samples:', length(muts_normal_all)) # 632
paste('Number of muts present in all tumour samples:', length(muts_tumour_all)) # 641
paste('Number of muts present in only tumour samples:', length(muts_tumour_only)) # 31
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_only)) # 0
paste('Number of muts present in only normal samples:', length(muts_normal_only)) # 16
paste('Number of muts present in all normal samples and max 1 tumour:', length(muts_normal_all_only)) # 0 # makes sense since tumour is from normal

# ONLY VAF
muts_all_vaf = twins_filtered_vaf[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist()
muts_normal_vaf = twins_filtered_vaf[sum_normal>=1, mut_ID] %>% unlist()
muts_tumour_vaf = twins_filtered_vaf[sum_tumour>=1, mut_ID] %>% unlist()
muts_normal_all_vaf = twins_filtered_vaf[sum_normal==12, mut_ID] %>% unlist()
muts_tumour_all_vaf = twins_filtered_vaf[sum_tumour==10, mut_ID] %>% unlist()
muts_tumour_only_vaf = twins_filtered_vaf[sum_normal==0, mut_ID] %>% unlist()
muts_tumour_all_only_vaf = Reduce(intersect, list(muts_tumour_all_vaf, muts_tumour_only_vaf)) %>% unlist()
muts_normal_only_vaf = twins_filtered_vaf[sum_tumour==0, mut_ID] %>% unlist()
muts_normal_all_only_vaf = Reduce(intersect, list(muts_normal_all_vaf, muts_normal_only_vaf)) %>% unlist()
muts_normal_all_nt_vaf = setdiff(muts_normal_all_vaf, muts_all_vaf)
muts_tumour_all_nt_vaf = setdiff(muts_tumour_all_vaf, muts_all_vaf)

paste('Mutations present in all samples:', length(muts_all_vaf)) # 616
paste('Number of muts present in min 1 normal sample:', length(muts_normal_vaf)) # 971
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour_vaf)) # 983
paste('Number of muts present in all normal samples:', length(muts_normal_all_vaf)) # 683
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_vaf)) # 685
paste('Number of muts present in only normal samples:', length(muts_normal_only_vaf)) # 19
paste('Number of muts present in only tumour samples:', length(muts_tumour_only_vaf)) # 31
paste('Number of muts present in all normal samples and no tumour:', length(muts_normal_all_only_vaf)) # 0
paste('Number of muts present in all tumour samples and no normal:', length(muts_tumour_all_only_vaf)) # 0
paste('Number of muts present in all normal samples but not all tumour:', length(muts_normal_all_nt_vaf)) # 67
paste('Number of muts present in all tumour samples but not all normal:', length(muts_tumour_all_nt_vaf)) # 69

# ONLY MTR
muts_all_mtr = twins_filtered_mtr[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist()
muts_normal_mtr = twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist()
muts_tumour_mtr = twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist()
muts_normal_all_mtr = twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()
muts_tumour_all_mtr = twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()
muts_tumour_only_mtr = twins_filtered_mtr[sum_normal==0, mut_ID] %>% unlist()
muts_tumour_all_only_mtr = Reduce(intersect, list(muts_tumour_all_mtr, muts_tumour_only_mtr)) %>% unlist()
muts_normal_only_mtr = twins_filtered_mtr[sum_tumour==0, mut_ID] %>% unlist()
muts_normal_all_only_mtr = Reduce(intersect, list(muts_normal_all_mtr, muts_normal_only_mtr)) %>% unlist()
muts_normal_all_nt_mtr = setdiff(muts_normal_all_mtr, muts_all_mtr)
muts_tumour_all_nt_mtr = setdiff(muts_tumour_all_mtr, muts_all_mtr)

paste('Mutations present in all samples:', length(muts_all_mtr)) # 581
paste('Number of muts present in min 1 normal sample:', length(muts_normal_mtr)) # 963
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour_mtr)) # 984
paste('Number of muts present in all normal samples:', length(muts_normal_all_mtr)) # 637
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_mtr)) # 649
paste('Number of muts present in only tumour samples:', length(muts_tumour_only_mtr)) # 39
paste('Number of muts present in only normal samples:', length(muts_normal_only_mtr)) # 18
paste('Number of muts present in all normal samples and no tumour:', length(muts_normal_all_only_mtr)) # 0
paste('Number of muts present in all tumour samples and no normal:', length(muts_tumour_all_only_mtr)) # 1
paste('Number of muts present in all normal samples but not all tumour:', length(muts_normal_all_nt_mtr)) # 56
paste('Number of muts present in all tumour samples but not all normal:', length(muts_tumour_all_nt_mtr)) # 68

######################################################################################################
# I think out of the ~60-70 mutations that are in all tumour but not all normal, we only placed 6 on the phylogenetic tree

sum(muts_tumour_all_nt %in% muts_tumour_spec) # 24
# okay so 24 of these are classified as tumour mutations
# therefore, 24 of those are assigned to tumour-specific mutations
# we know that other 8 mutations are shared b/n PD62341 normal and tumour 

# we can try to place the other mutations on this 
muts_tumour_all_nt_unplaced = setdiff(muts_tumour_all_nt, muts_tumour_spec)
muts_tumour_all_nt_unplaced = setdiff(muts_tumour_all_nt_unplaced, muts_PD62341) 
paste('Number of yet-unplaced mutations:', length(muts_tumour_all_nt_unplaced)) # 34

# always present in 10 tumour samples

# usually present in 10-12 tumour samples 
hist(twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced, sum_normal],
     xlab = 'Nr of normal samples',  main = '34 muts (unplaced tumour not all normal)')

twins_filtered_vaf[, missing_normal_sample := apply(twins_filtered_vaf[, ..samples_normal_vaf], 1, function(row) {
  cols_vaf01 = samples_normal_vaf[which(row < 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

twins_filtered_vaf[, present_normal_sample := apply(twins_filtered_vaf[, ..samples_normal_vaf], 1, function(row) {
  cols_vaf01 = samples_normal_vaf[which(row >= 0.1)]
  paste(cols_vaf01, collapse = ', ')
})]

# present in different nr of normal samples
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==4] # 1
# chr5_44907911_G_A q, aa, bb, h (could be contamination, but also 0.03 in ad and 0.06 in n, so maybe a late one)
# looks very believable, I'd be keen to buy this 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==5] # 0 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==6] # 0 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==7] # 0 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==8] # 1
# chr17_21307477_G_A missing in PD62341ad, h, n, PD63383w (is there at VAF 0.05-0.09)
# mapping doesn't give me a lot of confidence in this one 
twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==9] # 2
# chr17_18402053_G_C missing in PD63383w, ae, PD62341n (is there at VAF > 0.08 tho) # mapping isn't excellent, plausible
# chrX_798031_G_A missing in PD63383u, bb, PD62341h (is there at VAF 0.05 - 0.09)
# the thing with this one is that it looks believable but is in all normal samples at a really low VAF?

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==10] # 3
# chr11_4220977_T_C missing in PD62341q (0.03), PD63383ak (0.04) (mut only on poorly mapped reads)
# chr19_90678_G_A missing in PD62341n (0.087), PD63383bb (0.09) - likely present in all (looks super bad)
# chr5_176771316_G_A missing in PD62341n (0), PD62341aa (0.085) - I can see reads in PD62341n on Jbrowse tho, looks okay?

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==11] # 12
# chr11_4594537_A_G # poor mapping 
# chr17_36165401_A_G # poor mapping
# chr17_36388753_C_G # poor mapping
# chr1_157363251_C_T # looks real (0.08 in PD63383t)
# chr5_176771327_T_A # looks real (0.08 in PD62341n)
# chr6_32384367_G_A # looks quite plausible (0.04 in PD62341h)
# chr7_102544208_T_C # mapping is questionable
# chr7_55402890_G_T # poor mapping 
# chr7_64978551_T_C # looks okay (0.04 in PD62341h)
# chr7_76798258_T_C # looks okay (0.0976 in PD62341ad, likely everywhere, mapping could be better)
# chr8_7418832_G_C # looks kind of okay (0.087 in PD62341q)
# chr9_42329754_C_T # very questionable mapping 

twins_filtered_vaf[mut_ID %in% muts_tumour_all_nt_unplaced & sum_normal==12]
# 15 

######################################################################################################
# how many mutations do I have with VAF 0.1 - 0.2 in agg and 0.1 - 0.2 in both twins?
muts_low_vaf_both = twins_agg_vaf[vaf_normal_all > 0.1 & vaf_normal_all < 0.2 & vaf_normal_PD62341 > 0.1 & vaf_normal_PD62341 < 0.2 & vaf_normal_PD63383 > 0.1 & vaf_normal_PD63383 < 0.2, mut_ID] %>% unlist() # 59 
# wow there are actually 59 of those!

twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==4] # 1
# chr8_124634298_G_T (useless)
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==5] # 1
# chr14_19049856_C_T (useless)
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==6] # 0
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==7] # 6
# chr17_21439338_C_T # poor mapping 
# chr1_146568235_A_G # poor mapping 
# chr5_100144922_C_T # poor mapping 
# chr6_57048152_G_A # poor mapping 
# chr7_56823245_C_G # poor mapping 
# chrX_115889128_C_A # poor mapping 
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==8]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==9]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==10]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==11]
twins_filtered_vaf[mut_ID %in% muts_low_vaf_both & sum_normal==12]
# chr11_4220746_T_C # mapping could be better?
# chr12_44087006_G_A # double check but possible insertions
# chr19_40853394_C_A # poor mapping 
# chr3_190409217_C_T # looks real! > likely comes from the cluster on chr3
# chr6_94069453_A_T # insertions

# the real stuff seems to map back to clusters seen for mutations present in all samples
# perhaps they are just at lower VAF due to sampling (dropout in some samples etc.)

######################################################################################################
# are there any drivers of clonal hem in the set of normal mutations which would explain uneven transfusion b/n twins?

table(twins_dt[mut_ID %in% muts_normal_all & Gene != '-' & Effect != 'intronic', Gene])
# AGAP6, CCL3L3, CCL4L2, CROCC, CSH2, CYP2A6, DEFB4B, DUOX1, DUSP22 ,EIF5AL1, ERVV-2, GOLGA6L9, GSTT4
# GTF3C5, HERC2, KCNJ18, KIR2DL1, KRTAP10-4, LHFPL5, LILRA2, MAGEA8, MIR3680-1, NQO1, OR2T3, OR4F4, OR4K2
# OR4M1, POMZP3, POTEH, PRAMEF18, RNU1-2, SIMC1, SLC35E2B, TCAF1, THOC3, ZDHHC11B, ZNF705G

# clonal hem driver genes (from https://pmc.ncbi.nlm.nih.gov/articles/PMC11176083/)
# ZBTB33, ZNF318, ZNF234, SPRED2, SH2B3, SRCAP, SIK3, SRSF1, CHEK2, CCDC115, CCL22, BAX, YLPM1, MYD88, MTA2, MAGEC3 and IGLL5
# DNMT3A, TET2, ASXL1, PPM1D, SRSF2, SF3B1, GNB1, IDH1, IDH2, TP53, BRCC2, GNAS, JAK2, KDM6A, CBL, PHIP
# okay so it doesn't seem I have any CH drivers in my dataset 

CH_genes = c('ZBTB33', 'ZNF318', 'ZNF234', 'SPRED2', 'SH2B3', 'SRCAP', 'SIK3', 'SRSF1', 'CHEK2', 'CCDC115', 'CCL22', 'BAX', 'YLPM1', 'MYD88', 'MTA2', 'MAGEC3',
             'IGLL5', 'DNMT3A', 'TET2', 'ASXL1', 'PPM1D', 'SRSF2', 'SF3B1', 'GNB1', 'IDH1', 'IDH2', 'TP53', 'BRCC2', 'GNAS', 'JAK2', 'KDM6A', 'CBL', 'PHIP')

sum(twins_dt[Gene != '-' & Effect != 'intronic', Gene] %>% unlist() %in% CH_genes) # 42 apparently 

twins_dt_filters[Effect != 'intronic' & Gene %in% CH_genes & sum_req_filters6 <= 1, c('mut_ID', samples_vaf, 'Effect', 'Gene',
                                                                                      'f7_likelyGermline_PD63383', 'f7_likelyGermline_PD62341',                                                                                'f7_likelyGermline_bothTwins', 'f7_likelyGermline_aggTwins'), with=FALSE]
# mutations in: SH2B3, YLPM1, IDH2, SRCAP2, CCL22, SRSF1, PPM1D, TP53, ZNF234, ASXL1, CHEK2, CCDC115, IDH1, DNMT3A, JAK2, MAGEC3
# all look like they are germline regardless for how you test for this 



######################################################################################################
# Accounting for contamination (use purity estimates)

# VAF_observed = fraction_normal * VAF_normal + fraction_tumour * VAF_tumour 
# VAF_tumour = (VAF_observed - fraction_normal * VAF_normal) / fraction_tumour
# VAF_normal = (VAF_observed - fraction_tumour * VAF_tumour) / fraction_normal

# fraction_normal = fraction of normal cells, fraction_tumour = fraction of tumour cells 
# VAF_normal = VAF in normal cells, VAF_tumour = VAF in tumour cells 

# first, I need to obtain VAF values for aggregated normal and tumour samples
twins_agg_vaf = merge(twins_filtered_vaf[, c('mut_ID', samples_vaf), with=FALSE], merge(twins_filtered_mtr[, c('mut_ID', samples_mtr), with=FALSE], twins_filtered_dep[, c('mut_ID', samples_dep), with=FALSE], by = 'mut_ID'), by = 'mut_ID')

twins_agg_vaf[, dep_normal_all := rowSums(.SD), .SDcols = samples_normal_dep]
twins_agg_vaf[, dep_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_agg_vaf[, dep_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_agg_vaf[, dep_tumour_all := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_agg_vaf[, dep_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_agg_vaf[, dep_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]

twins_agg_vaf[, mtr_normal_all := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_agg_vaf[, mtr_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_agg_vaf[, mtr_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_agg_vaf[, mtr_tumour_all := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_agg_vaf[, mtr_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_agg_vaf[, mtr_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]

twins_agg_vaf[, vaf_normal_all := mtr_normal_all / dep_normal_all]
twins_agg_vaf[, vaf_normal_PD62341 := mtr_normal_PD62341 / dep_normal_PD62341]
twins_agg_vaf[, vaf_normal_PD63383 := mtr_normal_PD63383 / dep_normal_PD63383]
twins_agg_vaf[, vaf_tumour_all := mtr_tumour_all / dep_tumour_all]
twins_agg_vaf[, vaf_tumour_PD62341 := mtr_tumour_PD62341 / dep_tumour_PD62341]
twins_agg_vaf[, vaf_tumour_PD63383 := mtr_tumour_PD63383 / dep_tumour_PD63383]

twins_filtered_vaf_adj = data.table(twins_filtered_vaf)
twins_filtered_vaf_adj = merge(twins_filtered_vaf_adj, twins_agg_vaf[, c('mut_ID', 'vaf_tumour_all', 'vaf_normal_PD62341'), with=FALSE], by = 'mut_ID')

for (sample in samples_normal){
  sample_vaf = paste0(sample, '_VAF')
  sample_diff = paste0(sample, '_adj')
  sample_purity = as.numeric(purity_dt[sample==sample_vaf, purity_est])
  twins_filtered_vaf_adj[, (sample_diff) := (get(sample_vaf) - vaf_tumour_all * (1-sample_purity)) / sample_purity]
}

for (sample in samples_tumour){
  sample_vaf = paste0(sample, '_VAF')
  sample_diff = paste0(sample, '_adj')
  sample_purity = as.numeric(purity_dt[sample==sample_vaf, purity_est])
  twins_filtered_vaf_adj[, (sample_diff) := (get(sample_vaf) - vaf_normal_PD62341 * (1-sample_purity)) / sample_purity]
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
# Basic checks: number of mutations once adjusted for contamination 

# Create lists of mutations 
muts_all_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj==12 & sum_tumour_adj==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12 & sum_tumour==10, mut_ID] %>% unlist())) 
muts_normal_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal>=1, mut_ID] %>% unlist())) 
muts_tumour_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_tumour_adj>=1, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour>=1, mut_ID] %>% unlist()))
muts_normal_all_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj==12, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==12, mut_ID] %>% unlist()))
muts_tumour_all_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_tumour_adj==10, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==10, mut_ID] %>% unlist()))
muts_tumour_only_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_normal_adj==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_normal==0, mut_ID] %>% unlist()))
muts_tumour_all_only_adj = Reduce(intersect, list(muts_tumour_all_adj, muts_tumour_only_adj)) %>% unlist()
muts_normal_only_adj = Reduce(intersect, list(twins_filtered_vaf_adj[sum_tumour_adj==0, mut_ID] %>% unlist(), twins_filtered_mtr[sum_tumour==0, mut_ID] %>% unlist()))
muts_normal_all_only_adj = Reduce(intersect, list(muts_normal_all_adj, muts_normal_only_adj)) %>% unlist()
muts_normal_all_nt_adj = setdiff(muts_normal_all_adj, muts_all_adj)
muts_tumour_all_nt_adj = setdiff(muts_tumour_all_adj, muts_all_adj)

paste('Mutations present in all samples:', length(muts_all_adj)) # 344
paste('Number of muts present in min 1 normal sample:', length(muts_normal_adj)) # 945
paste('Number of muts present in min 1 tumour sample:', length(muts_tumour_adj)) # 984
paste('Number of muts present in all normal samples:', length(muts_normal_all_adj)) # 594
paste('Number of muts present in all tumour samples:', length(muts_tumour_all_adj)) # 413
paste('Number of muts present in only tumour samples:', length(muts_tumour_only_adj)) # 33
paste('Number of muts present in all tumour samples and only tumour samples:', length(muts_tumour_all_only_adj)) # 1
paste('Number of muts present in only normal samples:', length(muts_normal_only_adj)) # 12
paste('Number of muts in all tumour but not all normal:', length(muts_tumour_all_nt_adj)) # 69
paste('Number of muts in all normal but not all tumour:', length(muts_normal_all_nt_adj)) # 250
paste('Number of muts present in all normal samples and max 1 tumour:', length(muts_normal_all_only)) # 0 # makes sense since tumour is from normal

######################################################################################################
# what is the correlation between adjusted and un-adjusted VAF?
for (sample in samples_names){
  sample_vaf = paste0(sample, '_VAF')
  sample_adj = paste0(sample, '_adj')
  vaf = twins_filtered_vaf_adj[, ..sample_vaf] %>% unlist()
  adj = twins_filtered_vaf_adj[, ..sample_adj] %>% unlist()
  dt = data.table(cbind(vaf, adj))
  ggplot(dt, aes(x = vaf, y = adj))+
    geom_point(size=1.5)+
    theme_classic(base_size = 12)+
    labs(x = 'VAF', y = 'Adjusted VAF')+
    ggtitle(glue('{sample}'))
  ggsave(glue('Results/20241109_p5_vaf_vs_adj_{sample}.pdf'), height = 3, width = 3)
}

