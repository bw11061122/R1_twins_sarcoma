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
# Note: technically for ak I am estimating 0.1 but this gives rise to ridiculous values
# changed to .3 to see how it goes and will figure it out on Monday

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
    xlim(c(0, 1))+
    ylim(c(0, 1))+
    ggtitle(glue('{sample}'))+
    coord_equal(ratio=1)
  ggsave(glue('Results/20241109_p5_vaf_vs_adj_{sample}.pdf'), height = 3, width = 3)
}

######################################################################################################
# Accounting for twin-twin transfusion (spleen samples) - I want to have quantitative estimates of how much transfer there is 

muts_PD63383 = c("chr11_34011887_C_T", "chr4_75704880_G_A", "chr6_165179306_G_A",
                 "chr13_50815806_A_G", "chr4_74625500_G_T", "chr7_73831920_C_T",
                 "chr3_77633967_C_T", "chrX_115066661_C_T")
muts_PD62341 = c("chr14_105458006_C_A", "chr17_33422229_C_A", "chr15_49480646_T_A",
                 "chr16_5479739_C_T", "chr2_95662131_G_A", "chr3_50106043_C_T",
                 "chr3_62055057_C_G", "chr3_62055077_G_C", "chr20_44114996_C_T",
                 "chr20_44114996_C_T") 
mut_early = c(muts_PD62341, muts_PD63383)

mut_PD62341_dt = twins_vaf[mut_ID %in% muts_PD62341, 1:23]
mut_PD62341_melt = melt(mut_PD62341_dt, id.vars = 'mut_ID')
mut_PD62341_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD62341_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
mut_PD62341_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD62341_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD62341_melt[, mut_ID := factor(mut_ID, levels = 
                                      c('chr16_5479739_C_T','chr15_49480646_T_A',
                                        'chr17_33422229_C_A','chr14_105458006_C_A',  
                                        'chr3_50106043_C_T', 'chr2_95662131_G_A',
                                        "chr3_62055057_C_G", "chr3_62055077_G_C", 
                                        "chr20_44114996_C_T"))]


mut_PD63383_dt = twins_vaf[mut_ID %in% muts_PD63383, 1:23]
mut_PD63383_melt = melt(mut_PD63383_dt, id.vars = 'mut_ID')
mut_PD63383_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
mut_PD63383_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'))]
mut_PD63383_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
mut_PD63383_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
mut_PD63383_melt[, mut_ID := factor(mut_ID, levels = 
                                      c("chr4_75704880_G_A", "chr11_34011887_C_T",  "chr6_165179306_G_A",
                                        "chr4_74625500_G_T", "chrX_115066661_C_T", "chr13_50815806_A_G", 
                                        "chr7_73831920_C_T", "chr3_77633967_C_T"))]

# if this occurs, the spleen in PD62341 should be more similar to PD63383 than other PD62341 tissues
# at the same time, we expect spleen in PD63383 to be more similar to PD62341 than other PD63383 tissues 
# NB this is complicated because the tumour (which arose in PD62341) contaminated the skin sample (PD63383bb)

# check PD62341v and PD63383w in PD63383-specific mutations
mut_PD63383_melt[, sample_type := as.factor(paste(
  status, twin, sep = '_'))]

mut_PD63383_melt[, sample_type2 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(status, twin, sep = '_')
))]

mut_PD63383_melt[, sample_type3 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  sample == 'PD63383bb', 'normal_PD63383_skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(status, twin, sep = '_')
))]

mut_PD62341_melt[, sample_type2 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  !sample %in% c('PD62341v', 'PD63383w'), paste(status, twin, sep = '_')
))]

mut_PD62341_melt[, sample_type3 := as.factor(fcase(
  sample == 'PD62341v', 'normal_PD62341_spleen',
  sample == 'PD63383w', 'normal_PD63383_spleen',
  sample == 'PD63383bb', 'normal_PD63383_skin',
  !sample %in% c('PD62341v', 'PD63383w', 'PD63383bb'), paste(status, twin, sep = '_')
))]

ggplot(mut_PD63383_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_PD62341, '#C42F0F', col_PD63383, '#F99A49'))+
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  ggtitle(glue('PD63383-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_p3_vaf_dist_PD63383_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

ggplot(mut_PD62341_melt[status=='normal'], aes(x=mut_ID, y=value, color=sample_type2))+
  geom_point(size=2.5, position = position_jitterdodge(0.4))+
  scale_color_manual(values = c(col_PD62341, '#C42F0F', col_PD63383, '#F99A49'))+
  theme_classic(base_size = 14)+
  labs(x = 'Mutation', y = 'VAF', col = 'Sample category')+
  ggtitle(glue('PD62341-specific mutations'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_p3_vaf_dist_PD62341_muts_samples_normal_labelspleen.pdf'), width=7, height=4.5)

# Is there a way I can quantify this 
# Maybe first, compare mean VAF for 5 non-spleen samples and spleen 
means_PD62341muts_spleen_PD62341 = mut_PD62341_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD62341 = data.table(mut_PD62341_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD62341muts_spleen_PD63383 = mut_PD62341_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD62341muts_nonspleen_PD63383 = data.table(mut_PD62341_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD62341muts = merge(merge(means_PD62341muts_spleen_PD62341, means_PD62341muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD62341muts_spleen_PD63383, means_PD62341muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD62341muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))

# I would do the calculations for mutations where VAF in PD63383 non-spleen is 0 for now
# VAF_PD62341_spleen = VAF_PD62341_spleen * fraction_PD62341_spleen + VAF_PD63383 * fraction_PD63383 
# VAF_PD63383_spleen = VAF_PD62341 * fraction_PD62341 + VAF_PD63383_spleen * fraction_PD63383_spleen
means_PD62341muts_clean = means_PD62341muts[PD63383_nonspleen==0]
means_PD62341muts_clean[, PD62341_spleen_fPD62341 := PD62341_spleen / PD62341_nonspleen]
means_PD62341muts_clean[, PD63383_spleen_fPD62341 := PD63383_spleen / PD62341_nonspleen]

ggplot(means_PD62341muts_clean, aes(x=PD62341_nonspleen, y=PD62341_spleen))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'mean VAF PD62341 (non-spleen)', y = 'VAF PD62341 (spleen)')+
  ggtitle(glue('PD62341-specific mutations'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.5))+
  ylim(c(0, 0.5))
ggsave(glue('Results/20241109_spleen_vs_nonspleen_PD62341.pdf'), width=4, height=4)

# do the same for PD63383 specific mutations
means_PD63383muts_spleen_PD62341 = mut_PD63383_melt[sample == 'PD62341v', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD62341 = data.table(mut_PD63383_melt[sample %in% c('PD62341q','PD62341h', 'PD62341n', 'PD62341aa', 'PD62341ad'), mean(value), by = 'mut_ID'])
means_PD63383muts_spleen_PD63383 = mut_PD63383_melt[sample == 'PD63383w', c('mut_ID', 'value'), with=FALSE]
means_PD63383muts_nonspleen_PD63383 = data.table(mut_PD63383_melt[sample %in% c('PD63383t', 'PD63383u', 'PD63383ak','PD63383ae'), mean(value), by = 'mut_ID']) # ignore skin as this is contaminated
means_PD63383muts = merge(merge(means_PD63383muts_spleen_PD62341, means_PD63383muts_nonspleen_PD62341, by = 'mut_ID'),
                          merge(means_PD63383muts_spleen_PD63383, means_PD63383muts_nonspleen_PD63383, by = 'mut_ID'), by = 'mut_ID')
setnames(means_PD63383muts, c('value.x', 'V1.x', 'value.y', 'V1.y'), c('PD62341_spleen', 'PD62341_nonspleen', 'PD63383_spleen', 'PD63383_nonspleen'))
means_PD63383muts[, PD62341_spleen_fPD63383 := PD62341_spleen / PD63383_nonspleen] # contamination
means_PD63383muts[, PD63383_spleen_fPD63383 := PD63383_spleen / PD63383_nonspleen] # purity 

######################################################################################################
# CLASSES OF MUTATIONS 1: MUTATIONS PRESENT IN ALL NORMAL SAMPLES 
# Inspect mutations that are present in all samples
# Why were these not called as germline?

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_all] %>% unlist(), 
     xlab = 'VAF (agg normal samples)', main = '632 muts: all normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_all]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_PD62341.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD62341] %>% unlist(), 
     xlab = 'VAF (agg normal PD62341)', main = '632 muts: PD62341 normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD62341]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_PD63383.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD63383] %>% unlist(), 
     xlab = 'VAF (agg normal PD63383)', main = '632 muts: PD63383 normal samples', xlim = c(0, 1), breaks=30)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_normal_PD63383]), col = 'purple', lwd = 2.5)
dev.off()

# conclusion: these are clearly not germline based on their VAF 
# are these artefacts or real?

# plot VAF in PD62341 vs PD63383 
# if real / informative, likely enriched just in one twin
ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_normal_PD63383))+
  geom_point(size=1.5, alpha = 0.6)+
  theme_classic(base_size = 12)+
  labs(x = 'VAF (agg PD62341 normal)', y = 'VAF (agg PD63383 normal)')+
  xlim(c(0, 1))+
  ylim(c(0, 1))+
  ggtitle(glue('632 mutations (all normal samples)'))+
  coord_equal(ratio=1)
ggsave(glue('Results/20241109_vaf_normal_PD62341_vs_PD63383.pdf'), heigh = 4, width = 4)

# identify enriched mutations 

# log fold change 
twins_agg_vaf[, vaf_ratio := vaf_normal_PD62341 / vaf_normal_PD63383]
twins_agg_vaf[, vaf_ratio_log := log2(vaf_ratio)]

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_ratio_PD62341toPD63383.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio] %>% unlist(), 
     xlab = 'VAF ratio PD62341 / PD63383', main = '632 muts: VAF ratio', breaks=30)
dev.off()

pdf('Results/20241109_p4_hist_vaf_632muts_in_all_normal_samples_ratio_PD62341toPD63383_log.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_log] %>% unlist(), 
     xlab = 'log2(VAF ratio PD62341 / PD63383)', main = '632 muts: log2(VAF ratio)', breaks=30)
abline(v = -0.5, col = 'red', lwd = 2.5)
abline(v = 0.5, col = 'red', lwd = 2.5)
dev.off()

# get a p-value for comparison (KS test)
# Note: I think it doesn't hurt to do this test, BUT importantly with 6 samples per group you have very little power 
ks_test <- function(a, b) {ks.test(a, b, alternative = c("two.sided"))$p.value}
twins_agg_vaf[, ks_vaf := apply(.SD, 1, function(row) {
  a = as.numeric(row[samples_normal_PD62341_vaf])
  b = as.numeric(row[samples_normal_PD63383_vaf])
  ks_test(a, b)}), .SDcols = c(samples_normal_PD62341_vaf, samples_normal_PD63383_vaf)]
twins_agg_vaf[, ks_vaf_log := -1 * log10(ks_vaf)] # -log10(p value) for better display

twins_agg_vaf[, ks_vaf_adj := p.adjust(ks_vaf, method = 'BH')] # using Benjamini-Hochberg to adjust for multiple testing 
twins_agg_vaf[, ks_vaf_adj_log := -1 * log10(ks_vaf_adj)] # -log10(p value) for better display

paste('Minimum p-value (2-sided KS test) for normal PD62341 and PD63383:', min(twins_agg_vaf[mut_ID %in% muts_normal_all, ks_vaf])) # 0.0021645023
twins_agg_vaf[mut_ID %in% muts_normal_all & ks_vaf < 0.05, mut_ID] # 20
# "chr10_47676090_A_G" 
# "chr12_31851985_G_A" 
# "chr12_8792660_G_A"  
# "chr15_20220521_G_A"
# "chr15_20262122_C_T" 
# "chr15_22581980_C_A" 
# "chr15_81721403_G_A" 
# "chr15_81761276_G_A"
# "chr15_82483067_G_A" 
# "chr17_22069771_A_G" 
# "chr17_63871774_G_A" 
# "chr18_14868156_C_T"
# "chr19_53041299_G_C" 
# "chr1_16864604_G_A"  
# "chr1_16912588_C_T"  
# "chr1_38827952_C_A" 
# "chr3_198145129_C_A" 
# "chr3_75677255_A_G"  
# "chr7_143850635_G_T" 
# "chr8_123201027_T_C"
              
twins_agg_vaf[, mut_cat := factor(fcase(
  vaf_normal_PD62341 > 0.45 & vaf_normal_PD63383 > 0.45, 'putative germline', 
  vaf_ratio_log < -0.5 & ks_vaf < 0.05, 'PD63383 enriched',
  vaf_ratio_log > 0.5 & ks_vaf < 0.05, 'PD62341 enriched',
  ks_vaf > 0.05 | (vaf_ratio_log > -0.5 & vaf_ratio_log < 0.5), 'no difference'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_cat), aes(x = vaf_ratio_log, y = ks_vaf_log, col = mut_cat))+
  geom_point(size=1.5)+
  theme_classic(base_size = 12)+
  labs(x = 'log2fc (VAF PD62341 / VAF PD63383)', y = '-log10(p value), 2 sided KS test', col = 'Mutation category')+
  scale_color_manual(values = c('grey', col_PD62341, col_PD63383, 'darkred'))+
  ggtitle(glue('632 mutations (all normal samples)'))
ggsave(glue('Results/20241109_volcano_normal_PD62341_PD63383_col_cat.pdf'), width = 6, height = 4)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_normal_PD63383, col = mut_cat))+
  geom_point(size=1.5, alpha = 1)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF agg PD62341 normal', y = 'VAF agg PD63383 normal', col = 'Mutation category')+
  scale_color_manual(values = c('grey', col_PD62341, col_PD63383, 'darkred'))+
  ggtitle(glue('632 mutations (all normal samples)'))
ggsave(glue('Results/20241109_vaf_normal_PD62341_PD63383_col_cat.pdf'), width = 6, height = 4)

# inspect those mutations in detail 
twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat == 'putative germline'] # 5
# "chr12_131428827_C_A" # mapping is funny (checked on blat - not great)
# "chr15_85186362_A_G"  # again, mapping is funny
# "chr17_21783609_C_A" # mapping issues 
# "chr5_127735154_A_G" # mapping issues 
# "chr7_4912517_G_A" # mapping issues 
twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat == 'PD62341 enriched'] # 3
# "chr15_82483067_G_A" # kind of okay? double check
# "chr1_16864604_G_A"  # poor mapping (reads map somewhere else)
# "chr8_123201027_T_C" # does indeed look real
twins_agg_vaf[mut_ID %in% muts_normal_all & mut_cat == 'PD63383 enriched'] # 1
# "chr1_38827952_C_A" # looks real 

# plot the distribution of VAFs across samples for the good looking mutations
twins_agg_vaf_sub = twins_agg_vaf[, c('mut_ID', samples_normal_vaf), with=FALSE]
twins_agg_vaf_melt = melt(twins_agg_vaf_sub, id.vars = 'mut_ID')
twins_agg_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_agg_vaf_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341n', 'PD62341v', 'PD62341q', 'PD62341h',
  'PD63383ae', 'PD63383ak', 'PD63383u', 'PD63383bb', 'PD63383w', 'PD63383t'))]
twins_agg_vaf_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]

ggplot(twins_agg_vaf_melt[mut_ID=='chr8_123201027_T_C'], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutation: chr8:123201027, T>C'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr8_PD62341_enriched.pdf'), height = 3, width = 4.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr15_82483067_G_A'], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutation: chr15:82483067, G>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr15_PD62341_enriched.pdf'), height = 3, width = 4.5)

ggplot(twins_agg_vaf_melt[mut_ID=='chr1_38827952_C_A'], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutation: chr1:38827952, C>A'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(glue('Results/20241109_dist_samples_mut_chr1_PD63383_enriched.pdf'), height = 3, width = 4.5)

# adjusted VAF is nice but not exactly great so let's check we get a similar result with aggregate VAF (non-corrected)
purity_dt[, twin := tstrsplit(sample, '_VAF', fixed = TRUE, keep = 1)]
purity_dt[, twin := factor(fcase(twin %in% samples_PD62341, 'PD62341', twin %in% samples_PD63383, 'PD63383'))]
purity_normal_PD62341 = mean(purity_dt[status=='normal' & twin == 'PD62341', purity_est])
purity_normal_PD63383 = mean(purity_dt[status=='normal' & twin == 'PD63383', purity_est])
purity_tumour_PD62341 = mean(purity_dt[status=='tumour' & twin == 'PD62341', purity_est])
purity_tumour_PD63383 = mean(purity_dt[status=='tumour' & twin == 'PD63383', purity_est])
# not sure if I should even adjust this?

# Are all these mutations artefacts?
# Found mutations on chr12 that all look great - what is going on?
mut_chr12_cluster = c("chr12_31851749_C_T",  "chr12_31851985_G_A",  "chr12_31895265_G_A",
                      "chr12_31852349_C_G",  "chr12_31852949_G_A",  "chr12_31856122_C_T", 
                      "chr12_31867266_G_A",  "chr12_31869715_A_G",  "chr12_31872355_C_T", 
                      "chr12_31879431_C_T",  "chr12_31880061_G_A",  "chr12_31883693_C_T", 
                      "chr12_31885515_G_A",  "chr12_31885689_A_G",  "chr12_31887208_C_T" )

ggplot(twins_agg_vaf_melt[mut_ID %in% mut_chr12_cluster], aes(x = sample, y = value, col = twin))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Sample', y = 'VAF', col = 'Twin')+
  ggtitle(glue('Mutation on chr12 31,850,000-31,900,000'))+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylim(c(0, 0.6))
ggsave(glue('Results/20241109_dist_samples_mut_chr12.pdf'), width = 4.5, height = 3)

pdf('Results/20241109_p5_hist_chr12_PD62341_PD63383.pdf')
par(mfrow = c(2, 1))
hist(twins_agg[mut_ID %in% mut_chr12_cluster, PD63383_normal_vaf], xlab = 'VAF', xlim = c(0, 1), main = 'PD63383 chr12', breaks = 6)
hist(twins_agg[mut_ID %in% mut_chr12_cluster, PD62341_normal_vaf], xlab = 'VAF', xlim = c(0, 1), main = 'PD62341 chr12', breaks = 6)
dev.off()

# plot the coverage vs VAF for these mutations
twins_agg_vaf[, mut_chr12 := factor(fcase(
  mut_ID %in% mut_chr12_cluster, 'segment chr12',
  !mut_ID %in% mut_chr12_cluster, 'other'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_chr12), aes(x = dep_normal_PD62341, y = vaf_normal_PD62341, col = mut_chr12))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred'))+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD62341 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD62341.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_chr12), aes(x = dep_normal_PD63383, y = vaf_normal_PD63383, col = mut_chr12))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred'))+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD63383 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD63383.pdf'), height=5, width=7)

# similar thing on chr 14?
mut_chr14_cluster = c("chr14_19725920_G_T",  "chr14_19735467_A_C",  "chr14_19770060_C_T",  "chr14_19778921_G_T", 
                      "chr14_19803896_T_A",  "chr14_19814466_C_T",  "chr14_19823326_A_T", "chr14_19880026_C_G", 
                      "chr14_19881983_C_T",  "chr14_19887911_C_T",  "chr14_19895402_G_A",  "chr14_19899009_A_G", 
                      "chr14_19899033_A_G",  "chr14_19903456_A_G")


twins_agg_vaf[, mut_cluster := factor(fcase(
  mut_ID %in% mut_chr12_cluster, 'segment chr12',
  mut_ID %in% mut_chr14_cluster, 'segment chr14',
  !mut_ID %in% c(mut_chr12_cluster, mut_chr14_cluster), 'other'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_cluster), aes(x = dep_normal_PD62341, y = vaf_normal_PD62341, col = mut_cluster))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred', '#0db793'))+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD62341 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD62341_clusters.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all] %>% arrange(mut_cluster), aes(x = dep_normal_PD63383, y = vaf_normal_PD63383, col = mut_cluster))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  scale_color_manual(values = c('grey', 'darkred', '#0db793'))+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Location in the genome')+
  ggtitle(glue('PD63383 normal samples (632 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_632_PD63383_clusters.pdf'), height=5, width=7)

# can I plot the distribution of these mutations for each chromosome
twins_agg_vaf[, Chrom := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
twins_agg_vaf[, pos := tstrsplit(mut_ID, '_', fixed=TRUE, keep=2)]
twins_agg_vaf[, pos := as.numeric(pos)]

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

twins_agg_vaf = merge(twins_agg_vaf, chr_lengths_dt, by = 'Chrom')
twins_agg_vaf[, Chrom := factor(Chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                                              "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", 
                                              "chr20", "chr21", "chr22", "chrX"))]

for (chr in Chrom){
  dt = twins_agg_vaf[Chrom == chr & mut_ID %in% muts_normal_all]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = vaf_normal_PD62341))+
    geom_point(size=2.5)+
    theme_classic(base_size = 12)+
    scale_color_manual(values = '#7b13d4')+
    labs(x = 'Genomic position', y = 'VAF (PD62341 normal)')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241109_dist_across_genome_632_{chr}_PD62341.pdf'), height=3, width=4.5)
}

for (chr in Chrom){
  dt = twins_agg_vaf[Chrom == chr & mut_ID %in% muts_normal_all]
  length = as.numeric(dt[,Chrom_length] %>% unique()) 
  nr = dim(dt)[1]
  ggplot(dt, aes(x = pos, y = vaf_normal_PD63383))+
    geom_point(size=2.5)+
    theme_classic(base_size = 12)+
    scale_color_manual(values = '#7b13d4')+
    labs(x = 'Genomic position', y = 'VAF (PD63383 normal)')+
    ggtitle(glue('Chromosome: {chr}, {nr} mutations'))+
    xlim(c(0, length))+
    ylim(c(0, 0.7))
  ggsave(glue('Results/20241109_dist_across_genome_632_{chr}_PD63383.pdf'), height=3, width=4.5)
}

# compare to coverage across all mutations in the set
twins_agg_vaf[, mut_class := as.factor(fcase(
  mut_ID %in% muts_normal_all, 'all normal samples',
  mut_ID %in% muts_tumour_all_only, 'all tumour samples only',
  mut_ID %in% muts_tumour_all, 'all tumour samples',
  mut_ID %in% muts_tumour_only, 'only tumour samples',
  mut_ID %in% muts_normal_only, 'only normal samples',
  !mut_ID %in% c(muts_normal_all, muts_tumour_all, muts_all, muts_normal_only, muts_tumour_all_only, muts_tumour_only), 'other'))]

ggplot(twins_agg_vaf, aes(x = dep_normal_PD62341, y = vaf_normal_PD62341, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (PD62341)', y = 'VAF (PD62341 agg)', col = 'Mutation class')+
  ggtitle(glue('PD62341 normal samples (1002 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002_PD62341.pdf'), height=5, width=7)

ggplot(twins_agg_vaf, aes(x = dep_normal_PD63383, y = vaf_normal_PD63383, col = mut_class))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'Coverage (PD63383)', y = 'VAF (PD63383 agg)', col = 'Mutation class')+
  ggtitle(glue('PD63383 normal samples (1002 muts)'))+
  xlim(c(0, 500))+
  ylim(c(0, 0.7))
ggsave(glue('Results/20241109_cov_vs_vaf_mut_class_1002_PD63383.pdf'), height=5, width=7)

# Coverage histograms across different mutation classes
twins_dt[, sum_normal_PD63383_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_dt[, sum_normal_PD63383_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_dt[, sum_normal_PD62341_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_dt[, sum_normal_PD62341_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_dt[, sum_normal_PD63383_vaf := sum_normal_PD63383_mtr / sum_normal_PD63383_dep]
twins_dt[, sum_normal_PD62341_vaf := sum_normal_PD62341_mtr / sum_normal_PD62341_dep]
twins_dt[, sum_normal_mtr := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_dt[, sum_normal_dep := rowSums(.SD), .SDcols = samples_normal_dep]
twins_dt[, sum_normal_vaf := sum_normal_mtr / sum_normal_dep]

par(mfrow = c(1,1))

pdf('Results/20241109_p6_hist_PD62341_allmuts.pdf')
hist(twins_dt[sum_normal_dep < 600, sum_normal_PD62341_dep], 
     xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, all mutations (362,814)', breaks = 50)
abline(v = median(twins_dt[sum_normal_dep < 600, sum_normal_PD62341_dep]), col = 'red', lwd = 2)
dev.off()

pdf('Results/20241109_p6_hist_PD63383_allmuts.pdf')
hist(twins_dt[sum_normal_dep < 600, sum_normal_PD63383_dep], 
     xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, all mutations (362,814)', breaks = 50)
abline(v = median(twins_dt[sum_normal_dep < 600, sum_normal_PD63383_dep]), col = 'red', lwd = 2)
dev.off()

pdf('Results/20241109_p6_hist_PD62341_632muts.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all,dep_normal_PD62341] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, 632 normal')
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, dep_normal_PD62341]), col = 'blue', lwd = 2) # 264.5
dev.off()

pdf('Results/20241109_p6_hist_PD63383_632muts.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all,dep_normal_PD63383] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, 632 normal')
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all, dep_normal_PD63383]), col = 'blue', lwd = 2) # 264.5
dev.off()

pdf('Results/20241109_p6_hist_PD62341_1002muts.pdf')
hist(twins_agg_vaf[,dep_normal_PD62341] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD62341 normal)', main = 'PD62341, 1002 normal')
abline(v = median(twins_agg_vaf[,dep_normal_PD62341]), col = 'purple', lwd = 2) # 264.5
dev.off()

pdf('Results/20241109_p6_hist_PD63383_1002muts.pdf')
hist(twins_agg_vaf[,dep_normal_PD63383] %>% unlist(), breaks = 50, xlab = 'Depth (aggregated PD63383 normal)', main = 'PD63383, 1002 normal')
abline(v = median(twins_agg_vaf[,dep_normal_PD63383]), col = 'purple', lwd = 2) # 264.5
dev.off()

######################################################################################################
# Effects of those copy number changes on protein-coding genes

gene_effects_632 = twins_dt[mut_ID %in% muts_normal_all & Gene != '-', c('mut_ID', 'Gene', 'Effect')]
genes_632 = data.table(table(gene_effects_632[, 'Gene']))
genes_632_clusters = genes_632[N > 1]
# there are some genes which would have been amplified if my hypothesis is true
# I don't see anything that screams childhood cancer though 

######################################################################################################
# How do mutations present in all normal samples look like in tumour samples?

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_tumour_all))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD62341 normal)', y = 'VAF (tumour)')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD62341_normal_vs_tumour_632.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD63383, y = vaf_tumour_all))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD63383 normal)', y = 'VAF (tumour)')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD63383_normal_vs_tumour_632.pdf'), height=5, width=7)

twins_agg_vaf[, vaf_ratio_PD62341 := vaf_normal_PD62341 / vaf_tumour_all]
twins_agg_vaf[, vaf_ratio_PD62341_log := log2(vaf_ratio_PD62341)]
twins_agg_vaf[, vaf_ratio_PD63383 := vaf_normal_PD63383 / vaf_tumour_all]
twins_agg_vaf[, vaf_ratio_PD63383_log := log2(vaf_ratio_PD63383)]

pdf('Results/20241107_p1_hist_vaf_ratio_PD62341_normal_vs_tumour_632.pdf')
hist(twins_agg[mut_ID %in% muts_normal_all, vaf_ratio_PD62341] %>% unlist(), xlab = 'Ratio VAF PD62341 normal / tumour', main = '632 mutations (all normal)')
dev.off()

pdf('Results/20241107_p1_hist_vaf_ratio_PD62341_normal_vs_tumour_632_log.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_PD62341_log] %>% unlist(), xlab = 'log2(Ratio VAF PD62341 normal / tumour)', main = '632 mutations (all normal)')
dev.off()

pdf('Results/20241107_p1_hist_vaf_ratio_PD63383_normal_vs_tumour_632.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_PD63383] %>% unlist(), xlab = 'Ratio VAF PD63383 normal / tumour', main = '632 mutations (all normal)')
dev.off()

pdf('Results/20241107_p1_hist_vaf_ratio_PD63383_normal_vs_tumour_632_log.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all, vaf_ratio_PD63383_log] %>% unlist(), xlab = 'log2(Ratio VAF PD63383 normal / tumour)', main = '632 mutations (all normal)')
dev.off()

twins_agg_vaf[, mut_class_n2t_ratio_PD62341 := factor(fcase(
  vaf_ratio_PD62341_log < -0.5, 'enriched in tumour',
  vaf_ratio_PD62341_log > 0.5, 'enriched in normal (PD62341)',
  vaf_ratio_PD62341_log >= -0.5 & vaf_ratio_PD62341_log <= 0.5, 'no difference'))]

twins_agg_vaf[, mut_class_n2t_ratio_PD63383 := factor(fcase(
  vaf_ratio_PD63383_log < -0.5, 'enriched in tumour',
  vaf_ratio_PD63383_log > 0.5, 'enriched in normal (PD63383)',
  vaf_ratio_PD63383_log >= -0.5 & vaf_ratio_PD63383_log <= 0.5, 'no difference'))]

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD62341, y = vaf_tumour_all, col = mut_class_n2t_ratio_PD62341))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD62341 normal)', y = 'VAF (tumour)', col = 'Mutation category')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  scale_color_manual(values = c(col_PD62341, col_tumour, 'grey'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD62341_normal_vs_tumour_632_col.pdf'), height=5, width=7)

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all], aes(x = vaf_normal_PD63383, y = vaf_tumour_all, col = mut_class_n2t_ratio_PD63383))+
  geom_point(size=2.5, alpha = 0.6)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (PD63383 normal)', y = 'VAF (tumour)', col = 'Mutation category')+
  ggtitle(glue('Mutations present in all normal samples (632)'))+
  scale_color_manual(values = c(col_PD63383, col_tumour, 'grey'))+
  xlim(c(0, 0.7))+
  ylim(c(0, 0.7))+
  coord_equal(ratio = 1)
ggsave(glue('Results/20241109_vaf_PD63383_normal_vs_tumour_632_col.pdf'), height=5, width=7)

twins_agg_vaf[mut_ID %in% muts_normal_all & mut_class_n2t_ratio_PD62341 == 'enriched in tumour']
# chr1_16545599_G_C
# chr15_27127735_T_G
# chr4_189788936_T_C 

twins_agg_vaf[mut_ID %in% muts_normal_all & mut_class_n2t_ratio_PD63383 == 'enriched in tumour']
# chr1_16545599_G_C
# chr1_16864604_G_A
# chr1_86000_A_C

twins_agg_vaf[mut_ID %in% muts_normal_all & mut_class_n2t_ratio_PD63383 == 'enriched in normal (PD63383)' & mut_class_n2t_ratio_PD62341 == 'enriched in normal (PD62341)']
# "chr1_103717717_G_A" 
# "chr1_13225673_C_T"  
# "chr1_13225701_G_A"  
# "chr1_1646526_A_C"  
# "chr1_16652002_G_A"  
# "chr1_16772756_C_T"  
# "chr1_16889637_G_A"  
# "chr1_16892573_C_T" 
# "chr1_16905610_C_G"  
# "chr1_1697008_G_T"   
# "chr1_25409654_T_A"  
# "chr1_47072883_A_T" 
# "chr1_77402385_G_A"  
# "chr13_69946674_A_G" 
# "chr13_69946688_A_G" 
# "chr13_69946752_A_G"
# "chr18_44721988_C_T" 
# "chr18_55923631_A_G" 
# "chr18_56506391_G_C" 
# "chr9_39553748_C_T"

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
# CLASSES OF MUTATIONS 2: present in all normal but not all tumour samples

paste('Number of mutations present in all normal but not all tumour samples:', length(muts_normal_all_nt))
twins_agg_vaf[mut_ID %in% muts_normal_all_nt]

# Are the sets for un-adjusted and adjusted VAF very different?
setdiff(muts_normal_all_nt, muts_normal_all_nt_adj) # 11
setdiff(muts_normal_all_nt_adj, muts_normal_all_nt) # 21
Reduce(intersect, list(muts_normal_all_nt, muts_normal_all_nt_adj)) # 46

# plot a heatmap just so you know what's going on
mut_normal_nt = twins_filtered_vaf_adj[mut_ID %in% muts_normal_all_nt, c('mut_ID', samples_vaf), with=FALSE]
mut_normal_nt_mat = as.matrix(mut_normal_nt[, c(samples_vaf), with=FALSE])
rownames(mut_normal_nt_mat) = mut_normal_nt[,1] %>% unlist()  
colnames(mut_normal_nt_mat) = tstrsplit(colnames(mut_normal_nt_mat), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame(Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')), 
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_normal_nt_mat)
annotation_colors = list(Status = c(normal=col_normal, tumour=col_tumour), Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))

# heatmap
pdf('Results/20241109_p6_heatmap_all_normal_some_tumour.pdf')
pheatmap(mut_normal_nt_mat,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="57 muts: all normal, but not all tumour", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

# heatmap on VAF-adjusted values
mut_normal_nt_adj = twins_filtered_vaf_adj[mut_ID %in% muts_normal_all_nt, c('mut_ID', samples_normal_adj, samples_tumour_adj), with=FALSE]
mut_normal_nt_mat_adj = as.matrix(mut_normal_nt_adj[, c(samples_normal_adj, samples_tumour_adj), with=FALSE])
rownames(mut_normal_nt_mat_adj) = mut_normal_nt_adj[,1] %>% unlist()  
colnames(mut_normal_nt_mat_adj) = tstrsplit(colnames(mut_normal_nt_mat_adj), '_adj', fixed=TRUE, keep=1) %>% unlist()

pheatmap(mut_normal_nt_mat_adj,
         cellwidth=10, cellheight=2,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="57 muts: all normal, but not all tumour", 
         legend = T,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=11, cexCol=2) 

# how many tumour samples is the mutation usually absent from
pdf('Results/20241109_p7_hist_muts_all_normal_some_tumour_vaf.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_normal_all_nt, sum_tumour] %>% unlist(),
     xlab = 'Number of tumour samples with mutation', main = 'Muts in all normal but not all tumour (57), VAF')
dev.off()

pdf('Results/20241109_p7_hist_muts_all_normal_some_tumour_mtr.pdf')
hist(twins_filtered_mtr[mut_ID %in% muts_normal_all_nt, sum_tumour] %>% unlist(),
     xlab = 'Number of tumour samples with mutation', main = 'Muts in all normal but not all tumour (57), MTR')
dev.off()

# what chromosomes are these mutations on? could these be lost in some tumour clones?
chr_counts_normal_nt = data.table(table(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, Chrom]))
chr_counts_normal_nt[, V1 := factor(V1, levels = c('chr1', 'chr2', 'chr3', 'chr4', 'chr5',
                                                        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                                                        'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
                                                        'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                                                        'chr21', 'chr22', 'chrX'))]
ggplot(chr_counts_normal_nt, aes(x = V1, y = N))+
  geom_bar(stat = 'identity')+
  theme_classic(base_size = 14)+
  labs(x = 'Chromosome', y = 'Frequency')+
  ggtitle(glue('Muts in normal but not all tumour samples (57)'))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('Results/20241109_p4_muts_normal_nt_dist_chr.pdf', height = 5, width = 7.5)

muts_normal_all_nt
# "chr10_79512492_C_G" # slightly odd mapping   
# "chr12_129106160_C_T" # some mapping issues 
# "chr12_55594747_G_T" # variable mapping 
# "chr16_20575297_A_G" # variable mapping   
# "chr16_88009831_A_G" # looks okay
# "chr17_18843733_C_T" # some mapping issues 
# "chr22_24001246_G_A" # ay reads in 1 direction map terribly wrong
# "chr2_1440752_T_C" # looks great and such a fun region!!!   
# "chr2_20583359_T_C" # lots of insertions so cannot trust it  
# "chr3_190409217_C_T" # very interesting mutation pattern 
# "chr5_784531_G_A" # issues with read pair mates mapping     
# "chr6_43580159_A_G" # some mapping issues   
# "chr6_94069453_A_T" # insertions so am thinking not great? 
# "chr7_126168639_G_A" # mapping is a bit funny   
# "chr7_143742625_A_G" # mapping is a bit funny but kind of okay? 
# "chr9_133062780_G_C" # mapping issues 
# "chr9_64775433_A_G" # mapping issues yet again  
# "chr9_77762599_A_T" # mapping issues   
# "chrX_140138202_C_T" # extreme strand bias in wt reads  
# "chrX_50060796_A_G" # also poor mapping   
# "chrX_65021869_T_A" # issues with mate mapping 
# "chrX_66502763_G_A" # looks great

# chromosome 18
# "chr18_44721988_C_T" - looks great  
# "chr18_55923631_A_G" - can we please have a look at this one? 
# "chr18_56506391_G_C" - can we please have a look at this one? 

# chromosome 13
# "chr13_38383628_T_C" - looks fine 
# "chr13_69946674_A_G" - segments not correctly aligned  
# "chr13_69946688_A_G" - segments not correctly aligned  
# "chr13_69946752_A_G" - segments not correctly aligned  

# these mutations likely represent cases where a segment of the chr with the mutation was lost from some tumour samples
# but then, those mutations are not germline, so are they informative for the normal phylogeny?
pdf('Results/20241109_p6_hist_57muts_vaf.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_all], xlim = c(0, 0.6),
     xlab = 'VAF (agg normal)', main = '57 mutations (all normal, some tumour)', breaks = 20)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_all]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p6_hist_57muts_vaf_PD62341.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD62341], xlim = c(0, 0.6),
     xlab = 'VAF (agg normal)', main = '57 mutations (all normal, some tumour), PD62341', breaks = 20)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD62341]), col = 'purple', lwd = 2.5)
dev.off()

pdf('Results/20241109_p6_hist_57muts_vaf_PD63383.pdf')
hist(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD63383], xlim = c(0, 0.6),
     xlab = 'VAF (agg normal)', main = '57 mutations (all normal, some tumour), PD63383', breaks = 20)
abline(v = median(twins_agg_vaf[mut_ID %in% muts_normal_all_nt, vaf_normal_PD63383]), col = 'purple', lwd = 2.5)
dev.off()

# compare mutations between PD62341 and PD63383
twins_agg_vaf[, ratio_PD62341_PD63383_normal := vaf_normal_PD62341 / vaf_normal_PD63383]
twins_agg_vaf[, ratio_PD62341_PD63383_normal_log := log2(ratio_PD62341_PD63383_normal)]
twins_agg_vaf[, mut_class_ratio_n2n := factor(fcase(
  vaf_normal_PD62341 > 0.35 & vaf_normal_PD63383 > 0.35, 'putative germline',
  ratio_PD62341_PD63383_normal_log < -0.5, 'enriched in PD63383',
  ratio_PD62341_PD63383_normal_log > 0.5, 'enriched in PD62341',
  ratio_PD62341_PD63383_normal_log > -0.5 & ratio_PD62341_PD63383_normal_log < 0.5, 'no difference'))]

myColors <- c(col_PD62341, col_PD63383, 'grey', 'darkred')
names(myColors) <- levels(twins_agg_vaf[, mut_class_ratio_n2n])

ggplot(twins_agg_vaf[mut_ID %in% muts_normal_all_nt], aes(x = vaf_normal_PD62341, y = vaf_normal_PD63383, col = mut_class_ratio_n2n))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'VAF (agg PD62341 normal)', y = 'VAF (agg PD63383 normal)', col = 'Mutation category')+
  ggtitle(glue('57 mutations (all normal, not all tumour)'))+
  coord_equal(ratio = 1)+
  xlim(c(0, 0.6))+
  ylim(c(0, 0.6))+
  scale_color_manual(values = myColors)
ggsave('Results/20241109_vaf_PD62341_vs_PD63383_normal_57muts.pdf', height = 5, width = 7.5)

######################################################################################################
# CLASSES OF MUTATIONS 3: MUTATIONS PRESENT IN ONLY NORMAL SAMPLES 

# examine the smaller sets of mutations
muts_normal_only_adj # 12
# "chr11_34011887_C_T" # looks real
# "chr13_50815806_A_G" # looks real
# "chr1_103587565_A_C" # looks good
# "chr21_40193588_G_A" # looks real
# "chr3_165901319_C_A" # good enough
# "chr3_77633967_C_T" # looks okay
# "chr4_15905566_C_T" # looks real
# "chr4_74625500_G_T" # looks real 
# "chr4_75704880_G_A" # looks good
# "chr6_159851462_T_C" # absolutely not real
# "chr7_149688370_C_T" # looks okay
# "chrX_115066661_C_T" # looks real

muts_normal_only
# "chr11_34011887_C_T" # looks real
# "chr13_50815806_A_G" # looks real
# "chr1_103587565_A_C" # looks good
# "chr21_40193588_G_A" # looks real
# "chr4_15905566_C_T" # looks real
# "chr4_74625500_G_T" # looks real 
# "chr4_75704880_G_A" # looks good
# "chr6_159851462_T_C" # absolutely not real
# "chr7_149688370_C_T" # looks okay
# "chrX_115066661_C_T" # looks real

setdiff(muts_normal_only, muts_normal_only_adj)
# found only w/o adjusting for purity
# "chr17_20476029_C_T" # maps to several places
# "chr21_28681933_A_C" # not excellent
# "chr7_110349267_G_A" # bit questionable (maps to several places)
# "chr7_73831920_C_T" # looks okay

# Can I look at the VAF of each of those mutations
twins_vaf_adj_sub = twins_filtered_vaf_adj[, c('mut_ID', samples_normal_adj, samples_tumour_vaf), with=FALSE]
twins_vaf_adj_melt = melt(twins_vaf_adj_sub, id.vars = 'mut_ID')
twins_vaf_adj_melt[, sample := tstrsplit(variable, '_', keep = 1, fixed = TRUE)]
twins_vaf_adj_melt[, status := factor(fcase(
  sample %in% samples_normal, 'normal',
  sample %in% samples_tumour, 'tumour'))]
twins_vaf_adj_melt[, twin := factor(fcase(
  sample %in% samples_PD62341, 'PD62341',
  sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_adj_melt[, sample_type := paste(status, twin, sep = ', ')]
twins_vaf_adj_melt[, sample_type1 := factor(fcase(
  status == 'normal', paste0(twin, ', normal'),
  status == 'tumour', 'tumour'))]

ggplot(twins_vaf_adj_melt[mut_ID %in% muts_normal_only_adj], aes(x = value, y = mut_ID, col = sample_type1))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only normal samples, adj VAF): 12'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples_adj.pdf', height = 5, width = 7.5)

ggplot(twins_vaf_adj_melt[mut_ID %in% muts_normal_only], aes(x = value, y = mut_ID, col = sample_type1))+
  geom_point(size=2.5)+
  theme_classic(base_size = 14)+
  labs(x = 'Adjusted VAF', y = 'Mutation', col = 'Sample type')+
  ggtitle(glue('Mutations (only normal samples): 16'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))
ggsave('Results/20241109_p4_muts_only_normal_samples.pdf', height = 5, width = 7.5)

######################################################################################################
# CLASSES OF MUTATIONS 4: MUTATIONS PRESENT IN ALL TUMOUR BUT ONLY SOME NORMAL SAMPLES

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

######################################################################################################
# MUTATION CLASSES 5: MUTATIONS IN ONLY TUMOUR SAMPLES

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
# MUTATION CLASSES 4B: MUTATIONS ONLY IN PD63383 TUMOUR

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

muts_tumour_only

twins_filtered_vaf[sum_tumour > 0 & sum_normal == 0] # 31

twins_filtered_vaf[sum_tumour > 0 & sum_normal >= 1 & sum_normal <= 3] # 119
# which normal samples are these mutations present in?
colSums(twins_filtered_vaf[sum_tumour > 0 & sum_normal >= 1 & sum_normal <= 3, c(samples_normal_vaf), with=FALSE] >= 0.1)
# PD62341aa = 27
# PD62341h = 75
# PD63383bb = 86

colSums(twins_filtered_vaf_adj[sum_tumour_adj > 0 & sum_normal_adj >= 1 & sum_normal_adj <= 3, c(samples_normal_adj), with=FALSE] >= 0.1)
# PD62341aa = 12
# PD62341h = 53
# PD63383bb = 55

# identify mutations which are present in normal samples other than contaminated ones
contaminated_samples = c('PD62341aa', 'PD62341h', 'PD63383bb')
contaminated_samples_vaf = paste0(contaminated_samples, '_VAF')
contaminated_samples_adj = paste0(contaminated_samples, '_adj')

twins_filtered_vaf[, sum_normal_cont := rowSums(.SD >= 0.1), .SDcols = contaminated_samples_vaf]
twins_filtered_vaf[sum_tumour > 0 & sum_normal <= 3 & sum_normal == sum_normal_cont] # 128

twins_filtered_vaf_adj[, sum_normal_cont_adj := rowSums(.SD >= 0.1), .SDcols = contaminated_samples_adj]
twins_filtered_vaf_adj[sum_tumour_adj > 0 & sum_normal_adj <= 3 & sum_normal_adj == sum_normal_cont_adj] # 136

# check the other mutations (present in normal samples but not contaminated ones)
twins_filtered_vaf[sum_tumour > 0 & sum_normal <= 3 & sum_normal != sum_normal_cont, mut_ID] 
# "chr10_49977949_A_G" # poor mapping 
# "chr13_18565005_G_A" # poor mapping 
# "chr15_30811818_A_G" # poor mapping 
# "chr15_56691722_C_T" # looks real q, bb
# "chr18_71485446_G_A" # looks real q, bb, h
# "chr19_43132229_T_C" # looks real aa, PD63383u ???
# "chr19_7878006_A_G" # not real  
# "chr1_83136026_A_T" # poor mapping 
# "chr20_44114996_C_T" # looks real n, h (n after adjustment)
# "chr22_43290224_C_T" # please check this one! looks v interesting PD63383ae, ak
# "chr2_57814739_A_C" # looks real q, h (q after adjustment)
# "chr5_70784531_C_A" # poor mapping 
# "chr7_120677593_C_T" # looks real n, h, bb (n, h after adjustment)
# "chr7_139659050_G_A" # looks real q, h
# "chr7_63778593_T_G" # poor mapping 
# "chr8_131571989_A_T" # looks real q, h, bb (q after adjustment)
# "chr8_46678612_T_C" # looks real q, h, n (q, n after adjustment)
# "chr9_100061581_C_T" # double check, not convinced PD63383ak, w
# "chr9_100061582_G_A" # double check, not convinced PD63383ak, w
# "chr9_12824689_G_A" # looks good q, bb, h
# "chr9_41808224_G_T" # poor mapping  
# "chrX_124531709_C_T" # looks good q, h, bb

muts_tumour = twins_filtered_vaf[sum_tumour > 0 & sum_normal <= 3, mut_ID] %>% unlist() 

# how many tumour samples are these mutations usually in?
pdf('Results/20241109_p4_hist_nr_tumour_samples_150muts.pdf')
hist(twins_filtered_vaf[mut_ID %in% muts_tumour, sum_tumour], 
     xlab = 'Number of tumour samples with mutation', main = '150 mutations (min 1 tumour, max 3 normal)')
dev.off()



