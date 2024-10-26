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
twins_mtr_hq[,sum_mtr_mapped := rowSums(.SD), .SDcols = samples_mtr]

twins_dep_hq[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_dep_hq[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_dep_hq[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_dep_hq[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]
twins_dep_hq[,sum_dep_mapped := rowSums(.SD), .SDcols = samples_dep]

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
twins_mtr_lq[,sum_mtr_mapped := rowSums(.SD), .SDcols = samples_mtr]

twins_dep_lq[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_dep_lq[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_dep_lq[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_dep_lq[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]
twins_dep_lq[,sum_dep_mapped := rowSums(.SD), .SDcols = samples_dep]

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

# we want to check the ratio of mapped reads with low and high quality called mutations
# for each mutation, get a table of summed variant reads and coverage in HQ and LQ pileup
muts_hq = merge(twins_mtr_hq[, c('mut_ID','sum_mtr_mapped'), with=FALSE], twins_dep_hq[, c('mut_ID', 'sum_dep_mapped'), with=FALSE])
muts_lq = merge(twins_mtr_lq[, c('mut_ID','sum_mtr_mapped'), with=FALSE], twins_dep_lq[, c('mut_ID', 'sum_dep_mapped'), with=FALSE])

setnames(muts_hq, c('sum_mtr_mapped', 'sum_dep_mapped'), c('sum_mtr_hq', 'sum_dep_hq'))
setnames(muts_lq, c('sum_mtr_mapped', 'sum_dep_mapped'), c('sum_mtr_lq', 'sum_dep_lq'))

muts_hl = merge(muts_hq, muts_lq)
muts_hl[, ratio_mtr := sum_mtr_hq / sum_mtr_lq]
muts_hl[, ratio_dep := sum_dep_hq / sum_dep_lq]
muts_hl[, ratio_vaf := (sum_mtr_hq / sum_dep_hq) / (sum_mtr_lq / sum_dep_lq)]

hist(muts_hl[, ratio_mtr], main = 'MTR in HQ pileup / MTR in LQ pileup', xlab = 'MTR (HQ) / MTR (LQ)')
hist(muts_hl[, ratio_dep], main = 'DEP in HQ pileup / DEP in LQ pileup', xlab = 'DEP (HQ) / DEP (LQ)')
hist(muts_hl[, ratio_vaf], main = 'VAF in HQ pileup / VAF in LQ pileup', xlab = 'VAF (HQ) / VAF (LQ)')

# for how many mutations is the ratio below 0.8 for MTR?
muts_low_ratio = muts_hl[ratio_mtr < 0.8, mut_ID] %>% unlist()
paste('Nubmer of mutations with HQ / LQ MTR ratio < 0.8', length(muts_low_ratio)) # 31171

# can I check some of those on Jbrowse to see if what is picked up is low quality 
# chr10_79868989_T_C
# chr10_87249219_A_C
# chr1_16662296_G_T
# chr5_174644141_C_T
# chr8_17933926_C_G
# chr22_18803145_C_G
# chr1_26264685_C_A # good and very interesting pattern!!!
# chr2_235023547_G_A 
# chr3_29731967_G_A
# chr20_36233124_G_A
# throws away quite a lot of good stuff 

# okay, I think 0.8 is not specific enough, let's take 0.5
muts_low_ratio = muts_hl[ratio_mtr < 0.5, mut_ID] %>% unlist()
paste('Nubmer of mutations with HQ / LQ MTR ratio < 0.5', length(muts_low_ratio)) # 7871
# chr22_24649567_T_G
# chr14_45928737_G_T
# chr15_30373275_G_T
# chr10_50183083_A_G
# chr10_46340590_C_T
# chr11_89828206_C_G
# chr12_95840269_T_G
# chr7_76627756_A_G
# chr9_62121782_C_T
# chrX_74991577_G_A
# okay this seems better 

# check if mutations weve checked are of poor quality are picked up
poorly_mapped_muts <- c("chr15_23357707_C_T", "chr16_89042129_G_A", "chr2_113447401_A_T",
                        "chr2_131534951_C_T", "chr2_95748047_C_T", "chr8_133213958_T_C",
                        "chr8_33050650_A_C", "chr11_67710577_G_T", "chr15_30819247_G_A",
                        "chr12_75561694_T_C", "chr18_15378233_T_C", "chr17_18716394_G_T",
                        "chr19_39186678_T_G", "chr1_120476102_T_G", "chr1_148952577_G_A",
                        "chr2_113513495_A_G", "chr2_172537945_G_A", "chr3_195982104_C_T")
muts_low_ratio08 = muts_hl[ratio_mtr < 0.8, mut_ID] %>% unlist()
muts_low_ratio07 = muts_hl[ratio_mtr < 0.7, mut_ID] %>% unlist()
muts_low_ratio065 = muts_hl[ratio_mtr < 0.65, mut_ID] %>% unlist()
muts_low_ratio06 = muts_hl[ratio_mtr < 0.6, mut_ID] %>% unlist()
muts_low_ratio05 = muts_hl[ratio_mtr < 0.5, mut_ID] %>% unlist()

sum(poorly_mapped_muts %in% muts_low_ratio08) # 16
sum(poorly_mapped_muts %in% muts_low_ratio07) # 16
sum(poorly_mapped_muts %in% muts_low_ratio065) # 16
sum(poorly_mapped_muts %in% muts_low_ratio06) # 14
sum(poorly_mapped_muts %in% muts_low_ratio05) # 6

well_mapped_muts <- c("chr11_43029745_A_T", "chr12_19555916_C_T", "chr12_82832044_C_A",
                        "chr17_78256883_C_T", "chr22_17774668_G_A", "chr22_30553442_G_T",
                        "chr19_3961910_G_A", "chr7_115240063_T_C", "chr7_120677593_C_T",
                        "chrX_46375010_C_T", "chrX_56448796_C_A", "chr13_103318073_T_A",
                        "chr14_54164980_C_G", "chr17_70958446_C_T", "chr4_135461656_G_A",
                        "chr14_105458006_C_A", "chr3_62055057_C_G", "chr22_30553442_G_T")

sum(well_mapped_muts %in% muts_low_ratio08) # 5
sum(well_mapped_muts %in% muts_low_ratio07) # 0
sum(well_mapped_muts %in% muts_low_ratio065) # 0
sum(well_mapped_muts %in% muts_low_ratio06) # 0
sum(well_mapped_muts %in% muts_low_ratio05) # 0

# based on this, perhaps best to use threshold of 0.6 - 0.7 
muts_low_ratio = muts_hl[ratio_mtr < 0.7, mut_ID] %>% unlist()
paste('Number of mutations with HQ / LQ MTR ratio < 0.7', length(muts_low_ratio)) # 14,368

# what kinds of mutation types are these?
mut_types_lr = sapply(strsplit(muts_low_ratio, '_'), function(x) paste(x[3], x[4], sep = '>')) 
mut_types_counts = data.table(table(mut_types_lr))
setnames(mut_types_counts, 'mut_types_lr', 'mutation_class')
mut_types_counts[, mutation_class := fcase(
  mutation_class == 'A>T', 'T>A',
  mutation_class == 'A>C', 'T>G',
  mutation_class == 'G>C', 'C>G',
  mutation_class == 'G>T', 'C>A',
  mutation_class == 'A>G', 'T>C',
  mutation_class == 'G>A', 'C>T',
  mutation_class == 'C>T', 'C>T',
  mutation_class == 'C>A', 'C>A',
  mutation_class == 'C>G', 'C>G',
  mutation_class == 'T>A', 'T>A',
  mutation_class == 'T>C', 'T>C',
  mutation_class == 'T>G', 'T>G')] # don't change if a different mutation type 
mut_types_counts[, mutation_class := as.factor(mutation_class)]
ggplot(data=mut_types_counts, aes(x=mutation_class, y=N)) +
  geom_bar(stat='identity', fill = col_bar)+
  theme_bw(base_size = 12)+
  labs(x = 'Mutation type', y = 'Frequency', title = 'Mutation types in all samples')




