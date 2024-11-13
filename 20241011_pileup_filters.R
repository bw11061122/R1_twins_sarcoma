# Script to analyse the pileup (run 15/10/2024)
# 2024-10-11 - 2024-10-14
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# dataframe created in the script: 20241015_pileup_checks.R (pileup run 14/10/2024-15/10/2024)

# OUTPUT:
# 1 list of mutations which pass required filters and will be used for phylogeny reconstruction
# 2 dataframe with pass/fail status for each mutation for each filter 
# 3 complete dataframe (pileup data + pass/fail status for each mutation for each filter)
# 4 trinucleotide context plots: signatures of mutations retained after employing each filter 

# NOTE ON ADDING COLUMNS WITH FILTERS 
# 0 indicates that filter is passed; 1 indicates that filter is failed (sample does not meet the criteria and should be excluded) 

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
library(VGAM)
library(grid)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# INPUT

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Print the number of unique mutations identified in at least one sample (ie how many mutations do I have?)
paste('Number of unique mutations identified across samples:', dim(twins_dt)[1]) # 362,814
paste('Number of unique mutations identified across samples (check for duplicates):', length(twins_dt[, mut_ID] %>% unlist() %>% unique())) # 362,814

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
col_bar = '#e87811'

######################################################################################################
# SAMPLES

# Examine which samples are being analysed 
samples = grep("_MTR", names(twins_dt), value = TRUE) # MTR = number of variant reads 
samples = sapply(strsplit(samples,"_"), `[`, 1)
paste('Number of samples analysed:', length(samples)) # remove the PD38is_wgs sample: not from the samples of interest

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
# LIST OF APPLIED FILTERS 
# 1 Exclude mutations mapped to chromosome Y
# 2 Exclude mutations within 30 base pairs of an indel
# 3a Exclude mutations with median coverage < 15 in normal samples
# 3b Exclude mutations with median coverage > 60 in normal samples
# 4 Require mutation to be well mapped (MTR >= 4 AND VAF >= 0.1) in at least one sample
# 5 Exclude mutations exhibiting strand bias in mutant reads but not in wild type reads across all samples
# 6 Exclude mutations with substantial poor quality mapping (ratio of reads mapped in high vs low quality pileup < 0.7)
# 7 Exclude mutations identified as germline both normal samples of PD62341 AND normal samples of PD63383 in low-quality pileup 

######################################################################################################
# 1 FILTERING OUT READS MAPPING TO CHROMOSOME Y: EXCLUDE MUTATIONS MAPPED TO Y

# Add column to indicate read mapping to chromosome Y
twins_dt[, f1_mappedY := as.numeric(Chrom=='chrY')]
paste('Number of mutations mapped to chromosome Y:', dim(twins_dt[f1_mappedY==1])[1]) # 277

######################################################################################################
# 2 FILTERING BASED ON DISTANCE TO INDELS (DATA FROM PINDEL)
# Remove substitutions present next to any indels, as these are often poor quality (indels can lead to mapping issues) 

# Read indel dataframe (coordinates of low and high quality indels)
indel_all = data.table(read.csv('Data/PD62341_PD63382_unmatched_indels_pass_or_fail.txt', sep = '\t', header = FALSE))
setnames(indel_all, c('V1', 'V2', 'V3', 'V4'), c('Chrom', 'pos', 'ref', 'alt'))
indel_all[,mut_ID := paste(Chrom, pos, ref, alt, sep = '_')]
indel_all_id = indel_all[, mut_ID] %>% unlist()
paste("Number of indels identified in the dataset:", length(indel_all_id)) # 2,728,148

# Exclude reads where reference not mapped
indel_all_dt = indel_all[!is.na(ref)]
indel_all_dt[, length := str_count(alt) - str_count(ref)]

# Add start and end coordinates 30 bp around the deletion
# Identify boundaries of regions to screen for presence of substitutions to remove 
indel_all_dt[, start30 := pos - 30]
indel_all_dt[, end30 := pos + abs(length) + 30]

indel_all_dt = indel_all_dt[, c('Chrom', 'start30', 'end30'), with=FALSE]
indel_all_dt[, chrom_start := paste(Chrom, start30, sep = '_')]
indel_all_dt[, chrom_end := paste(Chrom, end30, sep = '_')]
chrom_start = indel_all_dt[, chrom_start] %>% unique() %>% unlist()
chrom_end = indel_all_dt[, chrom_end] %>% unique() %>% unlist()

twins_coord = twins_dt[, c('mut_ID', 'Chrom', 'Pos'), with=FALSE]
twins_coord[, coord := paste(Chrom, Pos, sep = '_')]

# Identify mutations close to indels present on each chromosome 
mut_indels = c()
for (chr in twins_coord[, Chrom] %>% unlist() %>% unique()){
  tc = twins_coord[Chrom==chr]
  ind = indel_all_dt[Chrom==chr]
  starts = ind[, start30] %>% unlist() %>% unique()
  ends = ind[, end30] %>% unlist() %>% unique()
  res = tc[ind, on = .(Pos >= start30, Pos <= end30), nomatch = 0]
  mut_indels = c(mut_indels, res[, mut_ID] %>% unlist())
}

# Add column to indicate mutation being present nearby an indel 
twins_dt[, f2_FailedIndelNearby30 := as.numeric(mut_ID %in% mut_indels)]
paste('Number of mutations near poor-quality indels:', dim(twins_dt[f2_FailedIndelNearby30==1])[1]) # 13,524

######################################################################################################
# 3 FILTERING BASED ON COVERAGE 

# Exclude mutations with median coverage < 15 or > 60 in normal samples
twins_dt[, median_dep_normal := apply(.SD, 1, function(x) median(x)), .SDcols = samples_dep] # calculate median depth across all normal samples 
twins_dt[,f3_lowDepthNormal := as.numeric(median_dep_normal < 15)] # 1 if median depth is below 15 
twins_dt[,f3_highDepthNormal := as.numeric(median_dep_normal > 60)] # 1 if median depth is above 60
paste('Number of sites with median coverage below 15 (normal samples):', dim(twins_dt[f3_lowDepthNormal==1])[1]) # 9,041
paste('Number of sites with median coverage above 60 (normal samples):', dim(twins_dt[f3_highDepthNormal==1])[1]) # 527

######################################################################################################
# 4 FILTERING BASED ON EVIDENCE IN MIN 1 SAMPLE

# Require that in at least sample, the mutation is supported by good evidence (MTR >= 4 AND VAF >= 0.1 in at least one sample)

# for each sample, check if the sample has MTR >= 4 and VAF >= 0.1
# 0 if both are fulfilled, 1 if either is not 
for(sample in samples_names){
  twins_sample = twins_dt[, c('mut_ID', glue('{sample}_MTR'), glue('{sample}_VAF'))]
  setnames(twins_sample, c(glue('{sample}_MTR'), glue('{sample}_VAF')), c('mtr', 'vaf'))
  twins_sample[, glue('f4_{sample}') :=  1 - as.numeric(mtr >= 4) * as.numeric(vaf >= 0.1)]
  twins_dt = merge(twins_dt, twins_sample[, c('mut_ID', glue('f4_{sample}')), with = FALSE], by = 'mut_ID')
}

f4_samples = paste0('f4_', samples_names)
twins_dt[, f4_sum := rowSums(.SD == 0), .SDcols = f4_samples]
twins_dt[, f4_mtrAndVaf := as.numeric(f4_sum == 0)] # 0 if VAF > 0.1 and MTR > 4 in at least one sample
paste("Number of mutations which are NOT well mapped by both MTR and VAF in at least 1 sample:", dim(twins_dt[f4_mtrAndVaf==1])[1]) # 5,501

######################################################################################################
# 5 FILTERING BASED ON STRAND BIAS (IN VARIANT READS, BUT NOT WT READS TO INCREASE SPECIFICITY)

# define the function to obtain p value from 2-sided exact binomial test 
bin_test_2sided <- function(a, b, p = 0.5) {
  if (b > 0) {binom.test(a, b, 0.5, alternative = c("two.sided"), conf.level = 0.95)$p.value}
  else {2}} # if no reads are mapped at all, return 2 so these cases can be easily identified and excluded

# Subset strand information (for normal and tumour samples)
forward_cols = grep("_FAZ|_FCZ|_FTZ|_FGZ", names(twins_dt), value = TRUE) 
reverse_cols = grep("_RAZ|_RCZ|_RTZ|_RGZ", names(twins_dt), value = TRUE) 

# 1 Strand bias across all reads  
twins_dt[,sum_FAZ := rowSums(.SD), .SDcols = patterns("_FAZ", cols = forward_cols)]
twins_dt[,sum_FTZ := rowSums(.SD), .SDcols = patterns("_FTZ", cols = forward_cols)]
twins_dt[,sum_FCZ := rowSums(.SD), .SDcols = patterns("_FCZ", cols = forward_cols)]
twins_dt[,sum_FGZ := rowSums(.SD), .SDcols = patterns("_FGZ", cols = forward_cols)]

twins_dt[,sum_RAZ := rowSums(.SD), .SDcols = patterns("_RAZ", cols = reverse_cols)]
twins_dt[,sum_RTZ := rowSums(.SD), .SDcols = patterns("_RTZ", cols = reverse_cols)]
twins_dt[,sum_RCZ := rowSums(.SD), .SDcols = patterns("_RCZ", cols = reverse_cols)]
twins_dt[,sum_RGZ := rowSums(.SD), .SDcols = patterns("_RGZ", cols = reverse_cols)]

twins_dt[,sum_forward := sum_FAZ + sum_FCZ + sum_FTZ + sum_FGZ]
twins_dt[,sum_reverse := sum_RAZ + sum_RCZ + sum_RTZ + sum_RGZ]
twins_dt[,sum_all := sum_forward + sum_reverse]

# Filter out mutations where sum of forward + reverse = 0
twins_dt[, p_val := mapply(bin_test_2sided, sum_forward, sum_all)]

# correct for multiple testing (Bonferroni correction for nr of tests)
twins_dt[, significance := as.factor(fcase( 
  p_val < 0.01/dim(twins_strands)[1], 'significant', # differs from 0.5 so biased and maybe we don't want it
  p_val >= 0.01/dim(twins_strands)[1], 'not_significant'
))]

# 1 test for strand bias in mutant reads 
# For each mutation, select relevant forward and reverse reads (corresponding to the Alt variant)
twins_dt[, forward_mut := sapply(1:.N, function(row) { # identify forward strands carrying the mutation
  alt = Alt[row]
  cols_to_sum = names(twins_dt)[grepl(paste0('F', alt, 'Z'), names(twins_dt))]
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)
})]  
twins_dt[, reverse_mut := sapply(1:.N, function(row) { # identify reverse strands carrying the mutation 
  alt = Alt[row]
  cols_to_sum = names(twins_dt)[grepl(paste0('R', alt, 'Z'), names(twins_dt))]
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)
})]

# Sum mutant reads (equivalent to N = number of trials) 
twins_dt[, sum_all_mut := forward_mut + reverse_mut]

# Calculate p value from exact binomial (forward reads = number of successes, all reads = number of trials)
twins_dt[, p_val_mut := mapply(bin_test_2sided, forward_mut, sum_all_mut)]

# correct for multiple testing (Bonferroni correction for nr of tests)
twins_dt[, significance_mut := as.factor(fcase( 
  p_val_mut < 0.01/dim(twins_dt)[1], 'significant', # support for H1; p < 0.5; include 
  p_val_mut >= 0.01/dim(twins_dt)[1], 'not_significant' # no support for H1; likely germline  
))] # really want to minimize false positives because there is a risk of throwing away high-quality mutations

# Add column with strand bias filter 
twins_strands[, f5_strandBiasMutOnly := as.numeric(significance==0 & significance_mut==1)] 
paste('Number of mutations which show strand bias (mutant reads only):', dim(twins_strands[f5_strandBiasMutOnly==1])[1]) # 1,539

######################################################################################################
# FILTERING BASED ON RATIO OF READS MAPPED IN HQ AND LQ PILEUP 

# read in files with high and low quality pileup 
twins_hq = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # -mq 30   
twins_lq = data.table(read.csv('Data/pileup_merged_20241025.tsv')) # -mq 5

twins_PDv38is = grep("PDv38is", names(twins_hq), value = TRUE)
twins_hq[, c(twins_PDv38is) := NULL]

twins_mtr_hq = twins_hq[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep_hq = twins_hq[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf_hq = twins_hq[,c('mut_ID', samples_vaf), with=FALSE]
twins_mtr_lq = twins_lq[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep_lq = twins_lq[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf_lq = twins_lq[,c('mut_ID', samples_vaf), with=FALSE]

twins_mtr_hq[,sum_mtr_mapped := rowSums(.SD), .SDcols = samples_mtr]
twins_mtr_lq[,sum_mtr_mapped := rowSums(.SD), .SDcols = samples_mtr]

twins_dep_hq[,sum_dep_mapped := rowSums(.SD), .SDcols = samples_dep]
twins_dep_lq[,sum_dep_mapped := rowSums(.SD), .SDcols = samples_dep]

muts_hq = merge(twins_mtr_hq[, c('mut_ID','sum_mtr_mapped'), with=FALSE], twins_dep_hq[, c('mut_ID', 'sum_dep_mapped'), with=FALSE])
muts_lq = merge(twins_mtr_lq[, c('mut_ID','sum_mtr_mapped'), with=FALSE], twins_dep_lq[, c('mut_ID', 'sum_dep_mapped'), with=FALSE])

setnames(muts_hq, c('sum_mtr_mapped', 'sum_dep_mapped'), c('sum_mtr_hq', 'sum_dep_hq'))
setnames(muts_lq, c('sum_mtr_mapped', 'sum_dep_mapped'), c('sum_mtr_lq', 'sum_dep_lq'))

muts_hl = merge(muts_hq, muts_lq)
muts_hl[, ratio_mtr := sum_mtr_hq / sum_mtr_lq]
muts_hl[, ratio_dep := sum_dep_hq / sum_dep_lq]
muts_hl[, ratio_vaf := (sum_mtr_hq / sum_dep_hq) / (sum_mtr_lq / sum_dep_lq)]

# mutations to exclude based on HQ / LQ ratio
muts_exclude_ratio = muts_hl[ratio_mtr < 0.7, mut_ID] %>% unlist() 
paste('Number of mutations with a low HQ / LQ ratio:', length(muts_exclude_ratio)) # 19,028

# identify mutations which are identified as germline in the low-quality pileup 
twins_normal_lq = twins_lq[, c('mut_ID', samples_normal_mtr, samples_normal_dep), with=FALSE]
twins_normal_lq[, sum_MTR_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_normal_lq[, sum_DEP_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_normal_lq[, sum_MTR_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_normal_lq[, sum_DEP_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_normal_lq[, sum_MTR_all := sum_MTR_PD62341 + sum_MTR_PD63383]
twins_normal_lq[, sum_DEP_all := sum_DEP_PD62341 + sum_DEP_PD63383]

# remove rows where there is no coverage (total or in either twin)
twins_normal_lq = twins_normal_lq[(sum_DEP_all!=0 & sum_DEP_PD62341!=0 & sum_DEP_PD63383!=0)]

# stat test 
twins_normal_lq[, p_val_PD62341 := mapply(bin_test_less, sum_MTR_PD62341, sum_DEP_PD62341)]
twins_normal_lq[, p_val_PD63383 := mapply(bin_test_less, sum_MTR_PD63383, sum_DEP_PD63383)]
twins_normal_lq[, p_val_all := mapply(bin_test_less, sum_MTR_all, sum_DEP_all)]

# Add column as a filter (NOTE: H0 is germline so if p > 0.01 you don't reject H0)
twins_normal_lq[, f7_likelyGermline_PD62341 := as.numeric(p_val_PD62341 >= 0.01/dim(twins_normal)[1])]
twins_normal_lq[, f7_likelyGermline_PD63383 := as.numeric(p_val_PD63383 >= 0.01/dim(twins_normal)[1])] 
twins_normal_lq[, f7_likelyGermline_aggTwins := as.numeric(p_val_all >= 0.01/dim(twins_normal)[1])]  
twins_normal_lq[, f7_likelyGermline_bothTwins := as.numeric(f7_likelyGermline_PD62341==1 & f7_likelyGermline_PD63383==1)]  
paste('Number of likely germline mutations (PD62341):', dim(twins_normal_lq[f7_likelyGermline_PD62341==1])[1]) # 352,719
paste('Number of likely germline mutations (PD63383):', dim(twins_normal_lq[f7_likelyGermline_PD63383==1])[1]) # 352,860
paste('Number of likely germline mutations (agg twins):', dim(twins_normal_lq[f7_likelyGermline_aggTwins==1])[1]) # 350,968
paste('Number of likely germline mutations (PD62341 + PD63383):', dim(twins_normal_lq[f7_likelyGermline_bothTwins==1])[1]) # 352,191

# write IDs of likely germline mutations to csv
muts_exclude_lqgermline = twins_normal_lq[f7_likelyGermline_bothTwins==1, mut_ID] %>% unlist
paste("Number of mutations identified as germline in LQ pileup:", length(muts_exclude_lqgermline)) # 352,191
length(setdiff(germline_mutations, muts_exclude_lqgermline)) # in germline but not LQ pileup excluded # 1762
length(setdiff(muts_exclude_lqgermline, germline_mutations)) # in LQ pileup but not identified as putative germline # 844 

# Filtering: remove mutations mapped in all LQ samples (but this could get rid of mutations clonal in 1 twin and subclonal in the other)
twins_mtr_lq[, sum_mtr_samples := rowSums(.SD >= 4), .SDcols = samples_mtr]
muts_exclude_lqpresent = twins_mtr_lq[sum_mtr_samples==22, mut_ID] %>% unlist
paste("Number of mutations which are present in all samples in LQ pileup:", length(muts_exclude_lqpresent)) # 328,308
length(setdiff(germline_mutations, muts_exclude_lqpresent)) # in germline but not LQ pileup excluded # 26,881
length(setdiff(muts_exclude_lqpresent, germline_mutations)) # 2080

######################################################################################################
# CREATE A DATAFRAME WHICH HAS ALL FILTERS TOGETHER

# Extract dataframes with filters and merge with the main dataframe
strand_filters = twins_strands[, c('mut_ID', 'f8_strandBias', 'f8_strandBiasMut', 'f8_strandBiasMutOnly'), with = FALSE]
twins_dt_filters = merge(twins_dt, strand_filters, by = 'mut_ID')

# add filter from beta-binomial (rho)
twins_dt_filters[, f6_lowQualRatio := as.numeric(mut_ID %in% muts_exclude_ratio)]
twins_dt_filters[, f7_lowQualGermline := as.numeric(mut_ID %in% muts_exclude_lqgermline)]

######################################################################################################
# CHECK NR OF MUTATIONS EXCLUDED IN EACH FILTER

paste('Number of all mutations examined:', dim(twins_dt_filters)[1]) # 362,046
# Note that some mutations have been excluded due to depth = 0 for all or some samples 

dim(twins_dt_filters[f1_mappedY==1])[1] # 147
dim(twins_dt_filters[f2_FailedIndelNearby30==1])[1] # 13,464
dim(twins_dt_filters[f3_lowDepthNormal==1])[1] # 8,733
dim(twins_dt_filters[f3_highDepthNormal==1])[1] # 500
dim(twins_dt_filters[f4_mtrAndVaf==1])[1] # 4,739
dim(twins_dt_filters[f5_strandBiasMutOnly==1])[1] # 1,539
dim(twins_dt_filters[f6_lowQualRatio==1])[1] # 13,619
dim(twins_dt_filters[f7_lowQualGermline==1])[1] # 351,714

######################################################################################################
# Specify required filters 
columns_req_filters = c('f1_mappedY', 'f2_FailedIndelNearby30','f3_lowDepthNormal', 'f3_highDepthNormal', 
                        'f4_mtrAndVaf', 'f5_strandBiasMutOnly', 'f6_lowQualRatio', 'f7_lowQualGermline')
twins_dt_filters[, sum_req_filters := rowSums(.SD), .SDcols = columns_req_filters] 
muts_included = twins_dt_filters[sum_req_filters==0, mut_ID] %>% unlist()
paste('Number of mutations that pass required filters (1):', dim(twins_dt_filters[sum_req_filters==0])[1]) # 1002  

######################################################################################################
# OUTPUT 1: SAVE FILTERED MUTATIONS TO A TXT FILE
paste('Number of mutations that passed required filters:', length(mut_included)) # 771
write.table(mut_included, 'Data/mutations_include_20241112_1002.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
# OUTPUT 2: SAVE THE DATAFRAME TO A FILE (INCUDING ALL CLASSES OF FILTERS)
write.csv(twins_dt_filters, 'Data/twins_dt_filters_20241112_1002.csv', quote = FALSE, row.names = F)

######################################################################################################
# OUTPUT 3: SAVE THE DATAFRAME TO A FILE (MUTATION + FILTERING STATUS ONLY)
twins_dt_filters_info = twins_dt_filters[, c('mut_ID', columns_required_filters, 'sum_req_filters'), with=FALSE]
write.csv(twins_dt_filters, 'Data/twins_status_filters_20241112_1002.csv', quote = FALSE, row.names = F)

######################################################################################################
# Plot distribution of mutational classes on full and filtered dataset
muts_dt = twins_dt_filters[, c('mut_ID', 'Ref', 'Alt', 'sum_req_filters'), with=FALSE]
muts_dt[, status := as.factor(fcase(
  sum_req_filters == 0, 'passed',
  sum_req_filters > 0, 'failed'
))]
muts_dt[, mut_type := paste(Ref, Alt, sep='>')]
muts_dt[, mut_type := as.factor(fcase(
  mut_type == 'A>T', 'T>A',
  mut_type == 'A>C', 'T>G',
  mut_type == 'G>C', 'C>G',
  mut_type == 'G>T', 'C>A',
  mut_type == 'A>G', 'T>C',
  mut_type == 'G>A', 'C>T'))]

mut_counts = data.table(sort(table(muts_dt[,mut_type])), decreasing=TRUE)
mut_counts[,V1:=as.factor(V1)]

# plot distribution of mutations
ggplot(data=mut_counts, aes(x=fct_inorder(V1), y=N)) +
  geom_bar(stat='identity', fill = col_bar)+
  theme_classic(base_size = 12)+
  labs(x = 'Mutation type', y = 'Frequency', title = 'Mutation types in all samples')

# plot distribution of mutations of failed and passed mutations
mut_counts_passed = data.table(sort(table(muts_dt[status == 'passed',mut_type])))
mut_counts_failed = data.table(sort(table(muts_dt[status == 'failed',mut_type])))

mut_counts_passed[, freq_norm := N/ length(mut_included)]
mut_counts_failed[, freq_norm := N/ (386412 - length(mut_included))]

mut_counts_passed[, QC := 'passed']
mut_counts_failed[, QC := 'failed']
mut_counts_qc = data.table(rbind(mut_counts_passed, mut_counts_failed))
mut_counts_qc[, V1 := as.factor(V1)]

ggplot(data=mut_counts_qc, aes(x=fct_inorder(V1), y=freq_norm, fill=QC)) +
  geom_bar(stat='identity', position = 'dodge')+
  scale_fill_manual(values = c('darkblue','lightblue'))+
  theme_bw(base_size = 12)+
  labs(x = 'Mutation type', y = 'Fraction of all mutations in the category', title = 'Mutation types in all samples')
ggsave('Results/20241105_p1_mut_filters_types.pdf', width = 6, height = 4.5)

# compare removing germline mutations (many likely will be real)
muts_dt = twins_dt_filters[f7_likelyGermline_bothTwins==0, c('mut_ID', 'Ref', 'Alt', 'sum_req_filters'), with=FALSE]
muts_dt[, status := as.factor(fcase(
  sum_req_filters == 0, 'passed',
  sum_req_filters > 0, 'failed'
))]
muts_dt[, mut_type := paste(Ref, Alt, sep='>')]
muts_dt[, mut_type := as.factor(fcase(
  mut_type == 'A>T', 'T>A',
  mut_type == 'A>C', 'T>G',
  mut_type == 'G>C', 'C>G',
  mut_type == 'G>T', 'C>A',
  mut_type == 'A>G', 'T>C',
  mut_type == 'G>A', 'C>T'))]

mut_counts = data.table(sort(table(muts_dt[,mut_type])), decreasing=TRUE)
mut_counts[,V1:=as.factor(V1)]

# plot distribution of mutations of failed and passed mutations
mut_counts_passed = data.table(sort(table(muts_dt[status == 'passed',mut_type])))
mut_counts_failed = data.table(sort(table(muts_dt[status == 'failed',mut_type])))

mut_counts_passed[, freq_norm := N/ length(mut_included)]
mut_counts_failed[, freq_norm := N/ (dim(muts_dt)[1] - length(mut_included))]

mut_counts_passed[, QC := 'passed']
mut_counts_failed[, QC := 'failed']
mut_counts_qc = data.table(rbind(mut_counts_passed, mut_counts_failed))
mut_counts_qc[, V1 := as.factor(V1)]

ggplot(data=mut_counts_qc, aes(x=fct_inorder(V1), y=freq_norm, fill=QC)) +
  geom_bar(stat='identity', position = 'dodge')+
  scale_fill_manual(values = c('darkblue','lightblue'))+
  theme_classic(base_size = 12)+
  labs(x = 'Mutation type', y = 'Fraction of all mutations in the category', title = 'Mutation types in all samples\nexcluded putative germline mutations')
ggsave('Results/20241105_p1_mut_filters_types_germline_excluded.pdf', width = 6, height = 4.5)

######################################################################################################
# Plot distribution of mutational classes in different trinucleotide contexts after applying specific filters 

# define the function to get the trinucleotide context
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

# Plot the distribution of mutation types on the unfiltered dataset (362,814 mutations)
mybed0 = twins_dt[,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins0 = get_trinucs(mybed0, BSgenome.Hsapiens.UCSC.hg38)
twins_dt$trins0=trins0

# plot the distribution of different mutations across different contexts 
mut_sign_counts = data.table(table(twins_dt[, trins0]))
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
  labs(x = 'Context', y = 'Count', title = 'All mutations (362,814)')+
  theme_minimal(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241105_p1_mut_trins_allmuts.pdf', width = 10, height = 5.5)

dtf = twins_dt_filters[sum_req_filters==0]
mnr = dim(dtf)[1]
mybed = twins_dt_filters[sum_req_filters==0,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins = get_trinucs(mybed, BSgenome.Hsapiens.UCSC.hg38)
dtf$trins=trins

# plot the distribution of different mutations across different contexts 
mut_sign_counts = data.table(table(dtf[, trins]))
setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = 'Final set (368)')+
  theme_minimal(base_size = 20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241105_p1_mut_trins_finalset.pdf', width = 10, height = 5.5)

# Create a series of plots that show what happens when you progressively apply different filters
filters  = grep('^f', names(twins_dt_filters), value = TRUE) # extract filter columns

for (f in filters){
  
  twins_dt_f = twins_dt_filters[get(f)!=1]
  mut_nr = dim(twins_dt_f)[1]
  mybed = twins_dt_filters[get(f)!=1,c('Chrom', 'Pos', 'Ref', 'Alt')]
  trins = get_trinucs(mybed, BSgenome.Hsapiens.UCSC.hg38)
  twins_dt_f$trins=trins
  
  # plot the distribution of different mutations across different contexts 
  mut_sign_counts = data.table(table(twins_dt_f[, trins]))
  setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
  mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
  mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
  mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]
  
  # aggregate by mutation class and context 
  ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
    geom_bar(stat = 'identity')+
    scale_fill_manual(values = colors_sign)+
    facet_grid(~mut_class, scales = "free_x")+
    guides(fill="none")+ # remove legend
    labs(x = 'Context', y = 'Count', title = glue('Excluded: {f}, n = {mut_nr}'))+
    theme_minimal(base_size = 15) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme( strip.background = element_blank())+ 
    theme(panel.spacing = unit(0, "lines"))+
    theme(strip.text.x = element_text(size = 13))+
    geom_hline(yintercept = 0, colour="black", size = 0.1)
  ggsave(glue('Results/20241105_p1_mut_trins_allmuts_{f}.pdf'), width = 10, height = 5.5)
}

# we want to show how signatures change after each filter is applied

# specify order of filters such that germline removed first 
filters = c("f7_likelyGermline_bothTwins", 'f4_mtr4_presentInAll', 'f5_presentAllVaf',
            "f6_mtrAndVaf", "f4_mtr4_presentInNone", "f4_mtr4_presentInOne", "f5_absentAllVaf",
            "f3_highDepthNormal", "f3_lowDepthNormal", "f9_lowQualRatio", "f9_lowQualPresent",
            "f2_FailedIndelNearby30", "f1_mappedY")

for (i in seq_along(filters)){
  
  subset_filters = filters[1:i]
  cond = paste(subset_filters, '!= 1', collapse = ' & ')
  twins_dt_f = twins_dt_filters[eval(parse(text = cond))]
  mut_nr = dim(twins_dt_f)[1]
  mybed = twins_dt_filters[eval(parse(text = cond)), c('Chrom', 'Pos', 'Ref', 'Alt')]
  trins = get_trinucs(mybed, BSgenome.Hsapiens.UCSC.hg38)
  twins_dt_f$trins=trins
  
  # plot the distribution of different mutations across different contexts 
  mut_sign_counts = data.table(table(twins_dt_f[, trins]))
  setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
  mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
  mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
  mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]
  
  # aggregate by mutation class and context 
  ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
    geom_bar(stat = 'identity')+
    scale_fill_manual(values = colors_sign)+
    facet_grid(~mut_class, scales = "free_x")+
    guides(fill="none")+ # remove legend
    labs(x = 'Context', y = 'Count', title = glue('n = {mut_nr}'))+
    theme_minimal(base_size = 15) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme( strip.background = element_blank())+ 
    theme(panel.spacing = unit(0, "lines"))+
    theme(strip.text.x = element_text(size = 13))+
    geom_hline(yintercept = 0, colour="black", size = 0.1)
  ggsave(glue('Results/20241105_p1_mut_trins_allmuts_filters_{i}.pdf'), width = 10, height = 5.5)
}

for (i in seq_along(columns_req_filters)){
  
  subset_filters = columns_req_filters[1:i]
  cond = paste(subset_filters, '!= 1', collapse = ' & ')
  twins_dt_f = twins_dt_filters[eval(parse(text = cond))]
  mut_nr = dim(twins_dt_f)[1]
  mybed = twins_dt_filters[eval(parse(text = cond)), c('Chrom', 'Pos', 'Ref', 'Alt')]
  trins = get_trinucs(mybed, BSgenome.Hsapiens.UCSC.hg38)
  twins_dt_f$trins=trins
  
  # plot the distribution of different mutations across different contexts 
  mut_sign_counts = data.table(table(twins_dt_f[, trins]))
  setnames(mut_sign_counts, c('V1', 'N'), c('trins', 'count'))
  mut_sign_counts[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
  mut_sign_counts[, mut_class := gsub("\\.", ">", mut_class)]
  mut_sign_counts[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]
  
  # aggregate by mutation class and context 
  ggplot(data=mut_sign_counts, aes(x=context, y=count, fill=mut_class)) +
    geom_bar(stat = 'identity')+
    scale_fill_manual(values = colors_sign)+
    facet_grid(~mut_class, scales = "free_x")+
    guides(fill="none")+ # remove legend
    labs(x = 'Context', y = 'Count', title = glue('n = {mut_nr}'))+
    theme_minimal(base_size = 15) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme( strip.background = element_blank())+ 
    theme(panel.spacing = unit(0, "lines"))+
    theme(strip.text.x = element_text(size = 13))+
    geom_hline(yintercept = 0, colour="black", size = 0.1)
  ggsave(glue('Results/20241105_p1_mut_trins_allmuts_req_filters_{i}.pdf'), width = 10, height = 5.5)
}

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# ALL DONE 



