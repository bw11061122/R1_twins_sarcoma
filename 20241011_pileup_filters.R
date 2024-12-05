###################################################################################################################################
# SCRIPT 1

# Script to obtain set of mutations (SNVs) to reconstruct the phylogeny
# 2024-10-11 - 2024-11-14
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# dataframe created in the script: 20241015_pileup_checks.R (pileup run 14/10/2024-15/10/2024)

# OUTPUT:
# 1 list of mutations which pass required filters and will be used for phylogeny reconstruction
# 2 dataframe with pass/fail status for each mutation for each filter 
# 3 complete dataframe (pileup data + pass/fail status for each mutation for each filter)
# 4 mutation type + tri-nucleotide context plots: signatures of mutations retained after employing each filter 

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
library(grid)
library(cowplot)
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# INPUT FILES

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = fread('Data/pileup_merged_20241016.tsv')

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
# 2 Exclude mutations which are on regions deemed as copy number altered by PURPLE
# 4a Exclude mutations with median coverage < 15 in normal samples
# 4b Exclude mutations with median coverage > 60 in normal samples
# 5 Require mutation to be well mapped (MTR >= 4 AND VAF >= 0.1) in at least one sample
# 6 Exclude mutations exhibiting strand bias in mutant reads but not in wild type reads across all samples
# 7 Exclude mutations with substantial poor quality mapping (ratio of reads mapped in high vs low quality pileup < 0.7)
# 8 Exclude mutations identified as germline both normal samples of PD62341 AND normal samples of PD63383 in low-quality pileup 

######################################################################################################
# 1 FILTERING OUT READS MAPPING TO CHROMOSOME Y: EXCLUDE MUTATIONS MAPPED TO Y

# Add column to indicate read mapping to chromosome Y
twins_dt[, f1_mappedY := as.numeric(Chrom=='chrY')]
paste('Number of mutations mapped to chromosome Y:', dim(twins_dt[f1_mappedY==1])[1]) # 277

######################################################################################################
# 2 FILTERING BASED ON DISTANCE TO INDELS (DATA FROM PINDEL)
# Remove substitutions present next to any indels, as these are often poor quality (indels can lead to mapping issues) 

# Read indel dataframe (coordinates of low and high quality indels)
indel_all = fread('Data/PD62341_PD63382_unmatched_indels_pass_or_fail.txt', sep = '\t', header = FALSE)
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
paste('Number of mutations near indels:', dim(twins_dt[f2_FailedIndelNearby30==1])[1]) # 13,524

######################################################################################################
# 3 FILTERING BASED ON COVERAGE 

# Exclude mutations with median coverage < 15 or > 60 in normal samples
twins_dt[, median_dep_normal := apply(.SD, 1, function(x) median(x)), .SDcols = samples_normal_dep] # calculate median depth across all normal samples 
twins_dt[,f3_lowDepthNormal := as.numeric(median_dep_normal < 15)] # 1 if median depth is below 15 
twins_dt[,f3_highDepthNormal := as.numeric(median_dep_normal > 60)] # 1 if median depth is above 60
paste('Number of sites with median coverage below 15 (normal samples):', dim(twins_dt[f3_lowDepthNormal==1])[1]) # 9,449
paste('Number of sites with median coverage above 60 (normal samples):', dim(twins_dt[f3_highDepthNormal==1])[1]) # 500

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
# 6 FILTERING BASED ON STRAND BIAS (IN VARIANT READS, BUT NOT WT READS TO INCREASE SPECIFICITY)

# define the function to obtain p value from 2-sided exact binomial test 
bin_test_2sided <- function(a, b, p = 0.5) {
  if (b > 0) {binom.test(a, b, 0.5, alternative = c("two.sided"), conf.level = 0.95)$p.value}
  else {2}} # if no reads are mapped at all, return 2 so these cases can be easily identified and excluded

# Subset strand information (for normal and tumour samples)
forward_cols = grep("_FAZ|_FCZ|_FTZ|_FGZ", names(twins_dt), value = TRUE) 
reverse_cols = grep("_RAZ|_RCZ|_RTZ|_RGZ", names(twins_dt), value = TRUE) 

# 1 Strand bias across mutant reads 

# For each mutation, identify relevant forward and reverse reads (corresponding to the Alt variant)
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

# For each mutation, identify the number of reads on forward and reverse strands
twins_dt[, forward_wt := sapply(1:.N, function(row) { # identify forward strands with the wt allele 
  ref = Ref[row]
  cols_to_sum = names(twins_dt)[grepl(paste0('F', ref, 'Z'), names(twins_dt))]
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)
})]  
twins_dt[, reverse_wt := sapply(1:.N, function(row) { # identify reverse strands with the wt allele  
  ref = Ref[row]
  cols_to_sum = names(twins_dt)[grepl(paste0('R', ref, 'Z'), names(twins_dt))]
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)
})]

# Sum mutant reads (equivalent to N = number of trials) 
twins_dt[, sum_all_mut := forward_mut + reverse_mut] # should be identical to MTR 
twins_dt[, sum_all_wt := forward_wt + reverse_wt] 

# P-value of exact binomial test: H0: p=0.5, H1: p=/=0.5
twins_dt[, p_val_mut := mapply(bin_test_2sided, forward_mut, sum_all_mut)]

# correct for multiple testing (Bonferroni correction for nr of tests)
twins_dt[, significance_mut := as.factor(fcase( 
  p_val_mut < 0.01/dim(twins_dt)[1], 'significant', # support for H1; p < 0.5; strand bias
  p_val_mut >= 0.01/dim(twins_dt)[1], 'not_significant' # no support for H1: no strand bias identified  
))] # really want to minimize false positives because there is a risk of throwing away high-quality mutations

# 2 Strand bias across all reads  
twins_dt[, sum_FAZ := rowSums(.SD), .SDcols = patterns("_FAZ", cols = names(twins_dt))]
twins_dt[, sum_FTZ := rowSums(.SD), .SDcols = patterns("_FTZ", cols = names(twins_dt))]
twins_dt[, sum_FCZ := rowSums(.SD), .SDcols = patterns("_FCZ", cols = names(twins_dt))]
twins_dt[, sum_FGZ := rowSums(.SD), .SDcols = patterns("_FGZ", cols = names(twins_dt))]

twins_dt[, sum_RAZ := rowSums(.SD), .SDcols = patterns("_RAZ", cols = names(twins_dt))]
twins_dt[, sum_RTZ := rowSums(.SD), .SDcols = patterns("_RTZ", cols = names(twins_dt))]
twins_dt[, sum_RCZ := rowSums(.SD), .SDcols = patterns("_RCZ", cols = names(twins_dt))]
twins_dt[, sum_RGZ := rowSums(.SD), .SDcols = patterns("_RGZ", cols = names(twins_dt))]

twins_dt[, sum_forward := sum_FAZ + sum_FCZ + sum_FTZ + sum_FGZ]
twins_dt[, sum_reverse := sum_RAZ + sum_RCZ + sum_RTZ + sum_RGZ]
twins_dt[, sum_all := sum_forward + sum_reverse]

# P-value of exact binomial test: H0: p=0.5, H1: p=/=0.5
twins_dt[, p_val := mapply(bin_test_2sided, sum_forward, sum_all)]

# Correct for multiple testing (Bonferroni correction for nr of tests)
twins_dt[, significance := as.factor(fcase( 
  p_val < 0.01/dim(twins_dt)[1], 'significant', # support for H1: strand bias identified 
  p_val >= 0.01/dim(twins_dt)[1], 'not_significant' # no support for H1: no strand bias
))]

# Fisher Exact test (based no GATK https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test)
fisher_test = function(row){
  mat = matrix(c(row[1], row[2], row[3], row[4]), nrow=2)
  fisher.test(mat)$p.value
}

twins_dt[, fishers_exact := apply(.SD, 1, fisher_test), .SDcols = c('forward_mut', 'reverse_mut', 'forward_wt', 'reverse_wt')]

twins_dt[, fishers_exact_significance := as.factor(fcase( 
  fishers_exact < 0.01/dim(twins_dt)[1], 'significant', # support for H1: strand bias identified 
  fishers_exact >= 0.01/dim(twins_dt)[1], 'not_significant' # no support for H1: no strand bias
))]

# Add column with strand bias filter 
twins_dt[, f5_strandBias := as.numeric(fishers_exact_significance == 'significant')] 
paste('Number of mutations which show strand bias:', dim(twins_dt[f5_strandBias==1])[1]) # 2,132

######################################################################################################
# 6 FILTERING BASED ON RATIO OF READS MAPPED IN HQ AND LQ PILEUP 

# Read in files with high and low quality pileup 
twins_lq = fread('Data/pileup_merged_20241025.tsv') # -mq 5

twins_dt[, sum_mtr_hq := rowSums(.SD), .SDcols = samples_mtr]
twins_lq[, sum_mtr_lq := rowSums(.SD), .SDcols = samples_mtr]
twins_dt = merge(twins_dt, twins_lq[, c('mut_ID', 'sum_mtr_lq'), with=FALSE], by = 'mut_ID')
twins_dt[, ratio_mtr := sum_mtr_hq / sum_mtr_lq]

# mutations to exclude based on HQ / LQ ratio
twins_dt[, f6_lowQualRatio := as.numeric(ratio_mtr < 0.7)]
paste('Number of mutations which low HQ / LQ ratio:', dim(twins_dt[f6_lowQualRatio==1])[1]) # 14,368

######################################################################################################
# 7 FILTERING BASED ON GERMLINE (LOW QUALITY PILEUP) 

bin_test_less <- function(a, b, p = 0.5) {
  if (b > 0) {binom.test(a, b, 0.5, alternative = c("less"), conf.level = 0.95)$p.value}
  else {2}} # if no reads are mapped at all, return 2 so these cases can be easily identified and excluded

# identify mutations which are identified as germline in the low-quality pileup 
twins_lq[, sum_mtr_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_lq[, sum_dep_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_lq[, sum_mtr_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_lq[, sum_dep_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_lq[, sum_mtr_all := sum_mtr_PD62341 + sum_mtr_PD63383]
twins_lq[, sum_dep_all := sum_dep_PD62341 + sum_dep_PD63383]

# stat test 
twins_lq[, p_val_PD62341 := mapply(bin_test_less, sum_mtr_PD62341, sum_dep_PD62341)]
twins_lq[, p_val_PD63383 := mapply(bin_test_less, sum_mtr_PD63383, sum_dep_PD63383)]
twins_lq[, p_val_all := mapply(bin_test_less, sum_mtr_all, sum_dep_all)]

# Add column as a filter (NOTE: H0 is germline so if p > 0.01 you don't reject H0)
twins_dt = merge(twins_dt, twins_lq[, c('mut_ID', 'p_val_PD62341', 'p_val_PD63383', 'p_val_all'), with=FALSE], by = 'mut_ID')
twins_dt[, likelyGermline_PD62341 := as.numeric(p_val_PD62341 >= 0.01/dim(twins_dt)[1])]
twins_dt[, likelyGermline_PD63383 := as.numeric(p_val_PD63383 >= 0.01/dim(twins_dt)[1])] 
twins_dt[, likelyGermline_aggTwins := as.numeric(p_val_all >= 0.01/dim(twins_dt)[1])]  
twins_dt[, f7_likelyGermline_bothTwins := as.numeric(likelyGermline_PD62341==1 & likelyGermline_PD63383==1)]  
paste('Number of likely germline mutations (PD62341):', dim(twins_dt[likelyGermline_PD62341==1])[1]) # 352,747
paste('Number of likely germline mutations (PD63383):', dim(twins_dt[likelyGermline_PD63383==1])[1]) # 352,889
paste('Number of likely germline mutations (agg twins):', dim(twins_dt[likelyGermline_aggTwins==1])[1]) # 350,995
paste('Number of likely germline mutations (PD62341 + PD63383):', dim(twins_dt[f7_likelyGermline_bothTwins==1])[1]) # 352,219

######################################################################################################
# CHECK NR OF MUTATIONS EXCLUDED IN EACH FILTER

paste('Number of all mutations examined:', dim(twins_dt)[1]) # 362,814
# Note that some mutations have been excluded due to depth = 0 for all or some samples 

dim(twins_dt[f1_mappedY==1])[1] # 277
dim(twins_dt[f2_FailedIndelNearby30==1])[1] # 13,524
dim(twins_dt[f3_lowDepthNormal==1])[1] # 9,449
dim(twins_dt[f3_highDepthNormal==1])[1] # 500
dim(twins_dt[f4_mtrAndVaf==1])[1] # 5,501
dim(twins_dt[f5_strandBias==1])[1] # 1,539
dim(twins_dt[f6_lowQualRatio==1])[1] # 14,368
dim(twins_dt[f7_likelyGermline_bothTwins==1])[1] # 352,219

######################################################################################################
# FILTER OUT MUTATIONS THAT DO NOT PASS FILTERS 

# Specify required filters 
columns_req_filters = c('f1_mappedY', 'f2_FailedIndelNearby30', 'f3_lowDepthNormal', 'f3_highDepthNormal', 
                        'f4_mtrAndVaf', 'f5_strandBias', 'f6_lowQualRatio', 'f7_likelyGermline_bothTwins')

twins_dt[, sum_req_filters := rowSums(.SD), .SDcols = columns_req_filters] 
muts_included = twins_dt[sum_req_filters==0, mut_ID] %>% unlist()
paste('Number of mutations that pass required filters (1):', dim(twins_dt[sum_req_filters==0])[1]) # 1069  

######################################################################################################
# OUTPUT 1: SAVE FILTERED MUTATIONS TO A TXT FILE
write.table(muts_included, 'Data/mutations_include_20241204_1069.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
# OUTPUT 1: SAVE PUTATIVE GERMLINE MUTATIONS TO A TXT FILE
muts_germline = twins_dt[sum_req_filters==1 & f7_likelyGermline_bothTwins==1, mut_ID] %>% unlist() # 332,974
write.table(muts_germline, 'Data/mutations_putativeGermline_20241204.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
# OUTPUT 1: SAVE PUTATIVE GERMLINE MUTATIONS TO A TXT FILE
muts_passed_qc = c(muts_included, muts_germline)
paste('Number of included mutations:', length(muts_passed_qc)) # 334,043
write.table(muts_passed_qc, 'Data/mutations_passedQuality_20241204.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
# OUTPUT 2: SAVE THE DATAFRAME TO A FILE (INCUDING ALL CLASSES OF FILTERS)
write.csv(twins_dt, 'Data/twins_dt_20241204_1069.csv', quote = FALSE, row.names = F)

######################################################################################################
# OUTPUT 3: SAVE THE DATAFRAME TO A FILE (MUTATION + FILTERING STATUS ONLY)
twins_dt_info = twins_dt[, c('mut_ID', columns_req_filters, 'sum_req_filters'), with=FALSE]
write.csv(twins_dt_info, 'Data/twins_status_filters_20241204_1069.csv', quote = FALSE, row.names = F)

######################################################################################################
# OUTPUT 4: PLOT DISTIRBUTON OF MUTATION CLASSES AFTER FILTERING 

# Plot distribution of mutational classes on full and filtered dataset
muts_dt = twins_dt[, c('mut_ID', 'Ref', 'Alt', 'sum_req_filters'), with=FALSE]
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
  mut_type == 'G>A', 'C>T',
  mut_type == 'T>A', 'T>A',
  mut_type == 'T>G', 'T>G',
  mut_type == 'C>G', 'C>G',
  mut_type == 'C>A', 'C>A',
  mut_type == 'T>C', 'T>C',
  mut_type == 'C>T', 'C>T'))]

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

mut_counts_passed[, freq_norm := N/ length(muts_included)]
mut_counts_failed[, freq_norm := N/ (386412 - length(muts_included))]

mut_counts_passed[, QC := 'passed']
mut_counts_failed[, QC := 'failed']
mut_counts_qc = data.table(rbind(mut_counts_passed, mut_counts_failed))
mut_counts_qc[, V1 := factor(V1, levels = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))]
# keep the order that is applied in plots of mutation signatures

num_mut = sum(mut_counts_qc[,N])
ggplot(data=mut_counts_qc, aes(x=V1, y=freq_norm, fill=QC)) +
  geom_bar(stat='identity', position = 'dodge')+
  scale_fill_manual(values = c('#e71919', '#5fd2ef'))+
  theme_classic(base_size = 10)+
  labs(x = 'Mutation type', y = 'Fraction of mutations in the category', fill = 'Quality control', 
       title = glue('All mutations ({num_mut})'))
ggsave('Results/20241114_p1_mut_types_QC_all.pdf', width = 5, height = 3.5)

# compare removing germline mutations (many likely will be real)
muts_dt2 = twins_dt[f7_likelyGermline_bothTwins==0, c('mut_ID', 'Ref', 'Alt', 'sum_req_filters'), with=FALSE]
muts_dt2[, status := as.factor(fcase(
  sum_req_filters == 0, 'passed',
  sum_req_filters > 0, 'failed'
))]
muts_dt2[, mut_type := paste(Ref, Alt, sep='>')]
muts_dt2[, mut_type := as.factor(fcase(
  mut_type == 'A>T', 'T>A',
  mut_type == 'A>C', 'T>G',
  mut_type == 'G>C', 'C>G',
  mut_type == 'G>T', 'C>A',
  mut_type == 'A>G', 'T>C',
  mut_type == 'G>A', 'C>T',
  mut_type == 'T>A', 'T>A',
  mut_type == 'T>G', 'T>G',
  mut_type == 'C>G', 'C>G',
  mut_type == 'C>A', 'C>A',
  mut_type == 'T>C', 'T>C',
  mut_type == 'C>T', 'C>T'))]

mut_counts2 = data.table(sort(table(muts_dt2[,mut_type])), decreasing=TRUE)
mut_counts2[,V1:=as.factor(V1)]

# plot distribution of mutations of failed and passed mutations
mut_counts_passed2 = data.table(sort(table(muts_dt2[status == 'passed',mut_type])))
mut_counts_failed2 = data.table(sort(table(muts_dt2[status == 'failed',mut_type])))

mut_counts_passed2[, freq_norm := N/ length(muts_included)]
mut_counts_failed2[, freq_norm := N/ (dim(muts_dt2)[1] - length(muts_included))]

mut_counts_passed2[, QC := 'passed']
mut_counts_failed2[, QC := 'failed']
mut_counts_qc2 = data.table(rbind(mut_counts_passed2, mut_counts_failed2))
mut_counts_qc2[, V1 := as.factor(V1)]
mut_counts_qc2[, V1 := factor(V1, levels = c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'))]
# keep the order that is applied in plots of mutation signatures

num_mut2 = sum(mut_counts_qc2[,N])
ggplot(data=mut_counts_qc2, aes(x=V1, y=freq_norm, fill=QC)) +
  geom_bar(stat='identity', position = 'dodge')+
  scale_fill_manual(values = c('#e71919', '#5fd2ef'))+
  theme_classic(base_size = 10)+
  labs(x = 'Mutation type', y = 'Fraction of mutations in the category', fill = 'Quality control', 
       title = glue('Excluded putative germline ({num_mut2})'))
ggsave('Results/20241114_p1_mut_types_QC_germline_excluded.pdf', width = 5, height = 3.5)

######################################################################################################
# OUTPUT 4: PLOT DISTIRBUTON OF MUTATIONS IN TRI-NUC CONTEXTS AFTER FILTERING 

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
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p1_mut_trins_allmuts.pdf', width = 8, height = 2.5)

# plot tri-nucleotide context for the filtered set of mutations  
dtf = twins_dt[sum_req_filters==0]
mnr = dim(dtf)[1]
mybed = twins_dt[sum_req_filters==0,c('Chrom', 'Pos', 'Ref', 'Alt')]
trins = get_trinucs(mybed, BSgenome.Hsapiens.UCSC.hg38)
dtf$trins=trins

# plot the distribution of different mutations across different contexts 
mut_sign_counts2 = data.table(table(dtf[, trins]))
setnames(mut_sign_counts2, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts2[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts2[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts2[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

ggplot(data=mut_sign_counts2, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  labs(x = 'Context', y = 'Count', title = glue('Final set ({mnr})'))+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241114_p1_mut_trins_finalset.pdf', width = 8, height = 2.5)

# Create a series of plots that show what happens when you progressively apply different filters
# filters in the order of which removes the most mutations first 
filters = c('f7_likelyGermline_bothTwins', 'f6_lowQualRatio', 'f2_FailedIndelNearby30',
            'f4_mtrAndVaf', 'f3_lowDepthNormal', 'f3_highDepthNormal', 'f5_strandBiasMutOnly', 'f1_mappedY')

for (f in filters){
  
  f_name = tstrsplit(f, '_', fixed = T, keep = 2)
  
  twins_dt_f = twins_dt[get(f)!=1]
  mut_nr = dim(twins_dt_f)[1]
  mybed = twins_dt[get(f)!=1,c('Chrom', 'Pos', 'Ref', 'Alt')]
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
    labs(x = 'Context', y = 'Count', title = glue('Excluded: {f} ({mut_nr})'))+
    theme_classic(base_size = 15) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme( strip.background = element_blank())+ 
    theme(panel.spacing = unit(0, "lines"))+
    theme(strip.text.x = element_text(size = 13))+
    geom_hline(yintercept = 0, colour="black", size = 0.1)
  ggsave(glue('Results/202411019_p1_mut_trins_allmuts_{f}.pdf'), width = 8, height = 2.5)
}

# show how signatures change after consecutive filters are applied
for (i in seq_along(filters)){
  
  subset_filters = filters[1:i]
  cond = paste(subset_filters, '!= 1', collapse = ' & ')
  twins_dt_f = twins_dt[eval(parse(text = cond))]
  mut_nr = dim(twins_dt_f)[1]
  mybed = twins_dt[eval(parse(text = cond)), c('Chrom', 'Pos', 'Ref', 'Alt')]
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
    theme_classic(base_size = 15) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme( strip.background = element_blank())+ 
    theme(panel.spacing = unit(0, "lines"))+
    theme(strip.text.x = element_text(size = 13))+
    geom_hline(yintercept = 0, colour="black", size = 0.1)
  ggsave(glue('Results/20241114_p1_mut_trins_allmuts_filters_{i}.pdf'), width = 8, height = 2.5)
}

# show the type of mutations that are thrown away by each filter 
for (i in seq_along(filters)){
  
  if (i > 1){
    
    previous_filters = filters[1:i-1]
    current_filter=filters[i]
    cond = paste(paste(previous_filters, '!= 1', collapse = ' & '), '&', paste(current_filter, '==1'))
    print(cond)
    twins_dt_f = twins_dt[eval(parse(text = cond))]
    mut_nr = dim(twins_dt_f)[1]
    mybed = twins_dt[eval(parse(text = cond)), c('Chrom', 'Pos', 'Ref', 'Alt')]
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
      theme_classic(base_size = 15) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      theme( strip.background = element_blank())+ 
      theme(panel.spacing = unit(0, "lines"))+
      theme(strip.text.x = element_text(size = 13))+
      geom_hline(yintercept = 0, colour="black", size = 0.1)
    ggsave(glue('Results/20241114_p1_mut_trins_allmuts_filters_removed_{i}.pdf'), width = 8, height = 2.5) }
}

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# ALL DONE 


