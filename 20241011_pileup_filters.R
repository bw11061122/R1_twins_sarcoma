# Script to analyse the pileup (run 15/10/2024)
# 2024-10-11 - 2024-10-14
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# dataframe created in the script: 20241015_pileup_checks.R (pileup run 14/10/2024-15/10/2024)

# OUTPUT: list of mutations which pass required filters and will be used for phylogeny reconstruction

# NOTE ON ADDING COLUMNS WITH FILTERS 
# I want to add a column which says 0 if sample is fine and 1 if sample is wrong 

# COLUMNS (from the header)
##FORMAT-:DEP	Total reads covering this position (for subs del positions should be ignored)
##FORMAT-:FAZ	Reads presenting A for this position, forward strand
##FORMAT-:FCZ	Reads presenting C for this position, forward strand
##FORMAT-:FGZ	Reads presenting G for this position, forward strand
##FORMAT-:FTZ	Reads presenting T for this position, forward strand
##FORMAT-:MDR	Variant allele read directions 0=no reads; 1=Forward; 2=Reverse; 3=Forward + Reverse
##FORMAT-:MTR	Reads reporting the variant allele
##FORMAT-:OFS	Original filter status as defined in input vcf FILTER field
##FORMAT-:RAZ	Reads presenting A for this position, reverse strand
##FORMAT-:RCZ	Reads presenting C for this position, reverse strand
##FORMAT-:RGZ	Reads presenting G for this position, reverse strand
##FORMAT-:RTZ	Reads presenting T for this position, reverse strand
##FORMAT-:VAF	Variant Allele Fraction (excludes ambiguous reads if any)
##FORMAT-:WDR	Reference allele read directions 0=no reads; 1=Forward; 2=Reverse; 3=Forward + Reverse
##FORMAT-:WTR	Reads reporting the reference allele
##FILTER-:1	New filter status 1=not called in any sample
##FILTER-:2	New filter status 2=called in at least one sample
##FILTER-:3	New filter status 3=called + passed
##FILTER-:BD	Location from bed file
##ORIGINAL_FILTER-:CR	Position falls within a centromeric repeat using the supplied bed file
##ORIGINAL_FILTER-:DTH	Less than 1/3 mutant alleles were >= 25 base quality
##ORIGINAL_FILTER-:GI	Position falls within a germline indel using the supplied bed file
##ORIGINAL_FILTER-:HSD	Position falls within a high sequencing depth region using the supplied bed file
##ORIGINAL_FILTER-:MN	More than 0.05 of mutant alleles that were >= 15 base quality found in the matched normal
##ORIGINAL_FILTER-:MNP Tumour sample mutant allele proportion - normal sample mutant allele proportion < 0.2
##ORIGINAL_FILTER-:MQ	Mean mapping quality of the mutant allele reads was < 21
##ORIGINAL_FILTER-:PH	Mutant reads were on one strand (permitted proportion on other strand: 0.04), and mean mutant base quality was less than 21
##ORIGINAL_FILTER-:PT	Mutant alleles all on one direction of read (1rd allowed on opposite strand) and in second half of the read. Second half of read contains the motif GGC[AT]G in sequenced orientation and the mean base quality of all bases after the motif was less than 20
##ORIGINAL_FILTER-:RP	Coverage was less than 8 and no mutant alleles were found in the first 2/3 of a read (shifted 0.08 from the start and extended 0.08 more than 2/3 of the read length)
##ORIGINAL_FILTER-:SE	Coverage is >= 10 on each strand but mutant allele is only present on one strand
##ORIGINAL_FILTER-:SR	Position falls within a simple repeat using the supplied bed file
##ORIGINAL_FILTER-:VUM	Position has >= 3 mutant allele present in at least 1 percent unmatched normal samples in the unmatched VCF

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
library(VGAM)
library(grid)

###################################################################################################################################
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

###################################################################################################################################
# Specify settings for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#d37e03"
col_PD63383 = "#d4b17e"
col_tumour_PD62341 = "#980505"
col_tumour_PD63383 = "#eb6767"
col_normal_PD62341 = "#0785a5"
col_normal_PD63383 = "#70c6db"
col_bar = '#8c08d3'

######################################################################################################
# SAMPLES
# Examine which samples are being analysed 
samples = colnames(twins_dt[,c(seq(2, 330, 15))])
samples = sapply(strsplit(samples,"_"), `[`, 1)
paste('Number of samples analysed:', length(samples)) # remove the PD38is_wgs sample: not from the samples of interest

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
# CHECKS 

# Determine the number of unique mutations identified in at least one sample (ie how many mutations do I have?)
paste('Number of unique mutations identified across samples:', dim(twins_dt)[1]) # 362,814
paste('Number of unique mutations identified across samples (check):', length(twins_dt[, mut_ID] %>% unlist() %>% unique())) # 362,814

# Check that class of each column makes sense
sapply(twins_dt, class) # all fine 

######################################################################################################
# 1 FILTERING OUT READS MAPPING TO CHROMOSOME Y: EXCLUDE MUTATIONS MAPPED TO Y

twins_dt[, f1_mappedY := as.numeric(Chrom=='chrY')]
paste('Number of mutations mapped to chromosome Y:', dim(twins_dt[f1_mappedY==1])[1]) # 277

######################################################################################################
# 2 FILTERING BASED ON DISTANCE TO INDELS (DATA FROM PINDEL)
# Remove substitutions present next to any indels, as these are often poor quality (indels can lead to mapping issues) 

# Read indels (coordinates of low and high quality indels)
indel_all = data.table(read.csv('Data/PD62341_PD63382_unmatched_indels_pass_or_fail.txt', sep = '\t', header = FALSE))
setnames(indel_all, c('V1', 'V2', 'V3', 'V4'), c('Chrom', 'pos', 'ref', 'alt'))
indel_all[,mut_ID := paste(Chrom, pos, ref, alt, sep = '_')]
indel_all_id = indel_all[, mut_ID] %>% unlist()
paste("Number of indels identified in the dataset:", length(indel_all_id)) # 2,728,148

# exclude reads where reference not mapped
indel_all_dt = indel_all[!is.na(ref)]
indel_all_dt[, length := str_count(alt) - str_count(ref)]

# add start and end coordinates depending on if you have insertion or deletion
# I discussed with Henry and he recommends looking at +/- 10-30 bp (I think it may be better to be less stringent)
indel_all_dt[, start30 := pos - 30]
indel_all_dt[, end30 := pos + abs(length) + 30]

indel_all_dt = indel_all_dt[, c('Chrom', 'start30', 'end30'), with=FALSE]
indel_all_dt[, chrom_start := paste(Chrom, start30, sep = '_')]
indel_all_dt[, chrom_end := paste(Chrom, end30, sep = '_')]
chrom_start = indel_all_dt[, chrom_start] %>% unique() %>% unlist()
chrom_end = indel_all_dt[, chrom_end] %>% unique() %>% unlist()

twins_coord = twins_dt[, c('mut_ID', 'Chrom', 'Pos'), with=FALSE]
twins_coord[, coord := paste(Chrom, Pos, sep = '_')]

mut_indels = c()
for (chr in twins_coord[, Chrom] %>% unlist() %>% unique()){
  tc = twins_coord[Chrom==chr]
  ind = indel_all_dt[Chrom==chr]
  starts = ind[, start30] %>% unlist() %>% unique()
  ends = ind[, end30] %>% unlist() %>% unique()
  res = tc[ind, on = .(Pos >= start30, Pos <= end30), nomatch = 0]
  mut_indels = c(mut_indels, res[, mut_ID] %>% unlist())
}

twins_dt[, f2_FailedIndelNearby30 := as.numeric(mut_ID %in% mut_indels)]
paste('Number of mutations near poor-quality indels:', dim(twins_dt[f2_FailedIndelNearby30==1])[1]) # 13,524

######################################################################################################
# 3 COVERAGE, NUMBER OF VARIANT READS, VAF

# check that DEPTH = MTR + WTR (mutant + wt reads)
paste('Fraction of DEP = MTR + WTR:', dim(twins_dt[PD62341ak_MTR+PD62341ak_WTR==PD62341ak_DEP])[1] / dim(twins_dt)[1]) # 1

# for each mutation, get data on MTR (nr mutant reads), DEP (total nr of reads), VAF
# create subsetted dataframes with this information only
twins_mtr = twins_dt[,c('mut_ID', samples_mtr), with=FALSE]
twins_dep = twins_dt[,c('mut_ID', samples_dep), with=FALSE]
twins_vaf = twins_dt[,c('mut_ID', samples_vaf), with=FALSE]

# determine min, max and mean nr of variant reads / depth / vaf
twins_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_mtr]
twins_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_mtr]
twins_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_mtr]
twins_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = samples_mtr]

twins_dep[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_dep]
twins_dep[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_dep]
twins_dep[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_dep]
twins_dep[,max_dep := apply(.SD, 1, max), .SDcols = samples_dep]

twins_vaf[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_vaf]
twins_vaf[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_vaf]
twins_vaf[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_vaf]
twins_vaf[,max_vaf := apply(.SD, 1, max), .SDcols = samples_vaf]

# determine this specifically in normal samples (in case of ploidy changes in the tumour)
twins_mtr[,mean_mtr_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr[,median_mtr_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr[,min_mtr_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_mtr]
twins_mtr[,max_mtr_normal := apply(.SD, 1, max), .SDcols = samples_normal_mtr]

twins_dep[,mean_dep_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_dep]
twins_dep[,median_dep_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_dep]
twins_dep[,min_dep_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_dep]
twins_dep[,max_dep_normal := apply(.SD, 1, max), .SDcols = samples_normal_dep]

twins_vaf[,mean_vaf_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf[,median_vaf_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf[,min_vaf_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = samples_normal_vaf]
twins_vaf[,max_vaf_normal := apply(.SD, 1, max), .SDcols = samples_normal_vaf]

######################################################################################################
# 3 FILTERING BASED ON COVERAGE (total number of reads covering a base)

# Filter out mutations which have consistently low or high coverage
# Based on the histogram (20241011_pileup_qcplots.R), the appropriate boundaries are perhaps 20 and 50
# I'm debating what exact coverage value would be best to use to max specificity / sensitivity 
# I decided to go less stringent than this to not exclude good mutations 
# Default cut-off in Sequoia (clonal samples) is median > 10x and < 500x, but it is perhaps better to use thresholds based on the coverage distirbution in my samples 
paste('Number of sites with median coverage below 15 (normal samples):', dim(twins_dep[median_dep_normal < 15])[1]) # 9,215
paste('Number of sites with median coverage above 60 (normal samples):', dim(twins_dep[median_dep_normal > 60])[1]) # 500

# Filter out mutations where there are no reads in at least 1 sample (cannot tell status in the sample)
twins_dep[,sum_no_reads := rowSums(.SD==0), .SDcols = cols_normal_dep]
paste('Number of mutations w/o any reads in at least one sample:', dim(twins_dep[sum_no_reads>0])[1]) # 2,720
# I am only looking at normal samples here because there could be weird copy number changes in tumour samples

# Add columns with specific filters
twins_dep[,f3_lowDepthNormal := as.numeric(median_dep_normal < 15)] # 1 if does not pass the filter 
twins_dep[,f3_highDepthNormal := as.numeric(median_dep_normal > 60)] # 1 if does not pass the filter
twins_dep[,f3_noReadsMapped := as.numeric(sum_no_reads > 0)] # 1 if does not pass the filter

######################################################################################################
# FILTERING BASED ON NUMBER OF VARIANT READS (MTR)

# 1 count samples with MTR > 4 and with MTR = 0
twins_mtr[,sum_mtr0 := rowSums(.SD==0), .SDcols = cols_mtr]
twins_mtr[,sum_mtr4 := rowSums(.SD>=4), .SDcols = cols_mtr]
twins_mtr[,sum_mtr0_normal := rowSums(.SD==0), .SDcols = cols_normal_mtr] # specific to normal samples (tumour lost chr1/18)
twins_mtr[,sum_mtr4_normal := rowSums(.SD>=4), .SDcols = cols_normal_mtr] # specific to normal samples (tumour lost chr1/18)

# 2 identify mutations present in all samples (by MTR = 4 cut-off)
paste('Number of mutations identified in all samples (MTR >= 4 in all):', dim(twins_mtr[sum_mtr4==length(samples)])[1]) # 322,563

# 3 identify mutations present in no samples (by MTR = 4 cut-off)
paste('Number of mutations identified in no samples (MTR < 4 in all):', dim(twins_mtr[sum_mtr4==0])[1]) # 5,310

# 4 identify mutations present in ONLY one sample (by MTR = 4 cut-off)
# Those probably won't be any useful but we may want to check they really are absent everywhere else
# A lot of those look a bit dodgy / in regions with truncations / very few reads mapped or poor mapping quality  
paste('Number of mutations identified in only one sample (at MTR >= 4):', dim(twins_mtr[sum_mtr4==1])[1]) # 1,704
paste('Number of mutations identified in only one NORMAL sample (at MTR >= 4):', dim(twins_mtr[sum_mtr4_normal==1])[1]) # 1,760

# 5 check for mutations absent in one sample
paste('Number of mutations absent in only one sample (at MTR >= 4):', dim(twins_mtr[sum_mtr0==1])[1]) # 10,389 - most on chr1/18 due to loss in tumour samples
paste('Number of mutations absent in only one NORMAL sample (at MTR >= 4):', dim(twins_mtr[sum_mtr0_normal==1])[1]) # 1,867

# Add filter columns
twins_mtr[,f4_mtr4_presentInAll := as.numeric(sum_mtr4==length(samples))] # present in all samples with MTR >= 4 # 322,563
twins_mtr[,f4_mtr4_presentInNone := as.numeric(sum_mtr4==0)] # no sample with MTR >= 4 # 5310
twins_mtr[,f4_mtr4_presentInOne := as.numeric(sum_mtr4==1)] # only 1 sample with MTR >= 4 # 1704
twins_mtr[,f4_mtr4_absentInOne := as.numeric(sum_mtr4==length(samples_names)-1)] # only 1 sample with MTR < 4 # 3504

# Get a list of mutations which pass the MTR filter (they are mapped with 4 reads in at least one case)
muts_mapped_mtr = twins_mtr[f4_mtr4_presentInNone==0, mut_ID] %>% unlist()

######################################################################################################
# FILTERING BASED ON VAF

# require VAF > 0.1 in at least one sample 
twins_vaf[, f5_absentAllVaf := as.numeric(max_vaf < 0.1)] 
paste("Number of mutations which don't reach VAF 0.1 in any sample:", dim(twins_vaf[f5_absentAllVaf==1])[1]) # 1,479

# exclude mutations with VAF > 0.4 in all samples (i.e., present everywhere at high levels)
twins_vaf[, f5_presentAllVaf := as.numeric(min_vaf > 0.4 | median_vaf_normal > 0.6)] 
paste("Number of mutations present in all samples (VAF > 0.4):", dim(twins_vaf[f5_presentAllVaf==1])[1]) # 23,037

######################################################################################################
# FILTERING BASED ON GOOD MAPPING IN AT LEAST 1 SAMPLE

# Require correct mapping (MTR >= 4 and VAF >= 0.1) in THE SAME sample - i.e., require a sample where everything is correct
twins_mtr_vaf = twins_dt[, c('mut_ID', samples_mtr, samples_vaf), with=FALSE]

# for each sample, check if the sample has MTR >= 4 and VAF >= 0.1
# 0 if both are fulfilled, 1 if either is not 

for(sample in samples_names){
    twins_sample = twins_mtr_vaf[, c('mut_ID', glue('{sample}_MTR'), glue('{sample}_VAF'))]
    setnames(twins_sample, c(glue('{sample}_MTR'), glue('{sample}_VAF')), c('mtr', 'vaf'))
    twins_sample[, glue('f6_{sample}') :=  1 - as.numeric(mtr >= 4) * as.numeric(vaf >= 0.1)]
    twins_mtr_vaf = merge(twins_mtr_vaf, twins_sample[, c('mut_ID', glue('f6_{sample}')), with = FALSE], by = 'mut_ID')
}

f6_samples = paste0('f6_', samples_names)
twins_mtr_vaf[, f6_sum := rowSums(.SD == 0), .SDcols = f6_samples]
twins_mtr_vaf[, f6_mtrAndVaf := as.numeric(f6_sum == 0)] # 0 if VAF > 0.1 and MTR > 4 in at least one sample
paste("Number of mutations which are NOT well mapped by both MTR and VAF in at least 1 sample:", dim(twins_mtr_vaf[f6_mtrAndVaf==1])[1]) # 5,501

######################################################################################################
# EXCLUDE LIKELY GERMLINE MUTATIONS BASED ON EXACT BINOMIAL

twins_normal = twins_dt[, c('mut_ID', samples_normal_mtr, samples_normal_dep), with=FALSE]
twins_normal[, sum_MTR_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_normal[, sum_DEP_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_normal[, sum_MTR_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_normal[, sum_DEP_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_normal[, sum_MTR_all := sum_MTR_PD62341 + sum_MTR_PD63383]
twins_normal[, sum_DEP_all := sum_DEP_PD62341 + sum_DEP_PD63383]

# remove rows where there is no coverage (total or in either twin)
twins_normal = twins_normal[(sum_DEP_all!=0 & sum_DEP_PD62341!=0 & sum_DEP_PD63383!=0)]

# function to run binomial test
bin_test <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=
                                                  c("two.sided"), conf.level = 0.95)$p.value}

# stat test 
twins_normal[, p_val_PD62341 := mapply(bin_test, sum_MTR_PD62341, sum_DEP_PD62341)]
twins_normal[, p_val_PD63383 := mapply(bin_test, sum_MTR_PD63383, sum_DEP_PD63383)]
twins_normal[, p_val_all := mapply(bin_test, sum_MTR_all, sum_DEP_all)]

# Add column as a filter (NOTE: H0 is germline so if p > 0.01 you don't reject H0)
twins_normal[, f7_likelyGermline_PD62341 := as.numeric(p_val_PD62341 >= 0.01/dim(twins_normal)[1])]
twins_normal[, f7_likelyGermline_PD63383 := as.numeric(p_val_PD63383 >= 0.01/dim(twins_normal)[1])] 
twins_normal[, f7_likelyGermline_bothTwins := as.numeric(p_val_all >= 0.01/dim(twins_normal)[1])]  
paste('Number of likely germline mutations (PD62341):', dim(twins_normal[f7_likelyGermline_PD62341==1])[1]) # 352,685
paste('Number of likely germline mutations (PD63383):', dim(twins_normal[f7_likelyGermline_PD63383==1])[1]) # 352,881
paste('Number of likely germline mutations (both twins):', dim(twins_normal[f7_likelyGermline_bothTwins==1])[1]) # 350,300

# write IDs of likely germline mutations to csv
germline_mutations = twins_normal[f7_likelyGermline_bothTwins==1, mut_ID] %>% unlist()
write.table(germline_mutations, 'Data/mutations_likely_germline_20241026.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
# FILTERING BASED ON STRAND BIAS (IN VARIANT READS)

# Subset strand information (only for normal samples: ploidy changes can affect biases)
forward_cols = grep("_FAZ|_FCZ|_FTZ|_FGZ", names(twins_dt), value = TRUE) 
forward_cols = grep('PD6', forward_cols, value=TRUE) # I don't want PD38is_wgs
forward_cols = grep("PD62341v|PD62341q|PD62341aa|PD62341ad|PD63383w|PD63383t|PD63383u|PD63383ae|PD63383ak|PD63383bb|PD62341h|PD62341n", forward_cols, value=TRUE) # only select normal samples
reverse_cols = grep("_RAZ|_RCZ|_RTZ|_RGZ", names(twins_dt), value = TRUE) 
reverse_cols = grep("PD62341v|PD62341q|PD62341aa|PD62341ad|PD63383w|PD63383t|PD63383u|PD63383ae|PD63383ak|PD63383bb|PD62341h|PD62341n", reverse_cols, value=TRUE) # only select normal samples
twins_strands = twins_dt[,c('mut_ID', 'Ref', 'Alt', forward_cols, reverse_cols), with=FALSE]
twins_strands[, (forward_cols) := lapply(.SD, as.numeric), .SDcols = forward_cols]
twins_strands[, (reverse_cols) := lapply(.SD, as.numeric), .SDcols = reverse_cols]

# for each mutation, select relevant forward and reverse reads (corresponding to the Alt variant)
twins_strands = twins_dt[,c('mut_ID', 'Ref', 'Alt', forward_cols, reverse_cols), with=FALSE]
twins_strands[, (forward_cols) := lapply(.SD, as.numeric), .SDcols = forward_cols]
twins_strands[, (reverse_cols) := lapply(.SD, as.numeric), .SDcols = reverse_cols]

twins_strands[, forward_mut := sapply(1:.N, function(row) {
  alt = Alt[row]
  cols_to_sum = names(twins_strands)[grepl(paste0('F', alt, 'Z'), names(twins_strands))]
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)
})]

twins_strands[, reverse_mut := sapply(1:.N, function(row) {
  alt = Alt[row]
  cols_to_sum = names(twins_strands)[grepl(paste0('R', alt, 'Z'), names(twins_strands))]
  sum(unlist(.SD[row, c(cols_to_sum), with=FALSE]), na.rm=TRUE)
})]

# aggregate across all samples 
twins_strands[, sum_all_mut := forward_mut + reverse_mut]

# Filter out mutations where sum of forward + reverse = 0
twins_strands = twins_strands[sum_all_mut!=0] 
twins_strands[, p_val_mut := mapply(bin_test, forward_mut, sum_all_mut)]

# correct for multiple testing (Bonferroni correction for nr of tests)
twins_strands[, significance_mut := as.factor(fcase( 
  p_val_mut < 0.01/dim(twins_strands)[1], 'significant', # differs from 0.5 so biased and maybe we don't want it
  p_val_mut >= 0.01/dim(twins_strands)[1], 'not_significant'
))] # really want to minimize false positives because there is a risk of throwing away high-quality mutations

# Add column as a filter 
twins_strands[, f8_strandBiasMut := as.numeric(significance_mut == 'significant')] 
paste('Number of mutations which show strand bias (across mutant reads):', dim(twins_strands[f8_strandBiasMut==1])[1]) # 4,764

# strand bias in normal reads 
twins_strands[,sum_FAZ := rowSums(.SD), .SDcols = patterns("_FAZ", cols = names(twins_strands))]
twins_strands[,sum_FTZ := rowSums(.SD), .SDcols = patterns("_FTZ", cols = names(twins_strands))]
twins_strands[,sum_FCZ := rowSums(.SD), .SDcols = patterns("_FCZ", cols = names(twins_strands))]
twins_strands[,sum_FGZ := rowSums(.SD), .SDcols = patterns("_FGZ", cols = names(twins_strands))]

twins_strands[,sum_RAZ := rowSums(.SD), .SDcols = patterns("_RAZ", cols = names(twins_strands))]
twins_strands[,sum_RTZ := rowSums(.SD), .SDcols = patterns("_RTZ", cols = names(twins_strands))]
twins_strands[,sum_RCZ := rowSums(.SD), .SDcols = patterns("_RCZ", cols = names(twins_strands))]
twins_strands[,sum_RGZ := rowSums(.SD), .SDcols = patterns("_RGZ", cols = names(twins_strands))]

twins_strands[,sum_forward := sum_FAZ + sum_FCZ + sum_FTZ + sum_FGZ]
twins_strands[,sum_reverse := sum_RAZ + sum_RCZ + sum_RTZ + sum_RGZ]
twins_strands[,sum_all := sum_forward + sum_reverse]

# Filter out mutations where sum of forward + reverse = 0
twins_strands = twins_strands[sum_all!=0] 
bin_test <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative = c("two.sided"), conf.level = 0.95)$p.value}
twins_strands[, p_val := mapply(bin_test, sum_forward, sum_all)]

# correct for multiple testing (Bonferroni correction for nr of tests)
twins_strands[, significance := as.factor(fcase( 
  p_val < 0.01/dim(twins_strands)[1], 'significant', # differs from 0.5 so biased and maybe we don't want it
  p_val >= 0.01/dim(twins_strands)[1], 'not_significant'
))]

# Add column as a filter 
twins_strands[, f8_strandBias := as.numeric(significance == 'significant')] 
paste('Number of mutations which show strand bias (across all reads):', dim(twins_strands[f8_strandBias==1])[1]) # 10,005

######################################################################################################
# FILTERING BASED ON RATIO OF READS MAPPED IN HQ AND LQ PILEUP 

# read in files with high and low quality pileup 
twins_hq = data.table(read.csv('Data/pileup_merged_20241016.tsv'))
twins_lq = data.table(read.csv('Data/pileup_merged_20241025.tsv'))

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

# identify mutations where reads are only present in the low-quality pileup 
numeric_columns = intersect(names(twins_mtr_hq)[sapply(twins_mtr_hq, is.numeric)], 
                            names(twins_mtr_lq)[sapply(twins_mtr_lq, is.numeric)])

lq_minus_hq = twins_mtr_lq[, ..numeric_columns] - twins_mtr_hq[, ..numeric_columns]
lq_minus_hq = cbind(twins_mtr_lq[, mut_ID], lq_minus_hq, twins_mtr_lq[, sum_mtr_mapped], twins_mtr_hq[, sum_mtr_mapped])
setnames(lq_minus_hq, c('V1', 'V3', 'V4'), c('mut_ID', 'sum_mtr_lq', 'sum_mtr_hq'))
setnames(lq_minus_hq, c(samples_mtr), c(paste0(samples_mtr, '_diff')))

# add count in high quality pileup
lq_minus_hq = cbind(lq_minus_hq, twins_mtr_hq[, 2:23, with=FALSE])

# identify samples where the MTR (HQ) is <= 2 and MTR (LQ) >= 4
for (sample in samples_names){
  lq_minus_hq_sample = lq_minus_hq[, c('mut_ID', glue('{sample}_MTR'), glue('{sample}_MTR_diff'))]
  setnames(lq_minus_hq_sample, c(glue('{sample}_MTR'), glue('{sample}_MTR_diff')), c('mtr', 'diff'))
  lq_minus_hq_sample[, glue('f9_{sample}') :=  as.numeric(mtr <= 2) * as.numeric(diff >= 4)]
  lq_minus_hq = merge(lq_minus_hq, lq_minus_hq_sample[, c('mut_ID', glue('f9_{sample}')), with = FALSE], by = 'mut_ID')
}

f9_samples = paste0('f9_', samples_names)
lq_minus_hq[, f9_sum := rowSums(.SD >= 1), .SDcols = f9_samples]
lq_minus_hq[, f9_lowQualPresent := as.numeric(f9_sum >= 1)] # assign 1 if in >= 1 sample presence / absence differs b/n LQ and HQ
muts_exclude_lqpresent = lq_minus_hq[f9_lowQualPresent==1, mut_ID] %>% unlist()
paste("Number of mutations with variable 1/0 by HQ/LQ in >= 1 samples:", length(muts_exclude_lqpresent)) # 8,851

######################################################################################################
# CREATE A DATAFRAME WHICH HAS ALL FILTERS TOGETHER

# Extract dataframes with filters and merge with the main dataframe
dep_filters = twins_dep[, c('mut_ID', "f3_highDepthNormal", 'f3_lowDepthNormal', 'f3_noReadsMapped'), with = FALSE]
mtr_filters = twins_mtr[, c('mut_ID', 'f4_mtr4_presentInAll', 'f4_mtr4_presentInNone', 'f4_mtr4_presentInOne'), with = FALSE]
vaf_filters = twins_vaf[, c('mut_ID', "f5_presentAllVaf", 'f5_absentAllVaf') , with = FALSE]
mtr_vaf_filters = twins_mtr_vaf[, c('mut_ID', 'f6_mtrAndVaf'), with=FALSE]
germline_filters = twins_normal[, c('mut_ID', "f7_likelyGermline_PD62341", "f7_likelyGermline_PD63383", "f7_likelyGermline_bothTwins"), with = FALSE]
strand_filters = twins_strands[, c('mut_ID', 'f8_strandBias', 'f8_strandBiasMut'), with = FALSE]

twins_dt_merge = merge(twins_dt, dep_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, mtr_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, vaf_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, mtr_vaf_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, germline_filters, by = 'mut_ID')
twins_dt_filters = merge(twins_dt_merge, strand_filters, by = 'mut_ID')

# add filter from beta-binomial (rho)
twins_dt_filters[, f9_lowQualRatio := as.numeric(mut_ID %in% muts_exclude_ratio)]
twins_dt_filters[, f9_lowQualPresent := as.numeric(mut_ID %in% muts_exclude_lqpresent)]

######################################################################################################
# CHECK NR OF MUTATIONS EXCLUDED IN EACH FILTER

paste('Number of all mutations examined:', dim(twins_dt_filters)[1]) # 361,699
# Note that some mutations have been excluded due to depth = 0 for all or some samples 

dim(twins_dt_filters[f1_mappedY==1])[1] # 147
dim(twins_dt_filters[f2_FailedIndelNearby30==1])[1] # 13,426
dim(twins_dt_filters[f3_lowDepthNormal==1])[1] # 8,635
dim(twins_dt_filters[f3_highDepthNormal==1])[1] # 500
dim(twins_dt_filters[f3_noReadsMapped==1])[1] # 2,060
dim(twins_dt_filters[f4_mtr4_presentInAll==1])[1] # 322,563
dim(twins_dt_filters[f4_mtr4_presentInNone==1])[1] # 4,324
dim(twins_dt_filters[f4_mtr4_presentInOne==1])[1] # 1,583
dim(twins_dt_filters[f5_presentAllVaf==1])[1] # 22,241
dim(twins_dt_filters[f5_absentAllVaf==1])[1] # 912
dim(twins_dt_filters[f6_mtrAndVaf==1])[1] # 4,506
dim(twins_dt_filters[f7_likelyGermline_PD62341==1])[1] # 352,521
dim(twins_dt_filters[f7_likelyGermline_PD63383==1])[1] # 352,719
dim(twins_dt_filters[f7_likelyGermline_bothTwins==1])[1] # 350,192
dim(twins_dt_filters[f8_strandBias==1])[1] # 10,005
dim(twins_dt_filters[f8_strandBiasMut==1])[1] # 4,764
dim(twins_dt_filters[f9_lowQualRatio==1])[1] # 13,385
dim(twins_dt_filters[f9_lowQualPresent==1])[1] # 8,154

######################################################################################################
# Specify required filters 
columns_req_filters = c('f1_mappedY', 'f2_FailedIndelNearby30', 
                        'f3_lowDepthNormal', 'f3_highDepthNormal', 'f3_noReadsMapped',
                        'f4_mtr4_presentInAll', 'f4_mtr4_presentInNone', 'f4_mtr4_presentInOne', 
                        'f5_presentAllVaf', 'f5_absentAllVaf', 'f6_mtrAndVaf',
                        'f7_likelyGermline_bothTwins', 'f8_strandBiasMut', 'f9_lowQualRatio', 'f9_lowQualPresent')
twins_dt_filters[, sum_req_filters := rowSums(.SD), .SDcols = columns_req_filters] 
paste('Number of mutations that pass required filters:', dim(twins_dt_filters[sum_req_filters==0])[1]) # 563  

cols_info = c('VariantID', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Filter', 'Gene', 'Transcript', 'RNA', 'CDS',
              'Protein', 'Effect')
twins_filtered = twins_dt_filters[sum_req_filters==0]
twins_filtered_mtr = twins_filtered[, c('mut_ID', cols_info, cols_mtr), with=FALSE]
twins_filtered_mtr[, (cols_mtr) := lapply(.SD, as.numeric), .SDcols = cols_mtr]
twins_filtered_mtr[, sum_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_mtr[, sum_normal := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_mtr[, sum_PD62341 := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_mtr[, sum_PD63383 := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_mtr[, sum_tumour_PD62341 := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_mtr[, sum_tumour_PD63383 := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_mtr[, sum_normal_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_mtr[, sum_normal_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]
mut_included = twins_filtered[, mut_ID] %>% unlist()
paste('Number of retained mutations:', length(mut_included)) # 413

######################################################################################################
# SAVE FILTERED MUTATIONS TO A TXT FILE
paste('Number of mutations that passed required filters:', length(mut_included)) # 413
write.table(mut_included, 'Data/mutations_include_20241030_413.txt', quote = FALSE, col.names = F, row.names = F)

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
  theme_bw(base_size = 12)+
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
ggsave('Results/20241030_p1_mut_filters_types.pdf', width = 6, height = 4.5)

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
  theme_bw(base_size = 12)+
  labs(x = 'Mutation type', y = 'Fraction of all mutations in the category', title = 'Mutation types in all samples\nexcluded putative germline mutations')
ggsave('Results/20241030_p1_mut_filters_types_germline_excluded.pdf', width = 6, height = 4.5)

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
  theme_minimal(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", size = 0.1)
ggsave('Results/20241030_p1_mut_trins_allmuts.pdf', width = 10, height = 5.5)


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
  ggsave(glue('Results/20241030_p1_mut_trins_allmuts_{f}.pdf'), width = 10, height = 5.5)
}

# we want to show how signatures change after each filter is applied
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
  ggsave(glue('Results/20241030_p1_mut_trins_allmuts_filters_{i}.pdf'), width = 10, height = 5.5)
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
  ggsave(glue('Results/20241030_p1_mut_trins_allmuts_req_filters_{i}.pdf'), width = 10, height = 5.5)
}

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# ALL DONE 



