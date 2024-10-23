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
# 1 FILTERING BASED ON QC (FILTER = PASS)
# Note: checking on JBrowse, if all other filters are passed, there are no reasons to exclude a mutation

# Add a column with a filter which states if the mutation passes the CaVEMan filter 
twins_dt[, f1_PASS:=as.numeric(Filter!='PASS')] # 0 if PASS, 1 if NOT PASS
paste('Number of mutations that DO NOT pass the quality filter:', dim(twins_dt[Filter!='PASS'])[1]) # 11,363

######################################################################################################
# 2 FILTERING OUT READS MAPPING TO CHROMOSOME Y: EXCLUDE MUTATIONS MAPPED TO Y

twins_dt[, f2_mappedY := as.numeric(Chrom=='chrY')]
paste('Number of mutations mapped to chromosome Y:', dim(twins_dt[f2_mappedY==1])[1]) # 277

######################################################################################################
# COVERAGE, NUMBER OF VARIANT READS, VAF

# first, I want to check that DEPTH = MTR + WTR (mutant + wt reads)
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
# FILTERING BASED ON COVERAGE (total number of reads covering a base)

# Filter out mutations which have consistently low or high coverage
# Based on the histogram (QC plots below) I decided to use filters 15 and 60
# Only retain mutations which are covered by at least 15 reads in all normal samples
# Only retain mutations which are covered by max 60 reads in all normal samples
# Default cut-off in Sequoia (clonal samples) is median > 10x and < 500x, but it is perhaps better to use thresholds based on the coverage distirbution in my samples 
paste('Number of sites with min coverage below 15 (normal samples):', dim(twins_dep[min_dep_normal < 15])[1]) # 22,199
paste('Number of sites with max coverage above 60 (normal samples):', dim(twins_dep[max_dep_normal > 60])[1]) # 19,038

paste('Number of sites with min coverage below 15 (normal samples):', dim(twins_dep[median_dep_normal < 20])[1]) # 14,500
paste('Number of sites with max coverage above 60 (normal samples):', dim(twins_dep[median_dep_normal > 50])[1]) # 1,342

paste('Number of sites with median coverage below 20 (all samples):', dim(twins_dep[median_dep < 20])[1]) # 13,691
paste('Number of sites with median coverage above 50 (all samples):', dim(twins_dep[median_dep > 50])[1]) # 1,432

# Filter out mutations where there are no reads in at least 1 sample (cannot tell status in the sample)
twins_dep[,sum_no_reads := rowSums(.SD==0), .SDcols = cols_dep]
paste('Number of mutations w/o any reads in at least one sample:', dim(twins_dep[sum_no_reads>0])[1]) # 3144

# Add columns with specific filters
twins_dep[,f3_lowDepth := as.numeric(median_dep < 20)] # 1 if does not pass the filter 
twins_dep[,f3_highDepth := as.numeric(median_dep > 50)] # 1 if does not pass the filter
twins_dep[,f3_lowDepth_normal := as.numeric(median_dep_normal < 20)] # 1 if does not pass the filter 
twins_dep[,f3_highDepth_normal := as.numeric(median_dep_normal > 50)] # 1 if does not pass the filter
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
twins_mtr[,f4_mtr4_presentInAll := as.numeric(sum_mtr4==length(samples))] # present in all 
twins_mtr[,f4_mtr4_presentInNone := as.numeric(sum_mtr4==0)] # no sample with MTR >= 4
twins_mtr[,f4_mtr4_presentInOne := as.numeric(sum_mtr4==1)] # only 1 sample with MTR >= 4
twins_mtr[,f4_mtr4_absentInOne := as.numeric(sum_mtr4==length(samples)-1)] # present in all but one sample  

######################################################################################################
# FILTERING BASED ON VAF

# require VAF > 0.1 in at least one sample 
twins_vaf[, f5_maxVAF01 := as.numeric(max_vaf < 0.1)] 
paste("Number of mutations which don't reach VAF 0.1 in any sample:", dim(twins_vaf[f5_maxVAF01==1])[1]) # 1,479

# mutations with VAF = NA (these happen is MTR = 0)
twins_vaf[, f5_VAFna := as.numeric(is.na(mean_vaf))] 
paste("Number of mutations with median VAF == NA (sth weird going on):", dim(twins_vaf[f5_VAFna==1])[1]) # 446

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
twins_normal[, f6_likelyGermline_PD62341 := as.numeric(p_val_PD62341 >= 0.01/dim(twins_normal)[1])]
twins_normal[, f6_likelyGermline_PD63383 := as.numeric(p_val_PD63383 >= 0.01/dim(twins_normal)[1])] 
twins_normal[, f6_likelyGermline_bothTwins := as.numeric(p_val_all >= 0.01/dim(twins_normal)[1])]  
paste('Number of likely germline mutations (PD62341):', dim(twins_normal[f6_likelyGermline_PD62341==1])[1]) # 352,685
paste('Number of likely germline mutations (PD63383):', dim(twins_normal[f6_likelyGermline_PD63383==1])[1]) # 352,881
paste('Number of likely germline mutations (both twins):', dim(twins_normal[f6_likelyGermline_bothTwins==1])[1]) # 350,300

######################################################################################################
# FILTERING BASED ON DISTANCE TO INDELS (DATA FROM PINDEL)
# I will want to remove substitutions present next to any indels because these can lead to mapping issues 

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
indel_all_dt[, start10 := (fcase( 
  length <= 0, pos + length - 10,
  length > 0, pos - 10
))] 
indel_all_dt[, end10 := (fcase( 
  length <= 0, pos + 10,
  length > 0, pos + length + 10
))]

indel_all_dt = indel_all_dt[, c('Chrom', 'start10', 'end10'), with=FALSE]
indel_all_dt[, chrom_start := paste(Chrom, start10, sep = '_')]
indel_all_dt[, chrom_end := paste(Chrom, end10, sep = '_')]
chrom_start = indel_all_dt[, chrom_start] %>% unique() %>% unlist()
chrom_end = indel_all_dt[, chrom_end] %>% unique() %>% unlist()

twins_coord = twins_dt[, c('mut_ID', 'Chrom', 'Pos'), with=FALSE]
twins_coord[, coord := paste(Chrom, Pos, sep = '_')]

mut_indels = c()
for (chr in twins_coord[, Chrom] %>% unlist() %>% unique()){
  tc = twins_coord[Chrom==chr]
  ind = indel_all_dt[Chrom==chr]
  starts = ind[, start10] %>% unlist() %>% unique()
  ends = ind[, end10] %>% unlist() %>% unique()
  res = tc[ind, on = .(Pos >= start10, Pos <= end10), nomatch = 0]
  mut_indels = c(mut_indels, res[, mut_ID] %>% unlist())
}

twins_dt[, f7_FailedIndelNearby10 := as.numeric(mut_ID %in% mut_indels)]
paste('Number of mutations near poor-quality indels:', dim(twins_dt[f7_FailedIndelNearby10==1])[1]) # 8,528

######################################################################################################
# FILTERING BASED ON STRAND BIAS (IN ALL READS)
# Subset strand information (only for normal samples: ploidy changes can affect biases)

forward_cols = grep("_FAZ|_FCZ|_FTZ|_FGZ", names(twins_dt), value = TRUE) 
forward_cols = grep('PD6', forward_cols, value=TRUE) # I don't want PD38is_wgs
forward_cols = grep("PD62341v|PD62341q|PD62341aa|PD62341ad|PD63383w|PD63383t|PD63383u|PD63383ae|PD63383ak|PD63383bb|PD62341h|PD62341n", forward_cols, value=TRUE) # only select normal samples
reverse_cols = grep("_RAZ|_RCZ|_RTZ|_RGZ", names(twins_dt), value = TRUE) 
reverse_cols = grep("PD62341v|PD62341q|PD62341aa|PD62341ad|PD63383w|PD63383t|PD63383u|PD63383ae|PD63383ak|PD63383bb|PD62341h|PD62341n", reverse_cols, value=TRUE) # only select normal samples
twins_strands = twins_dt[,c('mut_ID', forward_cols, reverse_cols), with=FALSE]
twins_strands[, (forward_cols) := lapply(.SD, as.numeric), .SDcols = forward_cols]
twins_strands[, (reverse_cols) := lapply(.SD, as.numeric), .SDcols = reverse_cols]

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
twins_strands[,ratio_forward_reverse := sum_forward / sum_reverse]

twins_strands[,chr := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
twins_strands[,chr := as.factor(chr)]

# Filter out mutations where sum of forward + reverse = 0
twins_strands = twins_strands[sum_all!=0] # NB: checked on JBrowse - all of those are useless
bin_test <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=
                                                  c("two.sided"), conf.level = 0.95)$p.value}
twins_strands[, p_val := mapply(bin_test, sum_forward, sum_all)]

# correct for multiple testing (Bonferroni correction for nr of tests)
twins_strands[, significance := as.factor(fcase( 
  p_val < 0.05/dim(twins_strands)[1], 'significant', # differs from 0.5 so biased and maybe we don't want it
  p_val >= 0.05/dim(twins_strands)[1], 'not_significant'
))]

# Add column as a filter 
twins_strands[, f8_strandBias := as.numeric(significance == 'significant')] 
paste('Number of mutations which show strand bias:', dim(twins_strands[f8_strandBias==1])[1]) # 11,467

######################################################################################################
# FILTERING BASED ON STRAND BIAS (IN VARIANT READS, MDR)
# Only using normal sample data here: changes in tumour ploidy can affect biases
cols_mdr = grep("MDR", names(twins_dt), value = TRUE)
cols_mdr = grep("PD62341v|PD62341q|PD62341aa|PD62341ad|PD63383w|PD63383t|PD63383u|PD63383ae|PD63383ak|PD63383bb|PD62341h|PD62341n", cols_mdr, value=TRUE) # only select normal samples
twins_mdr = twins_dt[,c('mut_ID', cols_mdr), with=FALSE]

# find mutations which are 3 (code for forward + reverse okay)
twins_mdr[, sum_MDR := rowSums(.SD==3), .SDcols = cols_mdr]
twins_mdr[, f8_MDR3all := as.numeric(sum_MDR!=12)] 
paste('Number of mutations with some strand bias:', dim(twins_mdr[f8_MDR3all==1])[1]) # 24,237

# let's say that we want '3' in >= 10 samples (that leaves 2 samples with strand bias so still should be okay)
twins_mdr[, f8_MDR3min10 := as.numeric(sum_MDR<10)] 
paste("Number of mutations with some strand bias in variant reads (in >= 10/12 samples):", dim(twins_mdr[f8_MDR3min10==1])[1]) # 15,290

######################################################################################################
# CREATE A DATAFRAME WHICH HAS ALL FILTERS TOGETHER

# Extract dataframes with filters and merge with the main dataframe
dep_filters = twins_dep[, c('mut_ID', "f3_lowDepth", "f3_highDepth", 'f3_lowDepth_normal', 'f3_highDepth_normal', "f3_noReadsMapped"), with = FALSE]
mtr_filters = twins_mtr[, c('mut_ID', 'f4_mtr4_presentInAll', 'f4_mtr4_presentInNone', 'f4_mtr4_presentInOne', 'f4_mtr4_absentInOne'), with = FALSE]
vaf_filters = twins_vaf[, c('mut_ID', "f5_maxVAF01", "f5_VAFna") , with = FALSE]
germline_filters = twins_normal[, c('mut_ID', "f6_likelyGermline_PD62341", "f6_likelyGermline_PD63383", "f6_likelyGermline_bothTwins"), with = FALSE]
strand_filters = twins_strands[, c('mut_ID', 'f8_strandBias'), with = FALSE]
mdr_filters = twins_mdr[, c('mut_ID', 'f8_MDR3all', 'f8_MDR3min10'), with = FALSE]

twins_dt_merge = merge(twins_dt, dep_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, mtr_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, vaf_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, germline_filters, by = 'mut_ID')
twins_dt_merge = merge(twins_dt_merge, strand_filters, by = 'mut_ID')
twins_dt_filters = merge(twins_dt_merge, mdr_filters, by = 'mut_ID')

######################################################################################################
# CHECK NR OF MUTATIONS EXCLUDED IN EACH FILTER

paste('Number of all mutations examined:', dim(twins_dt_filters)[1]) # 362,287
# Note that some mutations have been excluded due to depth = 0 for all or some samples 

dim(twins_dt_filters[f1_PASS==1])[1] # 10,941
dim(twins_dt_filters[f2_mappedY==1])[1] # 149
dim(twins_dt_filters[f3_lowDepth==1])[1] # 13,322
dim(twins_dt_filters[f3_lowDepth_normal==1])[1] # 14,207
dim(twins_dt_filters[f3_highDepth==1])[1] # 1,432
dim(twins_dt_filters[f3_highDepth_normal==1])[1] # 1,342
dim(twins_dt_filters[f3_noReadsMapped==1])[1] # 2,617
dim(twins_dt_filters[f4_mtr4_presentInAll==1])[1] # 322,563
dim(twins_dt_filters[f4_mtr4_presentInNone==1])[1] # 4,789
dim(twins_dt_filters[f4_mtr4_presentInOne==1])[1] # 1,698
dim(twins_dt_filters[f5_maxVAF01==1])[1] # 1,274
dim(twins_dt_filters[f5_VAFna==1])[1] # 241
dim(twins_dt_filters[f6_likelyGermline_PD62341==1])[1] # 352,685
dim(twins_dt_filters[f6_likelyGermline_PD63383==1])[1] # 352,881
dim(twins_dt_filters[f6_likelyGermline_bothTwins==1])[1] # 350,300
dim(twins_dt_filters[f7_FailedIndelNearby10==1])[1] # 8,506 # not much of a difference
dim(twins_dt_filters[f8_strandBias==1])[1] # 11,467
dim(twins_dt_filters[f8_MDR3all==1])[1] # 23,710
dim(twins_dt_filters[f8_MDR3min10==1])[1] # 14,763

######################################################################################################
# CHECK NR OF MUTATIONS THAT PASS ALL FILTERS

columns_all_filters = c('f1_PASS', 'f2_mappedY', 'f3_lowDepth', 'f3_lowDepth_normal', 'f3_highDepth', 'f3_highDepth_normal', 'f3_noReadsMapped', 
                       'f4_mtr4_presentInAll', 'f4_mtr4_presentInNone', 'f4_mtr4_presentInOne', 'f4_mtr4_absentInOne',
                       'f5_maxVAF01', 'f5_VAFna', 'f6_likelyGermline_PD62341', 'f6_likelyGermline_PD63383', 'f6_likelyGermline_bothTwins', 'f7_FailedIndelNearby10', 'f8_strandBias', 'f8_MDR3all', 'f8_MDR3min10')

twins_dt_filters[, sum_all_filters := rowSums(.SD), .SDcols = columns_all_filters] 
paste('Number of mutations that pass all filters:', dim(twins_dt_filters[sum_all_filters==0])[1]) 
# Not all filters should be required and we don't want to remove specific classes like mutations germline in one twin only 

columns_req_filters = c('f2_mappedY', 'f3_lowDepth', 'f3_highDepth', 'f3_noReadsMapped',
                        'f4_mtr4_presentInAll', 'f4_mtr4_presentInNone', 'f4_mtr4_presentInOne', 'f5_maxVAF01', 'f5_VAFna','f6_likelyGermline_bothTwins', 'f7_FailedIndelNearby10')
twins_dt_filters[, sum_req_filters := rowSums(.SD), .SDcols = columns_req_filters] 
paste('Number of mutations that pass required filters:', dim(twins_dt_filters[sum_req_filters==0])[1]) # 1,966  

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

######################################################################################################
# SAVE FILTERED MUTATIONS TO A TXT FILE
paste('Number of mutations that passed required filters:', length(mut_included)) # 1,966
write.table(mut_included, 'Data/mutations_include_20241023.txt', quote = FALSE, col.names = F, row.names = F)

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
# ALL DONE 