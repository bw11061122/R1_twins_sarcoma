# Script to analyse the pileup (run 15/10/2024)
# 2024-10-11 - 2024-10-14
# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# dataframe created in the script: 20241015_pileup_checks.R (pileup run 14/10/2024-15/10/2024)
# NOTE that for the moment, I am investigating if there is sth wrong in how the pileup was run
# I may alter the input df once everything is resolved

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
setwd('/Users/bw18/Desktop/1SB/Analysis')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

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
# Examine which samples are being analysed 
samples = colnames(twins_dt[,c(seq(2, 330, 15))])
samples = sapply(strsplit(samples,"_"), `[`, 1)
paste('Number of samples analysed:', length(samples)) # remove the PD38is_wgs sample: not from the samples of interest

# create lists of possible samples of interest
samples_names = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb",
                  "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq",
                  "PD62341b", "PD62341h", "PD62341n", "PD62341u")
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

# Determine the number of unique mutations identified in at least one sample (ie how many mutations do I have?)
paste('Number of unique mutations identified across samples:', dim(twins_dt)[1]) # 362,814
paste('Number of unique mutations identified across samples (check):', length(twins_dt[, mut_ID] %>% unlist() %>% unique())) # 362,814

# Check that class of each column makes sense
sapply(twins_dt, class) # all fine 

######################################################################################################
# subset dataframes with DEP, MTR and VAF

# for each mutation, get data on MTR (nr mutant reads), DEP (total nr of reads), VAF
# create subsetted dataframes with this information only
cols_mtr = grep("MTR", names(twins_dt), value = TRUE)
cols_dep = grep("DEP", names(twins_dt), value = TRUE)
cols_vaf = grep("VAF", names(twins_dt), value = TRUE)

twins_mtr = twins_dt[,c('mut_ID', cols_mtr), with=FALSE]
twins_dep = twins_dt[,c('mut_ID', cols_dep), with=FALSE]
twins_vaf = twins_dt[,c('mut_ID', cols_vaf), with=FALSE]

cols_normal_mtr = paste(samples_normal, 'MTR', sep='_')
cols_normal_dep = paste(samples_normal, 'DEP', sep='_')
cols_normal_vaf = paste(samples_normal, 'VAF', sep='_')

# determine min, max and mean nr of variant reads / depth / vaf
twins_mtr[,mean_mtr := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = cols_mtr]
twins_mtr[,median_mtr := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = cols_mtr]
twins_mtr[,min_mtr := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = cols_mtr]
twins_mtr[,max_mtr := apply(.SD, 1, max), .SDcols = cols_mtr]

twins_dep[,mean_dep := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = cols_dep]
twins_dep[,median_dep := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = cols_dep]
twins_dep[,min_dep := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = cols_dep]
twins_dep[,max_dep := apply(.SD, 1, max), .SDcols = cols_dep]

twins_vaf[,mean_vaf := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = cols_vaf]
twins_vaf[,median_vaf := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = cols_vaf]
twins_vaf[,min_vaf := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = cols_vaf]
twins_vaf[,max_vaf := apply(.SD, 1, max), .SDcols = cols_vaf]

# determine this specifically in normal samples (in case of ploidy changes in the tumour)
twins_mtr[,mean_mtr_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = cols_normal_mtr]
twins_mtr[,median_mtr_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = cols_normal_mtr]
twins_mtr[,min_mtr_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = cols_normal_mtr]
twins_mtr[,max_mtr_normal := apply(.SD, 1, max), .SDcols = cols_normal_mtr]

twins_dep[,mean_dep_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = cols_normal_dep]
twins_dep[,median_dep_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = cols_normal_dep]
twins_dep[,min_dep_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = cols_normal_dep]
twins_dep[,max_dep_normal := apply(.SD, 1, max), .SDcols = cols_normal_dep]

twins_vaf[,mean_vaf_normal := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = cols_normal_vaf]
twins_vaf[,median_vaf_normal := apply(.SD, 1, function(x) median(x[x>0])), .SDcols = cols_normal_vaf]
twins_vaf[,min_vaf_normal := apply(.SD, 1, function(x) min(x[x>0])), .SDcols = cols_normal_vaf]
twins_vaf[,max_vaf_normal := apply(.SD, 1, max), .SDcols = cols_normal_vaf]

# add extra calculations to MTR (MTR >= 4 for variant to be called as present)
twins_mtr[,sum_mtr0 := rowSums(.SD==0), .SDcols = cols_mtr]
twins_mtr[,sum_mtr4 := rowSums(.SD>=4), .SDcols = cols_mtr]
twins_mtr[,sum_mtr0_normal := rowSums(.SD==0), .SDcols = cols_normal_mtr] # specific to normal samples (tumour lost chr1/18)
twins_mtr[,sum_mtr4_normal := rowSums(.SD>=4), .SDcols = cols_normal_mtr] # specific to normal samples (tumour lost chr1/18)

##########################################################################################################
# CHECKS: check the number of mutations mapped in each sample
colSums(twins_dep[,2:23] >= 10) # mapped with enough coverage 
colSums(twins_mtr[,2:23] >= 4) # mutation identified as present

##########################################################################################################
# CHECKS: PLOTS OF COVERAGE AND NR OF MUTATIONS BN SAMPLES AND ACROSS CHROMOSOMES

# plot median coverage vs number of mutations mapped (MTR >= 4)
med_cov = apply(twins_dep[,2:23], 2, median)
sample_info = data.table(cbind(cols_dep, med_cov))
sample_info[,med_cov := as.numeric(med_cov)]
sample_info[,samples := tstrsplit(samples, '_', keep = 1)]

counts = c()
for (name in cols_mtr){
  counts = c(counts, sum(twins_mtr[,..name] >= 4))
}
sample_info = cbind(sample_info, counts)
sample_info[, status := as.factor(fcase( 
  samples %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  samples %in% samples_tumour, 'tumour'
))]
sample_info[, twin := as.factor(fcase( 
  samples %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  samples %in% samples_PD63383, 'PD63383'
))]
sample_info[, sample_type := as.factor(paste(status, twin, sep = '_'))]

ggplot(sample_info, aes(x=med_cov, y=counts, col=status))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'Median coverage', y='Number of mutations (MTR >= 4)', color = 'Sample')+
  ggtitle('Median coverage vs number of mutations mapped')+
  theme_bw(base_size=12)
ggsave('Results/20241020_qc1_medCoverage_countMTR4.pdf', width = 6, height = 4)

ggplot(sample_info, aes(x=med_cov, y=counts, col=sample_type))+
  geom_point(size=2)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  labs(x = 'Median coverage', y='Number of mutations (MTR >= 4)', color = 'Sample')+
  ggtitle('Median coverage vs number of mutations mapped')+
  theme_bw(base_size=12)
ggsave('Results/20241020_qc1_medCoverage_countMTR4_colSample.pdf', width = 6, height = 4)

# COVERAGE DISTRIBUTION
# histogram for all mutations across all samples 
pdf('Results/20241020_qc0_histogramCoverage_all.pdf')
hist(twins_dep[,2:23] %>% unlist(), xlab = 'Coverage', main = 'Distribution of coverage')
dev.off()

# I want all coverage to decide which coverage boundaries to take
covs = twins_dep[,2:23] %>% unlist() 
covs2 = to_vec(for(x in covs) if(x < 80) x) 
pdf('Results/20241020_qc0_histogram_Cov80.pdf')
hist(covs2, xlab = 'Coverage', main = 'Distribution of coverage (< 80)')
dev.off() 

pdf('Results/20241020_qc0_histogram_Cov80_thresholds.pdf')
hist(covs2, xlab = 'Coverage', main = 'Distribution of coverage (< 80)')
abline(v=20, col = 'red', lwd = 2)
abline(v=50, col='red', lwd = 2)
dev.off() 

rm(covs)
rm(covs2)

# how is coverage for different mutations distributed across samples?
pdf('Results/20241020_qc0_hist_total_depth_min.pdf')
hist(twins_dep[, min_dep] %>% unlist(), xlab = 'Minimum coverage for a mutation', main = 'Distribution of minimum coverage')
dev.off()
pdf('Results/20241020_qc0_hist_total_depth_max.pdf')
hist(twins_dep[, max_dep] %>% unlist(), xlab = 'Maximum coverage for a mutation', main = 'Distribution of maximum coverage')
dev.off()
pdf('Results/20241020_qc0_hist_total_depth_mean.pdf')
hist(twins_dep[, mean_dep] %>% unlist(), xlab = 'Mean coverage for a mutation', main = 'Distribution of mean coverage')
dev.off()

# does coverage differ across chromosomes for tumour and normal samples?
twins_dep[, chr := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
twins_dep_chr_agg = data.table(aggregate(twins_dep[chr!='chrY', 2:23], by=list(Category=twins_dep[chr!='chrY', chr]), FUN=sum))
twins_dep_chr_melt = melt(twins_dep_chr_agg, id.vars = 'Category')
twins_dep_chr_melt[, chr := tstrsplit(Category, "chr", fixed=TRUE, keep=2)]
twins_dep_chr_melt[, chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                   "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                   "20", "21", "22", "X"))]
twins_dep_chr_melt[, sample_name := tstrsplit(variable, "_DEP", fixed=TRUE, keep=1)]
twins_dep_chr_melt[, sample_name := as.factor(sample_name)]
twins_dep_chr_melt[, status := as.numeric(sample_name %in% samples_tumour)] # 1 maps to tumour
twins_dep_chr_melt[, status := ifelse(status==1, "Tumour", "Normal")]
twins_dep_chr_melt[, twin := as.factor(substr(sample_name, 1, 7))]

ggplot(data=twins_dep_chr_melt, aes(x=chr, y=(value), color = status)) +
  geom_point(stat='identity')+
  theme_bw()+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'Chromosome', y = 'Total nr of reads mapped to chromosome', color = 'Status')+
  ggtitle('Reads mapped across chromosomes\n(split by status: normal vs tumour)')
ggsave('Results/20241020_qc0_dist_of_coverage_on_chr.pdf', width = 5, height = 4)
# okay so there is a bit of a difference in chromosome 1 and 18 for samples where it was likely deleted
# there are two cancer samples that for reasons unknown have very high coverage (cf everything else)
# it seems that this is general across chromosomes - not due to chromosome duplication

# median coverage per chromosome
twins_dep_chr_agg = data.table(aggregate(twins_dep[chr!='chrY', 2:23], by=list(Category=twins_dep[chr!='chrY', chr]), FUN=median))
twins_dep_chr_melt = melt(twins_dep_chr_agg, id.vars = 'Category')
twins_dep_chr_melt[, chr := tstrsplit(Category, "chr", fixed=TRUE, keep=2)]
twins_dep_chr_melt[, chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                   "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                   "20", "21", "22", "X"))]
twins_dep_chr_melt[, sample_name := tstrsplit(variable, "_DEP", fixed=TRUE, keep=1)]
twins_dep_chr_melt[, sample_name := as.factor(sample_name)]
twins_dep_chr_melt[, status := as.numeric(sample_name %in% samples_tumour)] # 1 maps to tumour
twins_dep_chr_melt[, status := ifelse(status==1, "Tumour", "Normal")]
twins_dep_chr_melt[, twin := as.factor(substr(sample_name, 1, 7))]

ggplot(data=twins_dep_chr_melt, aes(x=chr, y=(value), color = status)) +
  geom_jitter(size = 1, alpha = 0.4, position = position_jitterdodge(jitter.width = 0.2))+
  theme_bw(base_size=13)+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'Chromosome', y = 'Median coverage', color = 'Status')+
  ggtitle('Median coverage across chromosomes (split by status: normal vs tumour)')
ggsave('Results/20241020_qc0_dist_of_medianCov_on_chr.pdf', width = 9.5, height = 3.5)

twins_dep_chr_melt_ordered = twins_dep_chr_melt[order(chr, value)]
twins_dep_chr_melt_ordered[chr=='1']
twins_dep_chr_melt_ordered[chr=='18']

######################################################################################################
# MTR PLOTS

# distribution of samples with reads (# samples with mutation at MTR = 4)
pdf('Results/20241020_qc2_dist_nr_samples_w_mut.pdf')
hist(twins_mtr[,sum_mtr4] %>% unlist(), xlab = 'Number of samples with variant read at MTR=4', main = 'Distribution of # samples with mutation')
dev.off()

pdf('Results/20241020_qc2_dist_nr_samples_w_mut_exclall.pdf')
hist(twins_mtr[!sum_mtr4 %in% c(0, 22),sum_mtr4] %>% unlist(), xlab = 'Number of samples with variant read at MTR=4', main = 'Distribution of # samples with mutation (excluded shared by all)')
dev.off()

# ggplot format
count_mut = data.table(table(twins_mtr[,sum_mtr4]))
count_mut[,V1 := as.numeric(V1)]
ggplot(data=count_mut, aes(x=V1, y=N)) +
  geom_bar(stat='identity', fill = col_bar)+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = col_bar)+
  labs(x = 'Number of samples a mutation with the mutation (MTR=4)', y = 'Frequency')+
  ggtitle('Mutation counts distribution')
ggsave('Results/20241020_qc2_dist_nr_samples_w_mut_barchart.pdf', width = 6, height = 4)

count_mut = data.table(table(twins_mtr[!sum_mtr4 %in% c(0, 22),sum_mtr4]))
count_mut[,V1 := as.numeric(V1)]
ggplot(data=count_mut, aes(x=V1, y=N)) +
  geom_bar(stat='identity', fill = col_bar)+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = col_bar)+
  labs(x = 'Number of samples a mutation with the mutation (MTR=4)', y = 'Frequency')+
  ggtitle('Mutation counts distribution (excluded shared by all)')
ggsave('Results/20241020_qc2_dist_nr_samples_w_mut_barchart_exclall.pdf', width = 6, height = 4)
# BTW - interesting! most mutations present in 15 samples (not 16/17)
# Can you please look at the samples which are missing these?
# is it because some samples miss chr1 / chr18 so they won't have these mutations?

# distribution in normal samples
count_mut = data.table(table(twins_mtr[!sum_mtr4_normal %in% c(0, 12),sum_mtr4_normal]))
count_mut[,V1 := as.numeric(V1)]
ggplot(data=count_mut, aes(x=V1, y=N)) +
  geom_bar(stat='identity', fill = col_normal)+
  theme_bw(base_size = 12)+
  scale_fill_manual(values = col_bar)+
  labs(x = 'Number of samples a mutation with the mutation (MTR=4)', y = 'Frequency')+
  ggtitle('Mutation counts distribution (normal samples, excl shared)')
ggsave('Results/20241020_qc2_dist_nr_samples_w_mut_barchart_normal.pdf', width = 6, height = 4)

# distribution in tumour samples
samples_tumour_mtr = paste(samples_tumour, 'MTR', sep='_')
twins_mtr[, sum_mtr4_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
count_mut = data.table(table(twins_mtr[!sum_mtr4_tumour %in% c(0, 10),sum_mtr4_tumour]))
count_mut[,V1 := as.numeric(V1)]
ggplot(data=count_mut, aes(x=V1, y=N)) +
  geom_bar(stat='identity', fill = col_tumour)+
  theme_bw(base_size = 12)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values = col_bar)+
  labs(x = 'Number of samples a mutation with the mutation (MTR=4)', y = 'Frequency')+
  ggtitle('Mutation counts distribution (tumour samples, excl shared)')
ggsave('Results/20241020_qc2_dist_nr_samples_w_mut_barchart_tumour.pdf', width = 6, height = 4)

# exclude chromosomes 1 and 18 (missing from some samples)
twins_mtr[, chr := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
count_mut = data.table(table(twins_mtr[!sum_mtr4_tumour %in% c(0, 10) & !chr %in% c('chr1', 'chr18'), sum_mtr4_tumour]))
count_mut[,V1 := as.numeric(V1)]
ggplot(data=count_mut, aes(x=V1, y=N)) +
  geom_bar(stat='identity', fill = col_tumour)+
  theme_bw(base_size = 12)+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_fill_manual(values = col_bar)+
  labs(x = 'Number of samples a mutation with the mutation (MTR=4)', y = 'Frequency')+
  ggtitle('Mutation counts distribution (tumour, excl chr1 and chr18)')
ggsave('Results/20241020_qc2_dist_nr_samples_w_mut_barchart_tumour_excl1_18.pdf', width = 6, height = 4)

# Can we do this for each chromosome for each sample?
# for each chromosome, count the nr of mutations at MTR > 4 in each sample 
# Note that I am getting rid of all stuff mapped to chromosome Y
twins_mtr_chr_agg = data.table(aggregate(twins_mtr[chr!='chrY', 2:23], by=list(Category=twins_mtr[chr!='chrY',chr]), FUN=sum))
twins_mtr_chr_melt = melt(twins_mtr_chr_agg, id.vars = 'Category')
twins_mtr_chr_melt[, chr := tstrsplit(Category, "chr", fixed=TRUE, keep=2)]
twins_mtr_chr_melt[, chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                   "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                   "20", "21", "22", "X"))]
twins_mtr_chr_melt[, sample_name := tstrsplit(variable, "_MTR", fixed=TRUE, keep=1)]
twins_mtr_chr_melt[, sample_name := as.factor(sample_name)]
twins_mtr_chr_melt[, status := as.numeric(sample_name %in% samples_tumour)] # 1 maps to tumour
twins_mtr_chr_melt[, status := ifelse(status==1, "Tumour", "Normal")]
twins_mtr_chr_melt[, twin := as.factor(substr(sample_name, 1, 7))]

ggplot(data=twins_mtr_chr_melt, aes(x=chr, y=value, color = status)) +
  geom_point(size = 1, alpha = 0.6, position = position_jitterdodge(jitter.width = 0.2))+
  theme_bw(base_size = 12)+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'Chromosome', y = 'Number of mutations', col = 'Status')+
  ggtitle('Mutation count across chromosomes')
ggsave('Results/20241020_qc3_dist_of_nr_mut_on_chr.pdf', width = 9, height = 4)
# I guess this could help you to identify samples which lost a chromosome - very clear 
# NOTE: it tends that PD62341ap has the highest mutation burden (esp chr 6, 8, 12, 14)
# Do we think this is due to differences in coverage? > YES (see plot above)
# But then why on specific chromosomes? Amplification?

# How is number of mutations related to coverage 
setnames(twins_mtr_chr_melt, 'value', 'nr_variant_reads')
setnames(twins_dep_chr_melt, 'value', 'coverage')

twins_chr_melt = cbind(twins_mtr_chr_melt[,c(1,3:7), with=FALSE], twins_dep_chr_melt[,3])
ggplot(data=twins_chr_melt, aes(x=coverage, y=nr_variant_reads, color = chr)) +
  geom_point(stat='identity')+
  theme_bw()+
  labs(x = 'Coverage', y = 'Number of reads reporting a variant')+
  ggtitle('Nr of variant reads vs coverage')
ggsave('Results/20241020_qc4_var_reads_vs_coverage_colChr.pdf', width = 5, height = 4)

ggplot(data=twins_chr_melt, aes(x=coverage, y=nr_variant_reads, color = status)) +
  geom_point(size = 1.5, alpha = 0.8)+
  theme_bw()+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'Coverage', y = 'Number of reads reporting a variant', color = 'Status')+
  ggtitle('Nr of variant reads vs coverage')
ggsave('Results/20241020_qc4_var_reads_vs_coverage_colStatus.pdf', width = 5, height = 4)

ggplot(data=twins_chr_melt, aes(x=coverage, y=nr_variant_reads, color = twin)) +
  geom_point(stat='identity')+
  theme_bw()+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  labs(x = 'Coverage', y = 'Number of reads reporting a variant')+
  ggtitle('Nr of variant reads vs coverage')
ggsave('Results/20241020_qc4_var_reads_vs_coverage_colTwin.pdf', width = 5, height = 4)

######################################################################################################
# EXTRA COVERAGE QC PLOTS 

# sum up coverage across chromosomes
twins_dep[, chr := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
twins_dep[, chr := as.factor(chr)]

# aggregate coverage by chromosomes 
twins_cov_agg = data.table(aggregate(twins_dep[chr!='chrY', 2:23], by=list(Category=twins_dep[chr!='chrY', chr]), FUN=sum))
twins_cov_melt = melt(twins_cov_agg, id.vars = 'Category')
twins_cov_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
twins_cov_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
twins_cov_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
twins_cov_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
twins_cov_melt[, chr := tstrsplit(Category, "chr", fixed=TRUE, keep=2)]
twins_cov_melt[, chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                               "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                               "20", "21", "22", "X"))]

chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", 
        "20", "21", "22", "X")
chr_lengths = as.numeric(c(248956422, 242193529, 198295559, 190214555,
                           181528259, 170805979, 159345973, 145138636,
                           139394717, 133797422, 135086622, 133275309,
                           114364328, 107043718, 101991189, 90338345,
                           83257441, 80373285, 58617616, 64444167,
                           46709983, 50818468, 150040895))

chr_lengths_dt = data.table(cbind(chr, chr_lengths))
chr_lengths_dt[, chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                               "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                               "20", "21", "22", "X"))]
twins_cov_melt = merge(twins_cov_melt, chr_lengths_dt, by = 'chr')

twins_cov_melt[, chr_lengths := as.numeric(chr_lengths)]
twins_cov_melt[, cov_length := value / chr_lengths]

ggplot(twins_cov_melt, aes(x = chr, y = value, color = twin))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  labs(x = 'Chromosome', y = 'Number of reads mapped to the chromosome')+
  ggtitle('Nr of mapped reads across chr')

ggplot(twins_cov_melt, aes(x = chr, y = value, color = status))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'Chromosome', y = 'Number of reads mapped to the chromosome')+
  ggtitle('Nr of mapped reads across chr')

ggplot(twins_cov_melt, aes(x = chr, y = value, color = sample_type))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  labs(x = 'Chromosome', y = 'Number of reads mapped to the chromosome')+
  ggtitle('Nr of mapped reads across chr')

# Coverage ie number of reads per chr length 
ggplot(twins_cov_melt, aes(x = chr, y = cov_length, color = twin))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  labs(x = 'Chromosome', y = 'Nr of mapped reads per bp')+
  ggtitle('Coverage across chr')

ggplot(twins_cov_melt, aes(x = chr, y = cov_length, color = status))+
  geom_point(size=1, position=position_dodge(.75))+
  scale_color_manual(values = c(col_normal, col_tumour))+
  theme_bw()+
  labs(x = 'Chromosome', y = 'Nr of mapped reads per bp')+
  ggtitle('Coverage across chr')

ggplot(twins_cov_melt, aes(x = chr, y = cov_length, color = sample_type))+
  geom_point(size=1, position=position_dodge(1))+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  theme_bw()+
  labs(x = 'Chromosome', y = 'Nr of mapped reads per bp')+
  ggtitle('Coverage across chr')

# which tumour samples have reduced coverage over chr1 and chr18:
twins_cov_sorted = twins_cov_melt[order(chr,value)] # order by chr and then value (nr mapped reads)
twins_cov_sorted[chr=='1']
# samples with lower coverage: PD63383ap, PD63383aq, PD62341ae, PD62341am, PD62341b
twins_cov_sorted[chr=='18']
# samples with lower coverage: PD63383ap, PD63383aq, PD62341ae, PD62341am, PD62341b
# good: these are the same samples! 

# samples with elevated coverage (across most chr, esp 2/3/5/6/7/8)
twins_cov_sorted[chr=='2'] # PD62341u_DEP, PD62341ap_DEP - for all chromosomes 

# aggregate by median
twins_mcov_agg = data.table(aggregate(twins_dep[chr!='chrY', 2:23], by=list(Category=twins_dep[chr!='chrY', chr]), FUN=median))
twins_mcov_melt = melt(twins_mcov_agg, id.vars = 'Category')
twins_mcov_melt[, sample := tstrsplit(variable, '_', fixed=TRUE, keep = 1)]
twins_mcov_melt[, status := as.factor(fcase( 
  sample %in% samples_normal, 'normal', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_tumour, 'tumour'
))]
twins_mcov_melt[, twin := as.factor(fcase( 
  sample %in% samples_PD62341, 'PD62341', # differs from 0.5 so biased and maybe we don't want it
  sample %in% samples_PD63383, 'PD63383'
))]
twins_mcov_melt[, sample_type := as.factor(paste(status, twin, sep = '_'))]
twins_mcov_melt[, chr := tstrsplit(Category, "chr", fixed=TRUE, keep=2)]
twins_mcov_melt[, chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                                "20", "21", "22", "X"))]
twins_mcov_melt = merge(twins_mcov_melt, chr_lengths_dt, by = 'chr')
twins_mcov_melt[, chr_lengths := as.numeric(chr_lengths)]

ggplot(twins_mcov_melt, aes(x = chr, y = value, color = twin))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(col_PD62341, col_PD63383))+
  labs(x = 'Chromosome', y = 'Median coverage')+
  ggtitle('Median coverage across chromosomes')

ggplot(twins_mcov_melt, aes(x = chr, y = value, color = status))+
  geom_point(size = 1, position = position_dodge(.8))+
  theme_bw(base_size=14)+
  scale_color_manual(values = c(col_normal, col_tumour))+
  labs(x = 'Chromosome', y = 'Median coverage')+
  ggtitle('Median coverage across chromosomes')

ggplot(twins_mcov_melt, aes(x = chr, y = value, color = sample_type))+
  geom_point(size = 1, alpha = 0.8, position = position_dodge(.5))+
  theme_bw(base_size = 14)+
  scale_color_manual(values = c(col_normal_PD62341, col_normal_PD63383, col_tumour_PD62341, col_tumour_PD63383))+
  labs(x = 'Chromosome', y = 'Number of reads mapped to the chromosome')+
  ggtitle('Median coverage across chromosomes')

twins_mcov_sorted = twins_mcov_melt[order(chr,value)] # order by chr and then value (nr mapped reads)
twins_mcov_sorted[chr=='1'] # low: PD63383aq, ap; PD62341ae, am, b
twins_mcov_sorted[chr=='18'] # low: PD63383aq, ap; PD62341ae, am, b
twins_mcov_sorted[chr=='2'] # high: PD63383: none, PD62341u, ap, aj, ag
twins_mcov_sorted[chr=='3'] # high: PD63383: none, PD62341u, ap, aj, ag

# from clinical data, low level duplication of chr 6, 8, 12, 13, 14, 16, 20, 21
twins_mcov_sorted[chr=='6'] # high: PD63383: none, PD62341u, ap, aj, ag, ak, also aa
twins_mcov_sorted[chr=='8'] # high: PD63383: none, PD62341u, ap, aj, ag, also aa
twins_mcov_sorted[chr=='12'] # high: PD63383: none, PD62341u, ap, aj, ak, also aa
twins_mcov_sorted[chr=='13'] # high: PD63383: none, PD62341u, ap, ak, ag
twins_mcov_sorted[chr=='14'] # high: PD63383: none, PD62341u, ap, aj, ak, also PD62341aa and q (contamination?)
twins_mcov_sorted[chr=='16'] # high: PD63383: none, PD62341u, ap, aj, also PD62341aa and q (contamination?) # the same pattern 
twins_mcov_sorted[chr=='20'] # high: PD63383: none, PD62341u, ap, aj, also PD62341aa and q (contamination?) # the same pattern
twins_mcov_sorted[chr=='21'] # high: PD63383: none, PD62341u, ap, aj, also PD62341aa and q (contamination?) # the same pattern

# check median coverage in each sample
apply(twins_dep[,2:23], 2, FUN=median)

######################################################################################################
######################################################################################################
# QC PLOTS ON FILTERED MUTATIONS

# read list of mutations IDs that passed all filters
mut_include = read.table('Data/mutations_include_20241018.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(mut_include)) # 2,559
twins_filtered_dt = twins_dt[mut_ID %in% mut_include]

######################################################################################################
# Compare classes of mutations pre- and post-filtering 

twins_dt[, mut_cat := paste(Ref, Alt, sep='>')]
twins_dt[, mut_cat := as.factor(fcase(
  mut_cat == 'A>T', 'T>A',
  mut_cat == 'A>C', 'T>G',
  mut_cat == 'G>C', 'C>G',
  mut_cat == 'G>T', 'C>A',
  mut_cat == 'A>G', 'T>C',
  mut_cat == 'G>A', 'C>T'
))]
twins_dt[, qc := as.factor(fcase(
  mut_ID %in% mut_include, 'passed',
  !mut_ID %in% mut_include, 'failed'
))]

mut_counts = data.table(sort(table(twins_dt[,mut_cat])), decreasing=TRUE)
mut_counts[,V1:=as.factor(V1)]

# plot distribution of mutations
ggplot(data=mut_counts, aes(x=fct_inorder(V1), y=N)) +
  geom_bar(stat='identity', fill = col_bar)+
  theme_bw(base_size = 12)+
  labs(x = 'Mutation type', y = 'Frequency', title = 'Mutation types in all samples')

# plot distribution of mutations of failed and passed mutations
mut_counts_passed = data.table(sort(table(twins_dt[qc == 'passed',mut_cat])))
mut_counts_failed = data.table(sort(table(twins_dt[qc == 'failed',mut_cat])))

mut_counts_passed[, freq_norm := N/ length(mut_include)]
mut_counts_failed[, freq_norm := N/ (386412 - length(mut_include))]

mut_counts_passed[, QC := 'passed']
mut_counts_failed[, QC := 'failed']
mut_counts_qc = data.table(rbind(mut_counts_passed, mut_counts_failed))
mut_counts_qc[, V1 := as.factor(V1)]

ggplot(data=mut_counts_qc, aes(x=fct_inorder(V1), y=freq_norm, fill=QC)) +
  geom_bar(stat='identity', position = 'dodge')+
  scale_fill_manual(values = c('darkblue','lightblue'))+
  theme_bw(base_size = 12)+
  labs(x = 'Mutation type', y = 'Fraction of all mutations in the category', title = 'Mutation types in all samples')
ggsave('Results/20241020_qc5_mut_filters_types.pdf', width = 6, height = 4.5)

######################################################################################################
# Coverage vs VAF for each sample

twins_filtered_dt[, loss := as.factor(fcase( 
  Chrom %in% c('chr1', 'chr18'), 'loss in tumour', # chr1 and chr18 segments lost in tumour samples
  !Chrom %in% c('chr1', 'chr18'), 'normal ploidy'
))]

# plot coverage vs VAF for each sample, color by ploidy status in the tumour
for (sample in samples_names){
  DEP = paste0(sample, '_DEP')
  VAF = paste0(sample, '_VAF')
  df = twins_filtered_dt[, c('mut_ID', VAF, DEP, 'loss'), with=FALSE]
  df[, chr := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
  ggplot(df, aes(x=get(names(df)[3]), y=get(names(df)[2]), color = loss))+
    geom_point()+
    theme_bw()+
    scale_color_manual(values = c(col_normal_PD62341, col_tumour_PD63383))+
    labs(x = 'Total number of reads covering a base', y = 'Fraction of reads reporting a variant allele', color = 'Chromosome')+
    ggtitle(glue('Nr of variant reads vs coverage: {sample}'))
  ggsave(glue('Results/20241020_p0_vaf_vs_cov_{sample}.pdf'), width=6, height=4.5)
}

# PLOT: coverage vs VAF (for each sample separately), color by duplication status (chr1 and chr18)
# Can we aggregate all tumour samples together?
samples_tumour_DEP = paste(samples_tumour, 'DEP', sep='_')
samples_tumour_VAF = paste(samples_tumour, 'VAF', sep='_')
twins_filtered_dt[, sum_dep_tumour := rowSums(.SD), .SDcols = samples_tumour_DEP]
twins_filtered_dt[, mean_vaf_tumour := apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = samples_tumour_VAF]
twins_filtered_dt[, chr := tstrsplit(mut_ID, '_', fixed=TRUE, keep=1)]
ggplot(twins_filtered_dt, aes(x=sum_dep_tumour, y=mean_vaf_tumour, color = loss))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c(col_normal_PD62341, col_tumour_PD63383))+
  labs(x = 'Total number of reads covering a base', y = 'Fraction of reads reporting a variant allele', color = 'Chromosome')+
  ggtitle(glue('Nr of variant reads vs coverage: all tumour samples aggregated'))
ggsave(glue('Results/20241020_p0_vaf_vs_cov_allTumour.pdf'), width=6, height=3.5)





