# Test scripts to get used to CaVEMan output
# 2024-10-01 
# Barbara Walkowiak bw18

# INPUT: dataframe with mutation calls (from CaVEMan) for tumour and normal samples
# dataframe received from Henry Lee-Six 01/10/2024 

# In this analysis, I want to:
## carry out basic checks on the dataset
## establish which twin the tumour first arose in 

# Load needed libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(glue)

# Read the file (received from Henry)
twins_sample = read.table('Data/twins_caveman_unmatched_20241001.txt')
paste('Number of mutations detected in any given sample:', dim(twins_sample)[1])
# rows: each row corresponds to a mutation detected in a given sample sample
# column: see explanation below 

# Specify column names to be more descriptive
setnames(twins_sample, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6','V7', 'V8', 'V9', 'V10', 'V11'),
         c('mutation', 'position', 'id', 'ref', 'var', 'filt', 'filt2', 'info', 'filters', 'normal', 'tumour'))
# Ensure data.table format 
twins_sample = data.table(twins_sample)

# COLUMN INFO: EXPLANATION
# DP: total depth 
# MP: sum of CaVEMan somatic genotype probabilities
# GP: Sum of CaVEMan germline genotype probabilities
# TG: Most probable genotype as called by CaVEMan
# TP: Probability of the most probable genotype as called by CaVEMan
# SG: 2nd most probable genotype as called by CaVEMan
# SP: Probability of the 2nd most probable genotype as called by CaVEMan
# DS: DBSnp ID of known SNP
# SNP: Position matches a dbSNP entry using the supplied bed file
# ASRD: A soft flag median (read length adjusted) alignment score of reads showing the variant allele
# CLPM: A soft flag median number of soft clipped bases in variant supporting reads
# ASMD: A soft flag median alignment score of reads showing the variant allele
# VD: Vagrent Default Annotation
# VW: Vagrent Most Deleterious Annotation
# VT: Variant type based on the Vagrent Default Annotation
# VC: Variant consequence based on the Vagrent Default Annotation

# FORMAT FORMAT
# GT: genotype 
# FAZ / FCZ / FGZ / FTZ: reads presenting A / C / G / T, forward strand
# RAZ / RCZ / RGZ / RTZ: reads presenting A / C / G / T, reverse strand
# PM: proportion of mutant allele 

# Select only specific columns of interest
cols <- c('mutation', 'position', 'ref', 'var', 'info', 'tumour')
twins_df = twins_sample[, cols, with=FALSE]

# Extract data on the sample (ie whc=ich twin it came from)
twins_df[, sample:= tstrsplit(mutation, ".", fixed=TRUE, keep=1)]
twins_df[, sample:=as.factor(sample)]

# Determine how many different samples were analysed in the dataset
paste('Number of different samples analysed:', length(twins_df[, sample] %>% unlist() %>% unique()))
table(twins_df[, sample]) # how many mutations are present in each of the samples

# Select samples for which we know which tissue / tumour they correspond to 
samples_known = list('PD62341ae', 'PD62341am', 'PD62341aa', 'PD62341q', 'PD62341v',
                     'PD63383ap', 'PD63383aq', 'PD63383bb', 'PD63383u', 'PD63383w')
twins_df = twins_df[sample %in% samples_known]
paste("Number of mutations in samples of interest:", dim(twins_df)[1])

# Add labels so we know which samples we are looking at
twins_df[,twin:=substr(sample, 1, 7)] # identify the twin 
twins_df[,twin:=as.factor(twin)]

twins_df[,sample_name:=as.factor(fcase( # add what sample that was (based on data from Henry)
  sample=='PD62341ae', 'tumour1_PD62341',
  sample=='PD62341am', 'tumour2_PD62341',
  sample=='PD62341aa', 'skin_PD62341',
  sample=='PD62341q', 'pancreas_PD62341',
  sample=='PD62341v', 'spleen_PD62341',
  sample=='PD63383ap', 'tumour1_PD63383',
  sample=='PD63383aq', 'tumour2_PD63383',
  sample=='PD63383bb', 'skin_PD63383',
  sample=='PD63383u', 'pancreas_PD63383',
  sample=='PD63383w', 'spleen_PD63383'
))] 

# add whether the sample is tumour or normal
twins_df[,sample_cat:=as.factor(fcase( # whether the sample is tumour or normal
  sample %in% c('PD62341ae', 'PD62341am','PD63383ap','PD63383aq'), 'tumour',
  sample %in% c('PD62341aa', 'PD62341q', 'PD62341v', 'PD63383bb', 'PD63383u', 'PD63383w'), 'normal'))]

# identify the actual mutation (this requires chromosome + position + change)
twins_df[, chr:=tstrsplit(mutation, ":", fixed=TRUE, keep=2)]
twins_df[, mut:=paste(chr, position, ref, var, sep = '_')]
twins_df[, mut_type:=factor(paste(ref, var, sep = '>'))]

#################################################################################################
#################################################################################################
# CHECKS: MUTATION TYPE AND MUTATION DISTRIBUTION ON CHROMOSOMES

# We can very quickly check the nr of mutations per chromosome and the types of mutations
# C>T and G>A are effectively the same mutation type, so it's good we have equal numbers of those
# Note that currently, we have all mutations, not just somatic / germline etc. 
mut_counts = data.table(sort(table(twins_df[,mut_type])), decreasing=TRUE)
mut_counts[,V1:=as.factor(V1)]
ggplot(data=mut_counts, aes(x=fct_inorder(V1), y=N)) +
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x = 'Mutation type', y = 'Frequency', title = 'Mutation types')
ggsave('Results/p1_mut_counts_by_type.pdf', width = 5, height = 4)

# Distribution of different mutation types between tumour vs normal samples
mut_counts_tumour = data.table(sort(table(twins_df[sample_cat=='tumour',mut_type])))
mut_counts_normal = data.table(sort(table(twins_df[sample_cat=='normal',mut_type])))
setnames(mut_counts_tumour, 'N', 'tumour')
setnames(mut_counts_normal, 'N', 'normal')
mut_counts_merged = merge(mut_counts_tumour, mut_counts_normal)
mut_counts_merged[,V1:=as.factor(V1)]
mut_counts_merged[,ratio:=(tumour/4)/(normal/6)]
ggplot(data=mut_counts_merged, aes(x=fct_inorder(V1), y=ratio)) +
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x = 'Mutation type', y = 'Frequency', title = 'Mutation types in tumour vs normal')
ggsave('Results/p1_mut_counts_by_type_tumour_vs_normal.pdf', width = 5, height = 4)
# seems like tumour samples don't have more mutations than normal
# but I guess this is because we have a TON that are germline cf everything else 
# we should check this but for mutations at lower VAF 

# We can also plot mutations by chromosome 
chr_counts = data.table(table(twins_df[,chr]))
chr_counts[,chr := tstrsplit(V1, "chr", fixed=TRUE, keep=2)]
chr_counts[,chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                             "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                             "20", "21", "22", "X", "Y"))]
ggplot(data=chr_counts, aes(x=chr, y=N)) +
  geom_bar(stat='identity')+
  theme_bw()+
  theme(text = element_text(size = 10))+
  labs(x = 'Chromosome', y = 'Frequency', title = 'Mutations by chromosome')
ggsave('Results/p1_mut_counts_by_chr.pdf', width = 5, height = 4)

chr_list = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
             "11", "12", "13", "14", "15", "16", "17", "18", "19", 
             "20", "21", "22", "X", "Y")
chr_lengths = as.numeric(c(248956422, 242193529, 198295559, 190214555,
                181528259, 170805979, 159345973, 145138636,
                139394717, 133797422, 135086622, 133275309,
                114364328, 107043718, 101991189, 90338345,
                83257441, 80373285, 58617616, 64444167,
                46709983, 50818468, 150040895, 57227415))

setnames(chr_counts, "chr", "chr_list")
chr_lengths_dt = data.table(cbind(chr_list, chr_lengths))
chr_dt = merge(chr_lengths_dt, chr_counts)
chr_dt[,chr_list := factor(chr_list, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                          "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                          "20", "21", "22", "X", "Y"))]

chr_dt[,chr_lengths := as.numeric(chr_lengths)]
chr_dt[,mut_density := N / chr_lengths * 10^6]
ggplot(data=chr_dt, aes(x=chr_list, y=mut_density)) +
  geom_bar(stat='identity')+
  theme_bw()+
  theme(text = element_text(size = 10))+
  labs(x = 'Chromosome', y = 'Mutation density (per Mb)', title = 'Mutations by chromosome')
ggsave('Results/p1_mut_counts_density.pdf', width = 5, height = 4)

# Mutation density - normal vs tumour samples
# Are mutations in tumours vs normal samples differently distributed across chromosomes?
chr_counts_tumour = data.table(table(twins_df[sample_cat=='tumour',chr]))
chr_counts_normal = data.table(table(twins_df[sample_cat=='normal',chr]))
setnames(chr_counts_tumour, 'V1', 'chr')
setnames(chr_counts_normal, 'V1', 'chr')

chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", 
        "20", "21", "22", "X", "Y")
chr_lengths = as.numeric(c(248956422, 242193529, 198295559, 190214555,
                           181528259, 170805979, 159345973, 145138636,
                           139394717, 133797422, 135086622, 133275309,
                           114364328, 107043718, 101991189, 90338345,
                           83257441, 80373285, 58617616, 64444167,
                           46709983, 50818468, 150040895, 57227415))
chr_data = data.table(cbind(chr, chr_lengths))

chr_counts_tumour[,chr := tstrsplit(chr, 'chr', keep=2)]
chr_counts_normal[,chr := tstrsplit(chr, 'chr', keep=2)]

chr_tumour = merge(chr_counts_tumour, chr_data)
chr_normal = merge(chr_counts_normal, chr_data)

chr_tumour[,chr_lengths := as.numeric(chr_lengths)]
chr_tumour[,mut_density_tumour := N / (chr_lengths * 4) * 10^6] # normalise to nr of samples = 4
chr_normal[,chr_lengths := as.numeric(chr_lengths)]
chr_normal[,mut_density_normal := N / (chr_lengths * 6) * 10^6] # normalise to nr of samples = 1

chr_merged = cbind(chr_tumour[,c("chr", "chr_lengths", "mut_density_tumour"), with=FALSE], chr_normal[,"mut_density_normal"])
chr_merged[,mut_ratio := mut_density_normal / mut_density_tumour]
chr_merged[,chr := factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                          "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                                          "20", "21", "22", "X", "Y"))]

ggplot(data=chr_merged, aes(x=chr, y=mut_ratio)) +
  geom_bar(stat='identity')+
  theme_bw()+
  theme(text = element_text(size = 10))+
  labs(x = 'Chromosome', y = 'Mutation density ratio (normal / tumour)', title = 'Ratio of mutation density')
ggsave('Results/p1_mut_counts_density_normal_vs_tumour.pdf', width = 5, height = 4)
# these results are consistent with the tumour losing chr1 and chr18 
# there doesn't seem to be anything doing on in chr 6 / any other chromosome 

# plot normal and tumour separately
ggplot(data=chr_merged, aes(x=chr, y=mut_density_tumour)) +
  geom_bar(stat='identity')+
  theme_bw()+
  theme(text = element_text(size = 10))+
  labs(x = 'Chromosome', y = 'Mutation density (per Mb)', title = 'Mutations by chromosome in tumour samples')
ggsave('Results/p1_mut_counts_density_tumour.pdf', width = 5, height = 4)

ggplot(data=chr_merged, aes(x=chr, y=mut_density_normal)) +
  geom_bar(stat='identity')+
  theme_bw()+
  theme(text = element_text(size = 10))+
  labs(x = 'Chromosome', y = 'Mutation density (per Mb)', title = 'Mutations by chromosome in normal samples')
ggsave('Results/p1_mut_counts_density_normal.pdf', width = 5, height = 4)

#################################################################################################
#################################################################################################

#################################################################################################
#################################################################################################
# CHECKS: VAF AND TOTAL COVERAGE DISTRIBUTION

# Extract VAF (proportion of reads with the variant)
# NOTE: most mutations are GERMLINE so we expect VAF to be distributed around 0.5
twins_df[,vaf:=tstrsplit(tumour, ':', keep=10)] 
twins_df[,vaf:=as.numeric(vaf)]

# Extract coverage data
twins_df[, total_depth := tstrsplit(info, 'DP=', keep=2)]
twins_df[, total_depth := tstrsplit(total_depth, ';MP', keep=1)]
twins_df[, total_depth := as.numeric(total_depth)]

# Determine the number of reads with the variant
twins_df[, variant_depth := as.integer(total_depth * vaf)]

# visualize the distribution of coverage (total depth)
pdf('Results/p2_hist_total_depth.pdf')
hist(twins_df[, total_depth], main = 'Distribution of coverage (total depth)', xlab = 'Total depth')
dev.off()
paste('Median coverage:', median(twins_df[, total_depth]))
paste('Number of mutations with total coverage below 50:', dim(twins_df[total_depth<50])[1])

pdf('Results/p2_hist_variant_depth.pdf')
hist(twins_df[,variant_depth], , main = 'Distribution of number of variant reads', xlab = 'Number of variant reads')
dev.off()
paste('Median nr of reads reporting a variant:', median(twins_df[, variant_depth]))
paste('Minimum number of reads reporting a variant:', min(twins_df[, variant_depth])) 
# we require 4 reads to call a variant - quite stringent (less likely to be an artefact, 
# although the error rate is quite high. 
# be careful of cases where you have 3 reads in a sample - are we confident a mutation is absent in this case?

# visualize the distribution of VAF values
pdf('Results/p2_hist_vaf.pdf')
hist(twins_df[,vaf], xlab = 'VAF', main = 'Histogram of VAF (entire dataset)')
dev.off()
# VAF if distributed around 0.5, which makes sense given that we are mostly looking at germline mutations
# For somatic mutations, we would expect those to be present at VAF < 0.5 (markedly < 0.5)
# Interestingly, you do get mutations with VAF 1 (across every chromosome)

# Compare VAF distribution between tumour and normal
pdf('Results/p2_hist_vaf_tumour_vs_normal.pdf')
par(mfrow = c(2,1))
hist(twins_df[sample_cat=='tumour',vaf], xlab = 'VAF', main = 'Histogram of VAF (tumour samples)')
hist(twins_df[sample_cat=='normal',vaf], xlab = 'VAF', main = 'Histogram of VAF (normal samples)')
dev.off()

pdf('Results/p2_hist_XY_checking_sex.pdf')
par(mfrow = c(2,1))
hist(twins_df[chr=='chrX', vaf], xlab = 'VAF (reads mapped to chr X)', main = 'VAF (mutations mapped to X)')
hist(twins_df[chr=='chrY', vaf], xlab = 'VAF (reads mapped to chr Y)',  main = 'VAF (mutations mapped to Y)')
dev.off()

#################################################################################################
# I think it is useful to filter out mutations which have consistently too high / too low coverage 
# these tend to correspond to mapping issues and may not be as reliable
twins_df_sub = twins_df[,c('mut', 'total_depth'), with=FALSE]
twins_df_sub[order(mut)]

max_depth = setDT(twins_df_sub)[, .SD[which.max(total_depth)], by=mut]
min_depth = setDT(twins_df_sub)[, .SD[which.min(total_depth)], by=mut]
setnames(max_depth, 'total_depth', 'max_depth')
setnames(min_depth, 'total_depth', 'min_depth')
depths_merged = merge(max_depth, min_depth)
depths_merged[,diff:=max_depth-min_depth]

pdf('Results/p2_hist_depth_difference.pdf')
hist(depths_merged[,diff], xlab = 'Max VAF - Min VAF', main = 'Difference between max and min coverage for a mutation')
dev.off()

paste('Median difference in max to min coverage for a mutation:', median(depths_merged[,diff]))

# Should we filter out mutations which have too high min coverage / too low max coverage?
# Perhaps we can remove those from our list of mutations to examine 
# NOTE: I don't know 100% which VAF threshold would be appropriate 
paste('Number of mutations with max coverage < 30:', dim(depths_merged[max_depth<30])[1])
paste('Number of mutations with min coverage > 100:', dim(depths_merged[min_depth>100])[1])

# Checked some of those on JBrowse and a lot of those look from not great to very dodgy 
depths_merged[min_depth>100, mut]
depths_merged[max_depth<30, mut]

mutations_exclude_vaf100 = depths_merged[min_depth>100, mut] %>% unlist() 
mutations_exclude_vaf30 = depths_merged[max_depth<30, mut] %>% unlist() 

#################################################################################################
#################################################################################################
# CHECKS: NUMBER OF MUTATIONS CALLED
# Identify patterns of mutation distribution across different samples 

# Determine how many mutations have been called (total and unique)
paste('Number of all mutations identified:', dim(twins_df)[1]) # 3,251,995
paste('Number of unique mutations identified:', length(twins_df[,mut] %>% unlist() %>% unique())) # 355,453

# Check that a specific mutation is only present once in a given sample (no duplicates w/in samples)
sample_names = c('tumour1_PD62341', 'tumour1_PD63383', 'tumour2_PD62341', 'tumour2_PD62341',
                    'skin_PD62341', 'pancreas_PD62341', 'spleen_PD62341',
                    'skin_PD63383', 'pancreas_PD63383', 'spleen_PD63383')

for (name in sample_names){
  if (length(twins_df[sample_name==name,mut] %>% unlist()) == length(twins_df[sample_name==name,mut] %>% unlist() %>% unique()))
  {
    print('No duplicate mutations')  
  }
  else
  {
    paste('Duplicate mutations identified in sample:', name)
  }
}

# In each sample, count the number of unique mutations 
mut_counts_by_sample = data.table(table(twins_df[,sample_name]))
mut_counts_by_sample

ggplot(data=mut_counts_by_sample, aes(x=fct_inorder(V1), y=N)) +
  geom_bar(stat='identity')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = 'Sample', y = 'Frequency', title = 'Mutation counts in samples')
ggsave('Results/p3_mutation_count_by_sample.pdf', width=6, height=4.5)

# Determine the distribution of mutations across samples
counts_mut = data.table(sort(table(twins_df[,mut]),decreasing=TRUE))
counts_mut2 = data.table(table(counts_mut[,N]))
ggplot(data=counts_mut2, aes(x=fct_inorder(V1), y=N)) +
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x = 'Number of samples the mutation is present in', y = 'Frequency', title = 'Mutation counts distribution')
ggsave('Results/p3_dist_nr_samples_with_mutation.pdf', width=6, height=4.5)
# most mutations are  present in all samples
# that makes sense if they are mostly germline / de novo in the zygote and present at high VAF - we will remove those

#################################################################################################
#################################################################################################
# PHYLOGENY RECONSTRUCTION
# ANALYSIS OF MUTATION STATUS 

# Aim: a table with mutation and status in each sample (we could do this by VAF and then convert to 0/1)
# The rows corresponds to unique mutations identified in any of the samples examined (10 samples total)
# Columns correspond to each of the 10 samples 
# values = VAF in the sample (> 0 means mutation is present)

# Create a dt with 1 column of all unique mutations across the whole dataset (10 samples)
all_mutations = data.table(twins_df[,mut] %>% unlist() %>% unique())
setnames(all_mutations, "V1", "mut")

# For each of those mutations, we want to know the vaf of this mutation in each sample 
sample_names = c('tumour1_PD62341', 'tumour2_PD62341',
                'skin_PD62341', 'pancreas_PD62341', 'spleen_PD62341',
                'tumour1_PD63383', 'tumour2_PD63383',
                'skin_PD63383', 'pancreas_PD63383', 'spleen_PD63383')

for (name in sample_names){
  
  # subset the df
  cols = c('mut', 'vaf')
  df = twins_df[sample_name==name,cols,with=FALSE] # subset the df 
  all_mutations = merge(all_mutations, df, all=TRUE)
  setnames(all_mutations, "vaf", glue("{name}"))
  
}

# replace NA with 0 
for (j in names(all_mutations)){
  
  set(all_mutations, which(is.na(all_mutations[[j]])), j, 0) 

}

# Add columns on the max, mean and min VAF
all_mutations[,max_vaf:= apply(.SD, 1, max), .SDcols = 2:11]
all_mutations[,mean_vaf:= apply(.SD, 1, mean), .SDcols = 2:11]
all_mutations[,min_vaf:= apply(.SD, 1, min), .SDcols = 2:11]

# min and mean VAF for samples where the variant read was identified
all_mutations[,mean_vaf_present:= apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = 2:11]
all_mutations[,min_vaf_present:= apply(.SD, 1, function(x) min(x[x>0])), .SDcols = 2:11]

# check the distribution of VAF values
hist(all_mutations[,max_vaf], xlab = 'Maximum VAF', main = 'Histogram of max VAF values')
hist(all_mutations[,mean_vaf], xlab = 'Mean VAF', main = 'Histogram of mean VAF values')
hist(all_mutations[,min_vaf], xlab = 'Minimum VAF', main = 'Histogram of min VAF values')

hist(all_mutations[,mean_vaf_present], xlab = 'Mean VAF', main = 'Histogram of mean VAF (in samples with the variant)')
hist(all_mutations[,min_vaf_present], xlab = 'Minimum VAF', main = 'Histogram of min VAF (in samples with the variant)')
# most of the minimum values are also around 0.4 - very likely to be present at true VAF 0.4

# determine the number of samples the mutation is present in
all_mutations[,sum_samples:= rowSums(.SD>0), .SDcols = 2:11]
all_mutations[,sum_normal:= rowSums(.SD>0), .SDcols = c(4:6, 9:11)]
all_mutations[,sum_tumour:= rowSums(.SD>0), .SDcols = c(2:3, 7:8)]
all_mutations[,sum_PD62341:= rowSums(.SD>0), .SDcols = c(2:6)]
all_mutations[,sum_PD63383:= rowSums(.SD>0), .SDcols = c(7:11)]
all_mutations[,sum_PD62341_tumour:= rowSums(.SD>0), .SDcols = c(2:3)]
all_mutations[,sum_PD63383_tumour:= rowSums(.SD>0), .SDcols = c(7:8)]
all_mutations[,sum_PD62341_normal:= rowSums(.SD>0), .SDcols = c(4:6)]
all_mutations[,sum_PD63383_normal:= rowSums(.SD>0), .SDcols = c(9:11)]

# check that this matches the distribution identified previously 
hist(all_mutations[,sum_samples], xlab = 'Number of samples with mutation', main = 'Distribution of # samples with mutation')

# check the distribution in normal and tumour samples 
hist(all_mutations[,sum_normal], xlab = 'Number of normal samples with the mutation', main = 'Normal samples')
hist(all_mutations[,sum_tumour], xlab = 'Number of tumour samples with the mutation', main = 'Tumour samples')

# compare distribution between twins 
hist(all_mutations[,sum_PD62341], xlab = 'Number of PD62341 samples with the mutation', main = 'PD62341')
hist(all_mutations[,sum_PD63383], xlab = 'Number of PD63383 samples with the mutation', main = 'PD63383')

# BINARY DATAFRAME
# Create a new data table, with mutation status encoded as binary (0 if absent, 1 if present)
all_mutations_binary = data.table(all_mutations)

# if the proportion of reads is greater than 0, convert to 1 
# only convert VAF values 
for (j in names(all_mutations_binary)[2:11]){
  set(all_mutations_binary, which(is.numeric(all_mutations_binary[[j]]) & all_mutations_binary[[j]] > 0), j, 1)
}

##################################################################################################################
##################################################################################################################
# SOME ATTEMPTS AT FILTERING

# Which mutations can we filter out?
# We know that most mutations of the ones (355,435) we see are germline, and only some are somatic

# 1 remove mutations present in all or only one sample - these are not informative  
all_mutations_filt = all_mutations[!sum_samples %in% c(1,10),] 
paste('Number of mutations present in 2-9 samples:', dim(all_mutations_filt)[1]) # should be 61,806

# 2 remove mutations with suspicious (too low / high coverage)
all_mutations_filt = all_mutations_filt[!mut %in% c(mutations_exclude_vaf100, mutations_exclude_vaf30)]
paste('Number of mutations which pass coverage filter:', dim(all_mutations_filt)[1]) # should be 61,345

# 3 retain mutations present at min VAF < 0.5 (very likely that in this case it is a germline mutation)
all_mutations_filt = all_mutations_filt[min_vaf_present < 0.5,] 
paste('Number of mutations with min VAF < 0.5:', dim(all_mutations_filt)[1]) # should be 58,958

# 3 retain mutations present at max VAF > 0.1 (if lower, likely to be artefactual)
all_mutations_filt = all_mutations_filt[max_vaf > 0.1,] 
paste('Number of mutations with max VAF > 0.1:', dim(all_mutations_filt)[1]) # should be 58,934

##################################################################################################################

# List of mutations present in a given nr of samples 
list_mutations_all = all_mutations[sum_samples==10,mut] %>% unlist() # mut in all samples
list_mutations_private = all_mutations[sum_samples==1,mut] %>% unlist() # mut in 1 sample
list_mutations_nall = all_mut_filt %>% select(mut) %>% unlist() # mut in 2-9 samples

##################################################################################################################
##################################################################################################################
##################################################################################################################
# PHYLOGENY BASED ON NORMAL SAMPLES V1 
# I will be using the dataframe with dodgy mutations filtered out 

# subset only normal samples
cols_tumour= c('tumour1_PD62341', 'tumour2_PD62341', 'tumour1_PD63383', 'tumour2_PD63383')
mut_filt_normal = all_mutations_filt[,.SD, .SDcols = !cols_tumour]
paste('Number of mutations considered to build the phylogeny:', dim(mut_filt_normal)[1]) # 58,934

mut_filt_normal[,max_vaf_normal:= apply(.SD, 1, max), .SDcols = 2:7]
mut_filt_normal[,mean_vaf_normal:= apply(.SD, 1, mean), .SDcols = 2:7]
mut_filt_normal[,min_vaf_normal:= apply(.SD, 1, min), .SDcols = 2:7]
mut_filt_normal[,mean_vaf_present_normal:= apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = 2:7]
mut_filt_normal[,min_vaf_present_normal:= apply(.SD, 1, function(x) min(x[x>0])), .SDcols = 2:7]

# identify mutations present in 2-5 samples (if in 1 or 6 samples, not informative for the phylogeny)
mut_filt_normal_informative = mut_filt_normal[sum_normal %in% c(2:5)]
paste('Number of informative mutations (2-5 samples):', dim(mut_filt_normal[sum_normal %in% c(2:5)])[1])

hist(mut_filt_normal_informative[,sum_samples]) 
# still are present mostly in 9 samples - maybe we should remove those as well?
# I would guess these are more likely to be missing due to an artefact rather than anything else really
# but then what are you going to do if this is present in some cancer samples but not other? did the cancer lose it or what?

# make lists of mutations matching specific branches on the phylogeny  
list2_mut_all_normal_samples = mut_filt_normal[sum_normal==6, mut] %>% unlist() # 24,441
list2_mut_PD62341_normal_samples = mut_filt_normal[sum_PD62341_normal==3 & sum_PD63383_normal==0, mut] %>% unlist() # 353
list2_mut_PD63383_normal_samples = mut_filt_normal[sum_PD62341_normal==0 & sum_PD63383_normal==3, mut] %>% unlist() # 177

# checking on Jbrowse how these look like
# I found a mutation that actually looks like it matches what it says on CaVEMan
# 1:43032865 (also note, this is located with a region with a lot of homo-polymers)
## the problem is, it seems to be absent due to low coverage in other samples (not enough variant reads to call it a mutation in those samples)
## and you cannot filter by this because it is not called with lower nr of reads in the first place anyway
twins_df[mut%in%list2_mut_PD62341_normal_samples & total_depth>50 & vaf<0.2, mut]
# 5:85566998 - what should we do about sth like this? (4 vs 3 - and also very weird genomic seq context)
# 1:5047338 - this is a similar one - could you confidently say this is absent?
# 6:11686704 - another one, also difficult sequence context 
# quite a lot of those are really poorly mapped

##################################################################################################################
##################################################################################################################
##################################################################################################################
# PHYLOGENY BASED ON NORMAL SAMPLES ONLY (V0)

# Subset to only include normal samples
cols_tumour= c('tumour1_PD62341', 'tumour2_PD62341', 'tumour1_PD63383', 'tumour2_PD63383')
all_mut_normal = all_mutations[,.SD, .SDcols = !cols_tumour]

# Add useful summary data (mean / max / min VAF)
# Note that these values change since now tumour samples are not considered
all_mut_normal[,max_vaf_normal:= apply(.SD, 1, max), .SDcols = 2:7]
all_mut_normal[,mean_vaf_normal:= apply(.SD, 1, mean), .SDcols = 2:7]
all_mut_normal[,min_vaf_normal:= apply(.SD, 1, min), .SDcols = 2:7]
all_mut_normal[,mean_vaf_present_normal:= apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = 2:7]
all_mut_normal[,min_vaf_present_normal:= apply(.SD, 1, function(x) min(x[x>0])), .SDcols = 2:7]

paste('Number of unique mutations present in at least 1 normal sample:', dim(all_mut_normal[sum_normal>0])[1]) # 350,779

# FILTERING (dataframe with exact VAF data)
# Remove mutations present in all or only 1 sample
all_mut_normal_filt = all_mut_normal[sum_normal %in% c(2:5)]
paste('Number of informative mutations (present in 2-5 samples):', dim(all_mut_normal_filt)[1])
# should be 32,026

# Plot the distribution of nr of samples with mutation
mut_normal_count = data.table(table(all_mut_normal_filt[,sum_normal]))
ggplot(data=mut_normal_count, aes(x=V1, y=N)) +
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x = 'Number of samples a mutation with the mutation', y = 'Frequency', 
       title = 'Mutation counts distribution in normal samples')

# THE PHYLOGENY WE THINK IS MOST LIKELY IS:
# 1 zygote > 2 twin1 + twin2 > 3 tissues in each twin split 
# NB we expect the same tree in twin1 and twin2 

# Identify mutations private to each twin 
# First, sum nr of samples from each twin with the mutation
all_mut_normal[,sum_PD62341_normal:= rowSums(.SD>0), .SDcols = 2:4]
all_mut_normal[,sum_PD63383_normal:= rowSums(.SD>0), .SDcols = 5:7]

# I want to check what the distribution is of the nr of samples b/n twins 
all_mut_normal[,PD62341_PD63383:= paste(sum_PD62341_normal, sum_PD63383_normal, sep = '_')]
count_normal_dist = data.table(table(all_mut_normal[,PD62341_PD63383]))
levels = c('1_0', '2_0', '3_0', '0_1', '0_2', '0_3', '1_1', '1_2', '2_1', '1_3', '3_1', '2_2', '2_3', '3_2', '0_0', '3_3')
count_normal_dist[,V1 := factor(V1, levels = levels)]
ggplot(data=count_normal_dist[!V1 %in% c('0_0', '3_3')], aes(x=V1, y=N)) +
  geom_bar(stat='identity')+
  theme_bw()+
  labs(x = 'Number of samples from twin1 / twin2 with mutation', y = 'Frequency', 
       title = 'Dist of samples with mutations between twins')

# Create lists of mutations with description in what samples they are present in
mut_list <- list()

for (dist in levels){
  df = all_mut_normal[PD62341_PD63383==dist]
  mut_list[[dist]] <- df[,mut] %>% unlist() %>% unique()
}

# Identify mutations present in all tissues of a specific twin 
mut_PD62341 = all_mut_normal[sum_PD62341_normal==3 & sum_PD63383_normal==0,mut] %>% unlist()
mut_PD63383 = all_mut_normal[sum_PD62341_normal==0 & sum_PD63383_normal==3,mut] %>% unlist()

paste('Number of mutations present in all tissues of PD62341 only:', length(mut_PD62341)) # 379
paste('Number of mutations present in all tissues of PD63383 only:', length(mut_PD63383)) # 201

# We could subset those and check how many of those have a reasonable (< 0.5) VAF
mut_PD62341_vafs = all_mut_normal[mut %in% mut_PD62341, c('mut', 'skin_PD62341', 'pancreas_PD62341', 'spleen_PD62341'), with=FALSE]

# Get average VAF 
# NOTE: mutations on chromosome Y have VAF = 1, this is because there is only one chr Y
# so if all cells have a mutation present on chromosome Y, VAF will be 1
# interestingly, VAF for chromosome X mutation is ~0.5 (> not in all cells?)
mut_PD62341_vafs[, mean_vaf:= apply(.SD, 1, mean), .SDcols = 2:4]

# We can do the same for mutations in PD63383
mut_PD63383_vafs = all_mut_normal[mut %in% mut_PD63383, c('mut', 'skin_PD63383', 'pancreas_PD63383', 'spleen_PD63383'), with=FALSE]
mut_PD63383_vafs[, mean_vaf:= apply(.SD, 1, mean), .SDcols = 2:4]

# PROPOSED PHYLOGENY
# 1: split into twins
# 2: twin1 splits into skin / pancreas-spleen; twin2 splits into skin / pancreas-spleen
# 3: twin1 pancreas-spleen splits into pancreas / spleen, twin2 pancreas-spleen splits into pancreas / spleen 

# we want to find mutations which would 'sit' on specific branches of these
# identify those mutations and write them to lists / check with VAF
# we require VAF > 0.1 to say a mutation is present in the sample 
all_mut_normal[,sum_PD62341_normal_01 := rowSums(.SD > 0.1), .SDcols = c(2:4)]
all_mut_normal[,sum_PD63383_normal_01 := rowSums(.SD > 0.1), .SDcols = c(5:7)]

# mutations present in all normal samples (most likely to be germline)
list_mut_PD62341_PD63383_all_01 = all_mut_normal[sum_PD62341_normal_01==3 & sum_PD63383_normal_01==3, mut] %>% unlist() # 308,801

# mutations present in one twin, but not the other
list_mut_PD62341_only_01 = all_mut_normal[sum_PD62341_normal_01==3 & sum_PD63383_normal_01==0, mut] %>% unlist() # 380
list_mut_PD63383_only_01 = all_mut_normal[sum_PD62341_normal_01==0 & sum_PD63383_normal_01==3, mut] %>% unlist() # 201

list_mut_PD62341_spleen_pancreas = all_mut_normal[spleen_PD62341 > 0 & pancreas_PD62341 > 0 & sum_PD62341_normal==2 & sum_PD63383_normal==0, mut] %>% unlist()
# 472
list_mut_PD63383_spleen_pancreas = all_mut_normal[spleen_PD63383 > 0 & pancreas_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==2, mut] %>% unlist()
# 347 

list_mut_PD62341_spleen = all_mut_normal[spleen_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 1796
list_mut_PD63383_spleen = all_mut_normal[spleen_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 1580

list_mut_PD62341_pancreas = all_mut_normal[pancreas_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 1453
list_mut_PD63383_pancreas = all_mut_normal[pancreas_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 1717

list_mut_PD62341_skin = all_mut_normal[skin_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 1648
list_mut_PD63383_skin = all_mut_normal[skin_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 1680

list_mut_absent_normal = all_mut_normal[sum_normal==0, mut] %>% unlist() # 4656

# which mutations are not congruent with this phylogeny?
# all mutations that are identified in either of the normal samples
list_mut_normal = all_mut_normal[sum_normal!=0, mut] %>% unlist() # 350779

# which mutations are congruent with the phylogeny?
list_mut_congruent = c(list_mut_PD62341_PD63383_all, list_mut_PD62341_only, list_mut_PD63383_only,
                       list_mut_PD62341_spleen_pancreas, list_mut_PD63383_spleen_pancreas,
                       list_mut_PD62341_spleen, list_mut_PD62341_pancreas, list_mut_PD62341_skin,
                       list_mut_PD63383_spleen, list_mut_PD63383_pancreas, list_mut_PD63383_skin) # 320,152
# NB note that this is dominated by mutations present in all samples

# which mutations are incongruent (all mutations - congruent mutations)
list_mut_incongruent = setdiff(list_mut_normal, list_mut_congruent) # 30,627 mutations not congruent with the phylogeny
# do we think all of these would be error / contamination etc?

# we can look at incongruent mutations - are those likely to be germline missed in some samples / artefacts?
hist(all_mut_normal[mut %in% list_mut_congruent, mean_vaf_present], main = 'VAF distribution of congruent mutations', xlab = 'mean VAF')
hist(all_mut_normal[mut %in% list_mut_incongruent, mean_vaf_present], main = 'VAF distribution of incongruent mutations', xlab = 'mean VAF')
# the incongruent mutations have a wider distribution 
# most of these looks quite germline to me? 
# could it be they are missing from a reasonable fraction of samples?
# could it be contamination with cells from a tissue that has a shared origin (then the VAF dist would be different)?

#######################################################################################################################
# LOOK AT MUTATION QUALITY AND INTERESTING MUTATIONS TO LOOK AT!

# mutations present in all samples
list_mut_PD62341_PD63383_all = all_mut_normal[sum_PD62341_normal==3 & sum_PD63383_normal==3, mut] %>% unlist() # 308,879

# mutations present in all samples from one twin but none from the other 
# note that this does not consider if these mutations are also in the tumour or not 
list_mut_PD62341_only = all_mut_normal[sum_PD62341_normal==3 & sum_PD63383_normal==0, mut] %>% unlist() # 379
list_mut_PD63383_only = all_mut_normal[sum_PD62341_normal==0 & sum_PD63383_normal==3, mut] %>% unlist() # 201

# mutations present in some samples from one twin but none from the other
list_mut_PD62341_some = all_mut_normal[sum_PD62341_normal>1 & sum_PD63383_normal==0, mut] %>% unlist() # 379
list_mut_PD63383_some = all_mut_normal[sum_PD62341_normal==0 & sum_PD63383_normal>1, mut] %>% unlist() # 201

# mutations present in single sample from one twin (do we believe these to be real?)
list_mut_PD62341_spleen = all_mut_normal[spleen_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 1796
list_mut_PD63383_spleen = all_mut_normal[spleen_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 1580

list_mut_PD62341_pancreas = all_mut_normal[pancreas_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 1453
list_mut_PD63383_pancreas = all_mut_normal[pancreas_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 1717

list_mut_PD62341_skin = all_mut_normal[skin_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 1648
list_mut_PD63383_skin = all_mut_normal[skin_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 1680

list_mut_absent_normal = all_mut_normal[sum_normal==0, mut] %>% unlist() # 4656
twins_df[mut %in% list_mut_absent_normal & total_depth>50 & chr!='chrY' & vaf > 0.1 & chr=='chr5'] 
# what I am finding is that a lot of those are 1) sequenced to very low depth (like < 30 reads total) or are also seen in JBrowse in the normal samples

# how many mutations we think are absent in the normal samples (given filtering for reasonable coverage)
dim(twins_df[mut %in% list_mut_absent_normal & total_depth>50 & total_depth < 250]) 
# that gets rid of quite a lot of those (from 4.6k to 3.8k)
# note that mutations absent from normal samples are generally called at very low depth (1,739 are at depth < 50)
# so that could suggest that a lot of those are missing because if they are at low VAF they were missed by chance / there are mapping issues 
# some of those (e.g., chr9_87927177_T_A) map absolutely terribly
# X:115726454 looks very concerning (only reads from one direction)

# MUTATIONS THAT ARE INTERESTING TO LOOK AT

# present in ALL normal samples from PD62341 twin only at VAF = 1
# interestingly, 3/4 are on chromosome 1, which we know is lost in the tumour - could that be contamination / artefact?
# need to double check the number of reads on Jbrowse and in the report I have from CaVEMan
# 1:146266994 - quite low mapping quality, so maybe it is not real
# 1:13203287 - low mapping quality, next to a poly-A tract, low depth, reads in 1 orientation
# 1:25287041 - this actually looks okay (good mapping quality, reads in both orientations, but very few reads in the tumour - only see 3 on Jbrowse)
# 9:78356 - low mapping quality, very few reads - probably would not believe it that much 
twins_df[mut %in% list_mut_PD62341_only_01 & vaf==1 & chr!='chrY', mut] %>% unique()

# Mutations present in ALL normal samples from PD63383 twin at VAF = 1
# generally, mutations at very high VAF are present at low depth (few reads mapping well)
twins_df[mut %in% list_mut_PD63383_only & vaf==1, mut] %>% unique()
# 2:27625977 - what does the star mean? is it a deletion (like the base is absent from the read?) but it looks like a convincing mut (also in tumour)
# X:52684628 - this looks rubbish (probably not real)
# 7:63512405 - low mapping quality

# look at mutations that are identified only in one twin and seem to be clonal in this twin (VAF ~ 0.5)
twins_df[mut %in% list_mut_PD63383_only & vaf > 0.5 & total_depth > 60]
# 15:24289784_G_T - why do we get both C and T mutations on Jbrowse?
# 9:15634271, 10:80873839 - it seems again that the mutation is present in other samples (ie normal PD62341) too

# are there any of those that are also found in the tumour?
table(twins_df[mut %in% list_mut_PD63383_only & vaf > 0.5 & total_depth > 60 & sample_cat=='tumour', mut])
# I looked at chr1_244568915_G_A (present in two tumour samples) and that looks okay but then again I don't believe this is restricted to PD63383 samples

# get out positions to see how often the same position is called as mutated in different ways
twins_df[, mut_pos := paste(chr, position, sep='_')]
list_unique_positions = twins_df[, mut_pos] %>% unlist() %>% unique()
paste('Number of unique POSITIONS that are mutated:', length(list_unique_positions)) 
# 355,434 - why then do I have unique mutations 355,435
# but generally seems that CaVEMan called only one mutation per position
# e.g., chr15_24289784 is always called as G to T even through on Jbrowse it shows at G to T or C 

# looked at mutations that are present in all normal samples from the PD63383 but not the other twin, and also tumour from PD62341
# question is how contaminated could the tumour samples be with normal tissues and these are actually normal P62341?
twins_df[mut %in% list_mut_PD63383_only & vaf > 0.5 & total_depth > 60 & sample_name=='tumour1_PD62341']
# chr8_32316023_G_A - looks ok
# chr1_244568915_G_A - looks ok
# chr5_112020282_A_G - this seems to be mutated into either G or C?
# chr4_73729105_A_C - look ok?

# are there any mutations present in all normal PD63383 and ALL tumour samples?
normal_1_tumour = data.table(table(twins_df[mut %in% list_mut_PD63383_only & sample_cat == 'tumour', mut]))
normal_1_tumour = normal_1_tumour[order(N, decreasing=TRUE),]
# then look at mutations in all tumour samples 
# chr10_87249115_A_G - looks a bit dodgy 
# chr15_39176051_C_T - looks real to me!
# 9:33579077 - not mapping correctly
# 8:38769026 - this is quite interesting in that the same site is often mutated to different bases (some to >T)

# Mutations which were not correctly called by CaVEMan?? these all look real
# chr15_57251022_G_T - looks real
# chr2_185542272_A_T
# chr7_20692412_C_T
# chr6_149925861_A_T
# chr7_20692412_C_T
# chr8_8258891_T_C
# I looked on tracks for Jbrowse for all samples (including PD62341 normal) and it seems to be present in all of them, so something must be wrong!

# NB I also looked at some mutations present in fewer tumour samples (2-3) and again some seem to show up on Jbrowse in all samples?
# no apparent issues with mapping quality etc.

# also looked at some mutations with very high coverage - these seem to be mostly mapping incorrectly
# 19:47928138 concerning mapping quality
# 11:94236821 not as bad as the one above in terms of mapping quality, but not great (what is the line on this?)

# do we have mutations present in the spleen of both twins but not anywhere else?
list_spleen_only = all_mutations[spleen_PD63383>0.1 & spleen_PD62341>0.1 & sum_PD63383==1 & sum_PD62341==1, mut] %>% unlist()
twins_df[mut %in% list_spleen_only & total_depth > 50]
# check several of those mutations
# 9:114160652, 13:55886023 - again, this is ~all samples (and also, in all samples look like they are real!)
# 1:75731354 - this is a good one to ask about - Jbrowse is telling me that it is real
# but in a tumour sample it would be called as absent bc only present in 2 reads - what then (esp if your coverage in this sample is lower?)
# ALSO - your coverage can be lower because there is a copy number difference - what then?

# 17:10697578 - NB this one is fun - there is two mutations in phase (A/A) and one anti-phase (T)
# what does this exactly mean? in each sample, you have a pool of reads w one mut and pool w another but they are never together 
# so they must have arisen in different cells - descendants in each sample?

# we could check for each mutation, how wide is the distribution of coverage values?

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# this is not extremely useful 

# can we look at the VAF distribution of congruent mutations
all_mut_normal_congruent = all_mut_normal[mut %in% list_mut_congruent] # 320,152
matrix_normal_congruent = as.matrix(all_mut_normal_congruent[,2:7])
heatmap3(matrix_normal_congruent, useRaster = TRUE, Colv=NA, Rowv=NA)

# how about we do this for mutations at low VAF which could be somatic
# will also remove those present in all samples as these are not informative
all_mut_normal_congruent02 = all_mut_normal_congruent[sum_normal %in% c(2:5) & max_vaf < 0.2]
# convert to binary
all_mut_normal_congruent02_binary = data.table(all_mut_normal_congruent02)
for (j in names(all_mut_normal_congruent02_binary)[2:7])
  {set(all_mut_normal_congruent02_binary, which(is.numeric(all_mut_normal_congruent02_binary[[j]]) & all_mut_normal_congruent02_binary[[j]] > 0), j, 1)}
matrix_normal_congruent02_binary = as.matrix(all_mut_normal_congruent02_binary[,2:7]) # 1,689 mutations
heatmap(matrix_normal_congruent02_binary, col = c('purple', 'green')) # okay but then I guess you gave it mutations that follow this pattern so what did you expect to get really
#######################################################################################################################

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# WHICH TWIN DID THE TUMOUR ARISE IN?
# Two plausible phylogenies:
# 1 tumour arose in twin 1 > carries mutations shared by all twin 1 samples, but not twin 2 samples
# 2 tumour arose in twin 2 > carries mutations shared by all twin 2 samples, but not twin 1 samples
all_mutations[,sum_normal:=rowSums(.SD>0), .SDcols = c(4:6, 9:11)]
all_mutations[,sum_tumour:=rowSums(.SD>0), .SDcols = c(2:3, 7:8)]
all_mutations[,sum_samples:=sum_normal+sum_tumour]

# exclude mutations present in all or 1 sample (not informative)
all_mutations_filt = all_mutations[sum_samples %in% c(2:9)]
paste('Number of informative mutations:', dim(all_mutations_filt)[1]) # 61,806

all_mutations_filt[,sum_PD62341_normal:= rowSums(.SD>0), .SDcols = 4:6]
all_mutations_filt[,sum_PD63383_normal:= rowSums(.SD>0), .SDcols = 9:11]
all_mutations_filt[,sum_PD62341_tumour:= rowSums(.SD>0), .SDcols = 2:3]
all_mutations_filt[,sum_PD63383_tumour:= rowSums(.SD>0), .SDcols = 7:8]

# more informative to know the mean / min VAF across samples the var was actually detected 
all_mutations_filt[,max_vaf:= apply(.SD, 1, function(x) max(x[x>0])), .SDcols = 2:11]
all_mutations_filt[,min_vaf:= apply(.SD, 1, function(x) min(x[x>0])), .SDcols = 2:11]
all_mutations_filt[,mean_vaf:= apply(.SD, 1, function(x) mean(x[x>0])), .SDcols = 2:11]

### HYPOTHESIS TESTING - CURRENT 
# we know that the cancer is SARCOMA
# SARCOMA mainly develops from mesodermal tissue (possibly ectoderm)
# skin = ectoderm, spleen = mesoderm, pancreas = endoderm

# We expect that the mutation will be shared by the tumour (from both twins) and normal samples from one, but not the other twin
# I guess we'd look for mutations that gave rise to the tumour - so they should be really high in the tumour samples
# on the other hand, maybe better not make assumptions about which normal cell the cancer came from
# note that it would make sense if it also shared mutations that were present in all cells 

mut_h1_PD62341 = all_mutations_filt[sum_PD63383_tumour==2 & sum_PD62341_tumour==2 & sum_PD63383_normal==0 & sum_PD62341_normal>0] # 189
mut_h1_PD63383 = all_mutations_filt[sum_PD63383_tumour==2 & sum_PD62341_tumour==2 & sum_PD63383_normal>0 & sum_PD62341_normal==0] # 115 
# okay so we have mutations that kind of match both phylogenies so what do we do now 

# with what samples do these tumours mostly share mutations?
dim(mut_h1_PD62341[skin_PD62341>0]) # 134
dim(mut_h1_PD62341[pancreas_PD62341>0]) # 119
dim(mut_h1_PD62341[spleen_PD62341>0]) # 132

dim(mut_h1_PD63383[skin_PD63383>0]) # 81
dim(mut_h1_PD63383[pancreas_PD63383>0]) # 51
dim(mut_h1_PD63383[spleen_PD63383>0]) # 52

dim(mut_h1_PD62341[skin_PD62341>0 & pancreas_PD62341>0 & spleen_PD62341>0]) # 60
dim(mut_h1_PD62341[skin_PD62341==0 & pancreas_PD62341>0 & spleen_PD62341>0]) # 32
dim(mut_h1_PD62341[skin_PD62341==0 & pancreas_PD62341==0 & spleen_PD62341>0]) # 13

# okay so across those there are 50-130s of those
# but realistically, most of those could be germline (just not picked up) or artefacts 
# germline will be present at high VAF and artefacts at rather low VAF
# we need to find mutations that look reliable (not germline, not artefacts) that match one but not the other phylogeny

### HYPOTHESIS TESTING 1
# HYPOTHESIS 1
# If the tumour arose in twin PD62431 before skin-spleen-pancreas split, 
# we can expect mutations present in:
# all samples from normal tissues of the twin (PD62431)
# all samples of the tumour of this twin (PD62431)
# no normal samples, but all tumour samples of the other twin (PD63383)
mut_ht1e = all_mutations_filt[sum_PD63383_normal==0 & sum_PD63383_tumour==2 & sum_PD62341_normal==3 & sum_PD62341_tumour==2]
mut_ht1e[,'hypothesis' := 'PD62341']
paste('Number of mutations supporting tumour arising from twin PD62431:', dim(mut_ht1e)[1])
list_mut_ht1e = mut_ht1e[,mut] %>% unlist()
# hypothesis = twin 1, early (before skin-pancreas-spleen split)

# HYPOTHESIS 2
# If the tumour arose in twin PD62431 before skin-spleen-pancreas split, 
# we can expect mutations present in:
# all samples from normal tissues of the twin (PD63383)
# all samples of the tumour of this twin (PD63383)
# no normal samples, but all tumour samples of the other twin (PD62431)
mut_ht2e = all_mutations_filt[sum_PD63383_normal==3 & sum_PD63383_tumour==2 & sum_PD62341_normal==0 & sum_PD62341_tumour==2]
mut_ht2e[,'hypothesis':= 'PD63383']
paste('Number of mutations supporting tumour arising from twin PD63383:', dim(mut_ht2e)[1])
list_mut_ht2e = mut_ht2e[,mut] %>% unlist()
# hypothesis = twin 2, early (before skin-pancreas-spleen split)
mut_hte = data.table(rbind(mut_ht1e, mut_ht2e))

mut_info = twins_df[, c("mut", "info"), with=FALSE]
mut_he_info = merge(mut_info, mut_hte)
mut_he_info = mut_he_info[,c('mut', 'info', 'mean_vaf', 'max_vaf', 'min_vaf', 'hypothesis')]

# Identify the genes affected by mutations and how they are affected
mut_he_info[,mut_effects:=tstrsplit(info, 'VD=', keep=2)]
mut_he_info[,gene_affected:=tstrsplit(mut_effects,"\\|", keep=1)] # this is how you split '|'
mut_he_info[,if_coding:=tstrsplit(mut_effects, '\\|', keep=6)]
mut_he_info[,if_coding:=tstrsplit(if_coding, ':', keep=1)]
mut_he_info[,if_coding:=as.factor(if_coding)]
mut_he_info[,mut_consequence:=tstrsplit(mut_effects, 'VC=', keep=2)]
mut_he_info[,mut_consequence:=tstrsplit(mut_consequence, ';VW', keep=1)]
mut_he_info[,mut_consequence:=as.factor(mut_consequence)]

# get rid of the 'info' column
mut_he_info = mut_he_info[,c("info", "mut_effects"):= NULL]
mut_he_info = unique(mut_he_info)

# OKAY, NOW WE TEST THAT THE TUMOUR CELLS AROSE AFTER SKIN / PANCREAS-SPLEEN LINEAGES SPLIT
# HYPOTHESIS 1
# If the tumour arose in twin PD62431 AFTER skin-spleen-pancreas split, 
# we can expect mutations present in:
# pancreas and spleen from normal tissues of the twin (PD62431)
# all samples of the tumour of this twin (PD62431)
# no normal samples, but all tumour samples of the other twin (PD63383)
mut_ht1l = all_mutations_filt[sum_PD63383_normal==0 & sum_PD63383_tumour==2 & sum_PD62341_normal==2 & sum_PD62341_tumour==2 & spleen_PD62341 > 0 & pancreas_PD62341 > 0]
mut_ht1l[,'hypothesis' := 'PD62341']
paste('Number of mutations supporting tumour arising from twin PD62431:', dim(mut_ht1l)[1])
list_mut_ht1l = mut_ht1l[,mut] %>% unlist()
# hypothesis = twin 1, early (before skin-pancreas-spleen split)

# HYPOTHESIS 2
# If the tumour arose in twin PD62431 AFTER skin-spleen-pancreas split, 
# we can expect mutations present in:
# pancreas and spleen of the twin (PD63383)
# all samples of the tumour of this twin (PD63383)
# no normal samples, but all tumour samples of the other twin (PD62431)
mut_ht2l = all_mutations_filt[sum_PD63383_normal==2 & sum_PD63383_tumour==2 & sum_PD62341_normal==0 & sum_PD62341_tumour==2 & spleen_PD63383 > 0 & pancreas_PD63383 > 0]
mut_ht2l[,'hypothesis':= 'PD63383']
paste('Number of mutations supporting tumour arising from twin PD63383:', dim(mut_ht2l)[1])
list_mut_ht2l = mut_ht2l[,mut] %>% unlist()
# hypothesis = twin 2, early (before skin-pancreas-spleen split)
mut_htl = data.table(rbind(mut_ht1l, mut_ht2l))

mut_info = twins_df[, c("mut", "info"), with=FALSE]
mut_hl_info = merge(mut_info, mut_htl)
mut_hl_info = mut_hl_info[,c('mut', 'info', 'mean_vaf', 'max_vaf', 'min_vaf', 'hypothesis')]

# Identify the genes affected by mutations and how they are affected
mut_hl_info[,mut_effects:=tstrsplit(info, 'VD=', keep=2)]
mut_hl_info[,gene_affected:=tstrsplit(mut_effects,"\\|", keep=1)] # this is how you split '|'
mut_hl_info[,if_coding:=tstrsplit(mut_effects, '\\|', keep=6)]
mut_hl_info[,if_coding:=tstrsplit(if_coding, ':', keep=1)]
mut_hl_info[,if_coding:=as.factor(if_coding)]
mut_hl_info[,mut_consequence:=tstrsplit(mut_effects, 'VC=', keep=2)]
mut_hl_info[,mut_consequence:=tstrsplit(mut_consequence, ';VW', keep=1)]
mut_hl_info[,mut_consequence:=as.factor(mut_consequence)]

# get rid of the 'info' column
mut_hl_info = mut_hl_info[,c("info", "mut_effects"):= NULL]
mut_hl_info = unique(mut_hl_info)

# select tumour-specific mutations
tumour_mut = all_mutations[sum_tumour==4 & sum_normal==0, c('mut', 'tumour1_PD62341', 'tumour2_PD62341',
                                                            'tumour1_PD63383', 'tumour2_PD63383'), with=FALSE]
heatmap(as.matrix(tumour_mut[, c(2:5), with=FALSE]))

#############################################################################################
#############################################################################################
#############################################################################################

# DRIVERS  

# How do we expect a driver mutation to behave?

# Based on pattern of mutation distribution
# Present in tumour, but not normal samples
# Present in a reasonable fraction of tumour cells (arose early)
# Present in tumour samples from both twins (because this is the same tumour really)

# Based on biology:
# Present in the coding region / regulatory region
# Affects the protein sequence / expression
# Affects function of a gene that plays a role in embryogenesis

# First, filter for mutations which are only present in tumour samples 
# In this way, we can reduce the size of the dt so this actually runs

list_mut_tumour_only_all = all_mutations[sum_tumour==4 & sum_normal==0, mut] %>% unlist()
list_mut_tumour_only = all_mutations[sum_tumour>=1 & sum_normal==0, mut] %>% unlist()
twins_df_mut_tumour_all = twins_df[mut %in% list_mut_tumour_only_all] 
twins_df_mut_tumour = twins_df[mut %in% list_mut_tumour_only] 
paste('Number of mutations observed in all tumour, but no normal samples:', dim(twins_df_mut_tumour_all)[1]) # 128 in all tumour samples
paste('Number of mutations observed in some tumour, but no normal samples:', dim(twins_df_mut_tumour)[1]) # 5756 in at least 1 tumour sample

# what is the VAF distribution of those mutations?
mut_status_tumour = all_mutations[mut %in% list_mut_tumour_only_all]
mut_status_tumour[, mean_vaf := apply(.SD, 1, mean), .SDcols = c(2:3, 7:8)]
hist(mut_status_tumour[, mean_vaf], xlab = 'Mean VAF', main = 'Mutations in all tumour samples only')
# Interestingly, also present at a really high VAF? Still germline?

# Extract information from the column 'info' from VD=
twins_df_mut_tumour[,mut_effects:=tstrsplit(info, 'VD=', keep=2)]
twins_df_mut_tumour[,gene_affected:=tstrsplit(mut_effects,"\\|", keep=1)] # this is how you split '|'
twins_df_mut_tumour[,if_coding:=tstrsplit(mut_effects, '\\|', keep=6)]
twins_df_mut_tumour[,if_coding:=tstrsplit(if_coding, ':', keep=1)]
twins_df_mut_tumour[,if_coding:=as.factor(if_coding)]
twins_df_mut_tumour[,mut_consequence:=tstrsplit(mut_effects, 'VC=', keep=2)]
twins_df_mut_tumour[,mut_consequence:=tstrsplit(mut_consequence, ';VW', keep=1)]
twins_df_mut_tumour[,mut_consequence:=as.factor(mut_consequence)]

# maybe filter for protein coding first
twins_df_mut_tumour_coding = twins_df_mut_tumour[if_coding=='protein_coding']
paste('Number of unique coding mutations:', length(twins_df_mut_tumour_coding[,mut] %>% unlist() %>% unique())) # 2055
paste('Number of affected genes:', length(twins_df_mut_tumour_coding[,gene_affected] %>% unlist() %>% unique())) # 1054

# we can check the number of substitution mutations
length(twins_df_mut_tumour_coding[mut_consequence=='missense',mut] %>% unlist() %>% unique()) 

# which of these are mutated in all tumour samples (more likely to be the early drivers)
paste('Number of unique coding mutations found in all tumour samples:', 
      length(twins_df_mut_tumour_coding[mut %in% list_mut_tumour_only_all,mut] %>% unlist() %>% unique())) # 52
paste('Number of genes affected in all tumour samples:', 
      length(twins_df_mut_tumour_coding[mut %in% list_mut_tumour_only_all,gene_affected] %>% unlist() %>% unique())) # 52

# Information from here: https://www.proteinatlas.org/ 

# "AC242842.3" - novel protein, identical to neuroblastoma breakpoint family, member 19 NBPF19
## well given the name (neuroblastoma) must be somehow associated with cancer? 
## interestingly, one sample from one twin also has a missense mutation
# "POLR3B" - RNA polymerase III subunit B     
# "GALNT9" - Polypeptide N-acylgalactosaminyltransferase      
# "LGMN" - Legumain (enzyme, metabolic protein)     
# "BBS4" - Bardet-Biedl Syndrome 4, disease related genes  
## interesting find - rare disease, but doesn't seem to be implicated in cancer (maybe renal cancer?)
# "DAPK3" - Death-associated protein kinase 3     
# "HECW2" - Enzyme, known to be disease-related      
# "BID" - BH3 interacting domain death agonist      
# "GAL3ST1" - Galactose-3-O-sulfotransferase  
# "PHF3" - PHD finger protein 3      
# "HIPK2" - homeodomain interacting protein kinase 2    
# "SHROOM2" - Shroom family member 2 (cytoplasmic, localizes to cell junctions) 
# "TENM1" - teneurin transmembrane protein, transporter   

# Note that none of those mutations is a substitution - mostly intronic
# However, this can affect expression or splicing so not impossible to have an effect 

#############################################################################################
#############################################################################################
#############################################################################################

# IS THERE A MUTATION INDICATIVE OF CANCER PREDISPOSITION?
# We expect such a mutation to be present in all normal samples
list_mut_normal_all = all_mutations[sum_normal==6, mut] %>% unlist() # note that this mutation can be present in the tumour (we expect it to be present)
twins_df_mut_normal_all = twins_df[mut %in% list_mut_normal_all] # 3,026,722 of those 

twins_df_mut_normal_all[,mut_effects:=tstrsplit(info, 'VD=', keep=2)]
twins_df_mut_normal_all[,gene_affected:=tstrsplit(mut_effects,"\\|", keep=1)] # this is how you split '|'
twins_df_mut_normal_all[,if_coding:=tstrsplit(mut_effects, '\\|', keep=6)]
twins_df_mut_normal_all[,if_coding:=tstrsplit(if_coding, ':', keep=1)]
twins_df_mut_normal_all[,if_coding:=as.factor(if_coding)]
twins_df_mut_normal_all[,mut_consequence:=tstrsplit(mut_effects, 'VC=', keep=2)]
twins_df_mut_normal_all[,mut_consequence:=tstrsplit(mut_consequence, ';VW', keep=1)]
twins_df_mut_normal_all[,mut_consequence:=as.factor(mut_consequence)]

# we can now search for protein coding missense and nonsense variants
sort(twins_df_mut_normal_all[if_coding=='protein_coding' & mut_consequence %in% c('missense', 'nonsense'), gene_affected] %>% unlist() %>% unique())
# there are 1,078 genes that are mutated in these samples with a non-sense / mis-sense outcome
# some interesting hits - TP53 target 5, IGF2R, SMARCA1 (NOT SMARCB1 though), SLC family genes, AKAP,
# ALG9, ATP subunits, BCL2L, CDK genes, DDX, FGF / FGFR, HLA-B, HOXA7, histones, ILs, KMT2A (supposedly rearranged in sarcoma)
# MYOs, NOTCH2, PARP, PIK3CD, SEMA3F, SLC family, TASOR, TOP2B, XRCC1, ZNFs (~30 of ZNFs)   
# it would be good to ask if there are specific genes used clinically to identify sarcoma predisposition

#############################################################################################
#############################################################################################
#############################################################################################

# HAVING FUN WITH BINOMIAL DISTRIBUTION: FILTERING OUT GERMLINE MUTATIONS

# how this is done in THHC 2019 (Science paper, supplementary methods section):
# We fitted a binomial distribution to the combined read counts of all normal samples from one
# patient per SNV site, with the total depth as the number of trials, and the total number of reads
# supporting the variant as number of successes. Germline and somatic variants were differentiated
# based on a one-sided exact binomial test. For this test, the null hypothesis is that the number of
# reads supporting the variants across copy number normal samples is drawn from a binomial
# distribution with p=0.5 (p=0.95 for copy number equal to one), and the alternative hypothesis
# drawn from a distribution with p<0.5 (or p<0.95). Resulting p-values were corrected for multiple
# testing with the Benjamini-Hochberg method (19) and a cut-off was set at q < 10-5 to minimize
# false positives as on average, roughly 40,000 variants were subjected to this statistical test.
# Variants for which the null hypothesis could be rejected were classified as somatic, otherwise as
# germline. When parental genomes were sequenced, de novo germline mutations were taken to be
# those variants classified as germline in the child, but absent in both parents.

# question to this: which samples should I use given twins are related and share some mutations?
# does it make sense to look at normal samples from either of the twins?
# is BH method markedly better than Bonferroni correction that I am using?

# TWIN 1 BINOMIAL 
# take normal samples from one of the twins
normal_PD62341 = all_mutations[,c('mut', 'skin_PD62341', 'pancreas_PD62341', 'spleen_PD62341'), with=FALSE]

# we need the total depth and the nr of reads supporting for each sample

# extract information on total depth and number of reads with the variant
PD62341_depths = twins_df[sample_name %in% c('skin_PD62341', 'pancreas_PD62341', 'spleen_PD62341'),
                        c("mut",'sample_name', 'total_depth', 'variant_depth'), with=FALSE]

# sum total depth and read depth for each mutation across all three samples 
PD62341_depths_agg = aggregate(PD62341_depths[,c("total_depth", "variant_depth")], by=list(Category=PD62341_depths[,mut]), FUN=sum)

# now, determine if there is evidence that true p (VAF) is < 0.5
bin_test <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=
                                            c("less"), conf.level = 0.95)$p.value}
PD62341_depths_agg$p_val = mapply(bin_test, PD62341_depths_agg$variant_depth, PD62341_depths_agg$total_depth)

# correct for multiple testing (Bonferroni correction for nr of tests)
PD62341_depths_agg = data.table(PD62341_depths_agg)
PD62341_depths_agg[, significance := as.factor(fcase( # add what sample that was (based on data from Henry)
  p_val < 0.05/dim(PD62341_depths_bin_agg)[1], 'significant',
  p_val >= 0.05/dim(PD62341_depths_bin_agg)[1], 'not_significant'
))]

list_putative_somatic_PD62341 = PD62341_depths_agg[significance=='significant', Category] %>% unlist() 
paste('Number of putative somatic mutations in PD62341:', length(list_putative_somatic_PD62341)) # 3368

# TWIN 2 BINOMIAL

# extract information on total depth and number of reads with the variant
PD63383_depths = twins_df[sample_name %in% c('skin_PD63383', 'pancreas_PD63383', 'spleen_PD63383'),
                          c("mut", 'sample_name', 'total_depth', 'variant_depth'), with=FALSE]

# sum total depth and read depth for each mutation across all three samples 
PD63383_depths_agg = aggregate(PD63383_depths[,c("total_depth", "variant_depth")], by=list(Category=PD63383_depths[,mut]), FUN=sum)

# now, determine if there is evidence that true p (VAF) is < 0.5
bin_test <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=
                                                  c("less"), conf.level = 0.95)$p.value}
PD63383_depths_agg$p_val = mapply(bin_test, PD63383_depths_agg$variant_depth, PD63383_depths_agg$total_depth)

# correct for multiple testing (Bonferroni correction for nr of tests)
PD63383_depths_agg = data.table(PD63383_depths_agg)
PD63383_depths_agg[, significance := as.factor(fcase( # add what sample that was (based on data from Henry)
  p_val < 0.05/dim(PD63383_depths_agg)[1], 'significant',
  p_val >= 0.05/dim(PD63383_depths_agg)[1], 'not_significant'
))]

list_putative_somatic_PD63383 = PD63383_depths_agg[significance=='significant', Category] %>% unlist() 
paste('Number of putative somatic mutations in PD63383:', length(list_putative_somatic_PD63383))

# determine the intersection between somatic identified in PD62341 and PD63383
paste('Number of shared somatic mutations:', length(Reduce(intersect, list(list_significant_somatic_PD62341, list_significant_somatic_PD63383))))
# I would probably be fine saying that it's fine if the mutation is deemed somatic in either of the twins
# on the other hand, given the nr of somatic mutations expected is lower than the intersection, maybe you'd be fine with the intersection

# okay so now we have a set of mutations that are likely not germline
# we can filter the set of all mutations to only include those and see what phylogenies we are getting out of it in this way
list_somatic_all = c(list_significant_somatic_PD62341, list_significant_somatic_PD63383) %>% unique()
list_somatic_intersect = Reduce(intersect, list(list_significant_somatic_PD62341, list_significant_somatic_PD63383))
paste('Number of putative somatic mutations (in either twin):', length(list_somatic_all)) # 3368 + 3020 - 1838 = 4550 
paste('Number of putative somatic mutations (in both twins):', length(list_somatic_intersect)) # 1838  

all_mutations_somatic = all_mutations[mut %in% list_somatic_all] # 4550
all_mutations_somatic_filt = all_mutations_somatic[!sum_samples %in% c(1, 10)] # 2303
# mutations present in only 1 or all samples are not informative

############################################################################################################
# phylogeny of normal samples from both twins 

list_som_mut_PD62341_PD63383_all = all_mutations_somatic[sum_PD62341_normal==3 & sum_PD63383_normal==3, mut] %>% unlist() # 1,625

list_som_mut_PD62341_only = all_mutations_somatic[sum_PD62341_normal==3 & sum_PD63383_normal==0, mut] %>% unlist() # 25
list_som_mut_PD63383_only = all_mutations_somatic[sum_PD62341_normal==0 & sum_PD63383_normal==3, mut] %>% unlist() # 23

list_som_mut_PD62341_spleen_pancreas = all_mutations_somatic[spleen_PD62341 > 0 & pancreas_PD62341 > 0 & sum_PD62341_normal==2 & sum_PD63383_normal==0, mut] %>% unlist()
# 47
list_som_mut_PD63383_spleen_pancreas = all_mutations_somatic[spleen_PD63383 > 0 & pancreas_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==2, mut] %>% unlist()
# 30 

list_som_mut_PD62341_spleen = all_mutations_somatic[spleen_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 149
list_som_mut_PD63383_spleen = all_mutations_somatic[spleen_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 194

list_som_mut_PD62341_pancreas = all_mutations_somatic[pancreas_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 251
list_som_mut_PD63383_pancreas = all_mutations_somatic[pancreas_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 163

list_som_mut_PD62341_skin = all_mutations_somatic[skin_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 273
list_som_mut_PD63383_skin = all_mutations_somatic[skin_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 231

# which mutations are not congruent with this phylogeny?
# all mutations that are identified in either of the normal samples
list_som_mut_normal = all_mutations_somatic[sum_normal!=0, mut] %>% unlist() # 4550

# which mutations are congruent with the phylogeny?
list_som_mut_congruent = c(list_som_mut_PD62341_PD63383_all, list_som_mut_PD62341_only, list_som_mut_PD63383_only,
                       list_som_mut_PD62341_spleen_pancreas, list_som_mut_PD63383_spleen_pancreas,
                       list_som_mut_PD62341_spleen, list_som_mut_PD62341_pancreas, list_som_mut_PD62341_skin,
                       list_som_mut_PD63383_spleen, list_som_mut_PD63383_pancreas, list_som_mut_PD63383_skin) # 3011

# which mutations are incongruent (all mutations - congruent mutations)
list_som_mut_incongruent = setdiff(list_som_mut_normal, list_som_mut_congruent) # 1539
# that's quite a lot! as a fraction of total, that's a third of all we think are somatic 

# we can look at incongruent mutations - are those likely to be germline missed in some samples / artefacts?
# clearly different distributions
hist(all_mutations_somatic[mut %in% list_som_mut_congruent, mean_vaf_present], main = 'VAF distribution of congruent mutations', xlab = 'mean VAF')
hist(all_mutations_somatic[mut %in% list_som_mut_incongruent, mean_vaf_present], main = 'VAF distribution of incongruent mutations', xlab = 'mean VAF')

#################################################################################################################
# try doing the same thing for intersect (mutations called as somatic in both twins?)
all_mutations[,sum_PD62341_normal := rowSums(.SD>0), .SDcols = c(4:6)]
all_mutations[,sum_PD62341_tumour := rowSums(.SD>0), .SDcols = c(2:3)]
all_mutations[,sum_PD63383_normal := rowSums(.SD>0), .SDcols = c(9:11)]
all_mutations[,sum_PD63383_tumour := rowSums(.SD>0), .SDcols = c(7:8)]

all_mutations_somatic2 = all_mutations[mut %in% list_somatic_intersect] # 4550
all_mutations_somatic2_filt = all_mutations_somatic2[!sum_samples %in% c(1, 10)] # 930

# 60% are present in all samples 
list_som2_mut_PD62341_PD63383_all = all_mutations_somatic2[sum_PD62341_normal==3 & sum_PD63383_normal==3, mut] %>% unlist() # 1,101

# you have none of those - but then maybe that makes sense as if these were present in all 3 they would be classified as germline in the sample?
list_som2_mut_PD62341_only = all_mutations_somatic2[sum_PD62341_normal==3 & sum_PD63383_normal==0, mut] %>% unlist() # 0
list_som2_mut_PD63383_only = all_mutations_somatic2[sum_PD62341_normal==0 & sum_PD63383_normal==3, mut] %>% unlist() # 0

# none of those either! 
list_som2_mut_PD62341_spleen_pancreas = all_mutations_somatic2[spleen_PD62341 > 0 & pancreas_PD62341 > 0 & sum_PD62341_normal==2 & sum_PD63383_normal==0, mut] %>% unlist()
list_som2_mut_PD63383_spleen_pancreas = all_mutations_somatic2[spleen_PD63383 > 0 & pancreas_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==2, mut] %>% unlist()

# no cases where a mutation is present in just one sample 
list_som2_mut_PD62341_spleen = all_mutations_somatic2[spleen_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 0
list_som2_mut_PD63383_spleen = all_mutations_somatic2[spleen_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 0
list_som2_mut_PD62341_pancreas = all_mutations_somatic2[pancreas_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 0
list_som2_mut_PD63383_pancreas = all_mutations_somatic2[pancreas_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 0
list_som2_mut_PD62341_skin = all_mutations_somatic2[skin_PD62341 > 0 & sum_PD62341_normal==1 & sum_PD63383_normal==0, mut] %>% unlist() # 0
list_som2_mut_PD63383_skin = all_mutations_somatic2[skin_PD63383 > 0 & sum_PD62341_normal==0 & sum_PD63383_normal==1, mut] %>% unlist() # 0

# which mutations are congruent with the phylogeny?
list_som2_mut_congruent = c(list_som2_mut_PD62341_PD63383_all, list_som2_mut_PD62341_only, list_som2_mut_PD63383_only,
                           list_som2_mut_PD62341_spleen_pancreas, list_som2_mut_PD63383_spleen_pancreas,
                           list_som2_mut_PD62341_spleen, list_som2_mut_PD62341_pancreas, list_som2_mut_PD62341_skin,
                           list_som2_mut_PD63383_spleen, list_som2_mut_PD63383_pancreas, list_som2_mut_PD63383_skin) # 1101

# which mutations are incongruent (all mutations - congruent mutations)
list_som2_mut_incongruent = setdiff(list_som2_mut_normal, list_som2_mut_congruent) # 737 
# this is 40% - essentially, either present in all samples or in a pattern that doesn't make much sense 

# clearly different distributions - maybe the incongruent ones are errors?
hist(all_mutations_somatic2[mut %in% list_som2_mut_congruent, mean_vaf_present], main = 'VAF distribution of congruent mutations', xlab = 'mean VAF')
hist(all_mutations_somatic2[mut %in% list_som2_mut_incongruent, mean_vaf_present], main = 'VAF distribution of incongruent mutations', xlab = 'mean VAF')

# could the incongruent mutations be errors? 
# we required at least 4 reads to call them, so the chances these are sequencing artefacts is (to me) quite low - ?
# I don't think the distribution of read depths is massively different to all the other mutations
hist(twins_df[mut %in% list_som2_mut_incongruent, total_depth])

############################################################################################################
# using this to do the tumour phylogeny 

# we can again check which phylogeny is better matched looking at mutations
mut_h1_PD62341_somatic = all_mutations_somatic_filt[sum_tumour==4 & sum_PD63383_normal==0 & sum_PD62341_normal>0] # 13
mut_h1_PD63383_somatic = all_mutations_somatic_filt[sum_tumour==4 & sum_PD63383_normal>0 & sum_PD62341_normal==0] # 24 

# okay so we have mutations supporting both so this is not terribly helpful
# I guess all 24 / all 13 from a given phylogeny would have to be artefactual?
# what I guess is a bit odd is that in both cases, these shared mutations are with the skin
# some sarcomas are believed to arise from ectodermal tissue - maybe this one too?
# MAYBE the stuff from the skin is contamination (some skin sampled together with the tumour)?? 
# especially that it is often at lower VAF? - would that be more likely?

mut_h1_PD62341_somatic3 = all_mutations_somatic_filt[sum_tumour>=3 & sum_PD63383_normal==0 & sum_PD62341_normal>0] # 38
mut_h1_PD63383_somatic3 = all_mutations_somatic_filt[sum_tumour>=3 & sum_PD63383_normal>0 & sum_PD62341_normal==0] # 46 

############################################################################################################


