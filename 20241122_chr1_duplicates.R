###################################################################################################################################
# SCRIPT 6

# Script to try to identify the reason for inflated VAF in chr1:38827952;C>A 
# We think this is an early embryonic mutation that came from a cell lineage that contributed unevenly to both twins 
# 2024-11-23
# Barbara Walkowiak bw18

# INPUT: 
# 1 sam file with reads spanning chr1_38827952_C_A (.sam file for each sample)

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
library(Biostrings)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

###################################################################################################################################
# INPUT FILES 

# read sam files covering this mutation for each sample 
# each sam file has all reads spanning the position in a given sample, and associated read parameters  
sam = data.table()
for (sample in samples_names){
  sam_dt = data.table(read.table(paste0('Data/chr1_388_reads_', sample, '_short.txt'), sep='\t', fill = TRUE, header=FALSE))
  sam_dt = sam_dt[, c(1:10)]
  sam_dt[, sample := sample]
  sam = rbind(sam, sam_dt)
}

###################################################################################################################################
# A possible explanation for the inflation of VAF is preferential PCR duplication of reads with the mutant allele 
# Alternatively, wt reads are mapping less well than mutant reads or are excluded for other reasons 

# Are mutant reads over-represented in reads flagged as duplicates?
sam[, dist_to_mut := 38827952 - V4]
sam[, mut_site := substr(V10, dist_to_mut+1, dist_to_mut+1)]

table(sam[, mut_site]) # 0.565 is contributed by wt reads, oh that's not great right 
# mutation is C TO A, therefore C is ref and A is mutant 
# calculate VAF of A across all reads
paste('VAF A (all reads):', round(339/(339+564+2+4), 4)) # 0.3729

table(sam[V2 >= 1024, mut_site]) # 0.565 is contributed by wt reads, oh that's not great right 
paste('VAF A (duplicate reads):', round(26 / (26+37+1), 4)) # 0.4062

table(sam[V2 >= 512, mut_site]) # 0.565 is contributed by wt reads, oh that's not great right 
paste('VAF A (duplicate reads / QC failed):', round(35 / (35+49+1+2), 4)) # 0.4023

table(sam[V2 < 512, mut_site])
paste('VAF A (okay reads):', round(304 / (304+515+1+2), 4)) # 0.3698

###################################################################################################################################
# Dataframe with final values from the pileup 
# get a dataframe with MTR, DEP, VAF values for each one (agg + forward + reverse separately)
chr1_388_dt = data.table()
for (sample in samples_names){
  sample_vaf = paste0(sample, '_VAF')
  sample_mtr = paste0(sample, '_MTR')
  sample_dep = paste0(sample, '_DEP')
  sample_forward_vaf = paste0('forward_vaf_', sample)
  sample_reverse_vaf = paste0('reverse_vaf_', sample)
  sample_forward_mut = paste0('forward_mut_', sample)
  sample_reverse_mut = paste0('reverse_mut_', sample)
  sample_forward_wt = paste0('forward_wt_', sample)
  sample_reverse_wt = paste0('reverse_wt_', sample)
  sample_dt = twins_filtered_dt[mut_ID == 'chr1_38827952_C_A', c(sample_mtr, sample_vaf, sample_dep, 
                                                                 sample_forward_vaf, sample_reverse_vaf,
                                                                 sample_forward_mut, sample_reverse_mut,
                                                                 sample_forward_wt, sample_reverse_wt), with=FALSE]
  setnames(sample_dt, c(sample_mtr, sample_vaf, sample_dep, sample_forward_vaf, sample_reverse_vaf,
                        sample_forward_mut, sample_reverse_mut,sample_forward_wt, sample_reverse_wt),
           c('MTR', 'VAF', 'DEP', 'forward_vaf', 'reverse_vaf', 'forward_mut', 'reverse_mut', 'forward_wt', 'reverse_wt'))
  sample_dt[, forward_dep := forward_mut + forward_wt]
  sample_dt[, reverse_dep := reverse_mut + reverse_wt]
  sample_dt[, vaf_ratio := log2(forward_vaf / reverse_vaf)]
  sample_dt[, dep_ratio := log2(forward_dep / reverse_dep)]
  sample_dt[, sample := sample]
  chr1_388_dt = rbind(chr1_388_dt, sample_dt)
}

# check the number of reads mapped in each sample
table(sam[, sample])

sam_all = data.table(table(sam[, sample])) # duplicates are flagged with 1024 
setnames(sam_all, c('V1', 'N'), c('sample', 'all_reads'))
sam_dup = data.table(table(factor(sam[V2 >= 1024, sample], levels = samples_names))) # duplicates are flagged with 1024 
setnames(sam_dup, c('V1', 'N'), c('sample', 'dup_reads'))

chr1_388_dt = merge(chr1_388_dt, sam_all, by = 'sample')
chr1_388_dt = merge(chr1_388_dt, sam_dup, by = 'sample')
chr1_388_dt[, diff_dep := all_reads - DEP]

colSums(chr1_388_dt[, c('MTR', 'DEP', 'forward_mut', 'reverse_mut', 'forward_dep', 'reverse_dep'), with=FALSE])
# MTR = 275
# DEP = 751 
# VAF across all samples: 0.3662
# Forward MTR = 151
# Forward DEP = 399
# VAF across forward reads = 151/399 = 0.3784
# Reverse MTR = 124
# Reverse DEP = 352
# VAF across reverse reads = 124/352 = 0.3523

###################################################################################################################################
# Plots to check relationship between the nr of reads marked as duplicates in a sample and relevant parameters 

plot(chr1_388_dt[, dup_reads], chr1_388_dt[,MTR])
plot(chr1_388_dt[, dup_reads], chr1_388_dt[,VAF], xlab = 'Nr duplicate reads', ylab = 'VAF')
plot(chr1_388_dt[, dup_reads], chr1_388_dt[,DEP])

plot(chr1_388_dt[, dup_reads], chr1_388_dt[,forward_mut])
plot(chr1_388_dt[, dup_reads], chr1_388_dt[,reverse_mut])

plot(chr1_388_dt[, dup_reads], chr1_388_dt[,forward_dep])
plot(chr1_388_dt[, dup_reads], chr1_388_dt[,reverse_dep])

plot(chr1_388_dt[, dup_reads], chr1_388_dt[,forward_vaf], xlab = 'Nr duplicate reads', ylab = 'VAF forward strand')
plot(chr1_388_dt[, dup_reads], chr1_388_dt[,reverse_vaf], xlab = 'Nr duplicate reads', ylab = 'VAF reverse strand')
plot(chr1_388_dt[, dup_reads], chr1_388_dt[,vaf_ratio], xlab = 'Nr duplicate reads', ylab = 'VAF forward / reverse')

plot(chr1_388_dt[, dup_reads], chr1_388_dt[,diff_dep])

# run correlation (few data points tho)
cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,MTR])
cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,VAF])
cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,DEP])

cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,forward_mut])
cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,reverse_mut])

cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,forward_dep])
cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,reverse_dep])

cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,forward_vaf]) # -0.23
cor.test(chr1_388_dt[, dup_reads], chr1_388_dt[,reverse_vaf]) # -0.25 

# there is a negative correlation between nr of duplicated reads and VAF 
# seems like it could be related that samples with more dup reads also have higher depth

###################################################################################################################################
# why is VAF ratio so off in some samples?
plot(chr1_388_dt[, dup_reads], chr1_388_dt[,vaf_ratio])
plot(chr1_388_dt[, DEP], chr1_388_dt[,vaf_ratio])
plot(chr1_388_dt[, forward_vaf], chr1_388_dt[,vaf_ratio])
plot(chr1_388_dt[, reverse_vaf], chr1_388_dt[,vaf_ratio])

plot(chr1_388_dt[, dep_ratio], chr1_388_dt[,VAF])
plot(chr1_388_dt[, dep_ratio], chr1_388_dt[,vaf_ratio])

# is it possible that some of the mutant reads are duplicates but not reported as such?
chr1_388_dt[VAF > 0.5]
# in all but one case, VAF is > 0.5 on both strands and relatively close to one another 
# VAF is > 0.5 for 5 samples of PD63383 (other than bb, which we know is contaminated) (0.43 in PD6338bb)

# for each sample, find out how many reads map to the exact same place
# and how many of those are marked as duplicates
flags = data.table()

for (s in samples_names){
  
  t = data.table(table(sam[sample==s, V4]))
  t[,V1:= as.numeric(V1)]
  pos_over1 = t[N>= 2, V1] %>% unlist()
  flags = rbind(flags, sam[sample==s & V4 %in% pos_over1])

}

# okay, now you have a dataframe where you have the original read and possible duplicates 
dim(flags)[1] # 313 reads, so this is original reads + all their possible duplicates
dim(flags[V2 >= 1024])[1] # 61 
# okay so only 61 of those are marked as possible duplicates 

# okay, let's get rid of reads that are marked as duplicates already 
flags_rm_marked_fq = flags[V2 < 512] # remove reads marked with 1024 (duplicates) and 512 (failed QC)
# this leaves us with 249 reads, some of which could still be duplicates 

# remove reads which are now single (in each sample)
flags_rm2 = data.table(flags_rm_marked_fq %>% group_by(sample) %>% filter(duplicated(V4) | duplicated(V4, fromLast=TRUE)))
# leaves us with 200 reads to inspect which could be possible duplicates 

# for each of those reads, check if they are wt of mutant
flags_rm2[, dist_to_mut := 38827952 - V4]
flags_rm2[, mut_site := substr(V10, dist_to_mut+1, dist_to_mut+1)]

table(flags_rm2[, mut_site]) # 0.565 is contributed by wt reads, oh that's not great right 

###################################################################################################################################
# Quality checks
# Analyse the distribution of flags (which could suggest specific issues in this region)
table(sam[, V2])
# 73 = 64 + 8 + 1 (3) - mate is unmapped
# 83 = 64 + 16 + 2 + 1 (203) - read is on the reverse strand
# 99 = 64 + 32 + 2 + 1 (233) - mate is on the reverse strand
# 121 = 64 + 32 + 16 + 8 + 1 (2) - mate is unmapped 
# 145 = 128 + 16 + 1 (1) - second read in the pair, on the reverse strand
# 147 = 128 + 16 + 2 + 1 (186) - 145 + mapped in a proper pair 
# 163 = 128 + 32 + 2 + 1 (194) - 147 but mate on the reverse strand 
# 595 = 512 + 64 + 16 + 2 + 1 (4) - fails QC
# 611 = 512 + 64 + 32 + 2 + 1 (2) - fails QC
# 659 = 512 + 128 + 16 + 2 + 1 (8) - fails QC 
# 675 = 512 + 128 + 32 + 2 + 1 (9) - failed QC 
# 1097 = 1024 + 64 + 8 + 1 - PCR duplicate 
# 1107 = 1024 + 64 + 16 + 2 + 1 - PCR duplicate 
# 1123 = 1024 + 64 + 32 + 2 + 1 - PCR duplicate 
# 1171 = 1024 + 128 + 16 + 2 + 1 - PCR duplicate 
# 1187 = 1024 + 128 + 32 + 2 + 1 - PCR duplicate 
# 1635 = 1024 + 512 + 64 + 32 + 2 + 1 - PCR duplicate + fails QC 
# 1683 = 1024 + 512 + 128 + 16 + 2 + 1 - PCR duplicate + fails QC 

# the 4 good flags are 83, 99, 147, 163 
sam[mut_site=='C' & V2 %in% c(83, 99, 147, 163)] # 511
sam[mut_site=='A' & V2 %in% c(83, 99, 147, 163)] # 302

sam[V2 %in% c(83, 99, 147, 163)] 
# 816 - so I am not sure why the total depth across samples is 751 
# I really don't know how one could resolve this atm 

# But then VAF 302/816 is not very different from VAF 275/751 (it is actually slightly higher)
# Therefore, throwing away good reads doesn't sound like a good explanation for this to me 

sam_hq = sam[V2 %in% c(83, 99, 147, 163)] 
table(sam[, sample])
table(sam_hq[, sample])

# can I identify all reads which seem okay (correct flags) and report the variant allele?
sam_mut_hq = sam[mut_site=='A' & V2 %in% c(83, 99, 147, 163)] # 302
sam_mut_pos = data.table(table(sam_mut_hq[, V4]))
hist(sam_mut_pos[, N], xlim = c(1,8))

p_duplicates = sam_hq[sam_hq[, .I[.N > 1], by = .(V4, sample, mut_site)]$V1] # 109 
table(p_duplicates[, mut_site]) # 0.39 so not much of a difference 

reads_site = sam_hq[sam_hq[, .I[.N > 1], by = .(sample, mut_site)]$V1] # 812

for (s in samples_names){
  dt = reads_site[sample==s]
  hist(dt[, dist_to_mut] %>% unlist(), main = glue('{s}'))
}

p_dup_seq_only = sam_hq[sam_hq[, .I[.N > 1], by = .(sample, V10, mut_site)]$V1] # 62 
table(p_dup_seq_only[, mut_site]) # 0.39 so not much of a difference 

p_dup_seq = sam_hq[sam_hq[, .I[.N > 1], by = .(V4, V10, sample, mut_site)]$V1] # 62 
table(p_dup_seq[, mut_site]) # 0.37 so not much of a difference 

p_dup_mate = sam_hq[sam_hq[, .I[.N > 1], by = .(V4, V9, sample, mut_site)]$V1] # 0 
table(p_dup_mate[, mut_site]) # no reads that would map like this  

p_dup_mate_seq = sam_hq[sam_hq[, .I[.N > 1], by = .(V4, V10, V9, sample, mut_site)]$V1] # 0 
table(p_dup_mate_seq[, mut_site]) # no reads that would map like this  

p_dup_id = p_duplicates[, V1] %>% unlist() # get rid of any remotely possible duplicates 

# could it be that these reads are not identified as PCR duplicates due to a shift or sth?

# can the CIGAR strings tell me anything?
# CIGAR positions: MNDI
# M = exact match of M positions
# N = alignment gap (next x positions on reference that don't match)
# D = deletion (next x positions on reference that don't match)
# I = insertion (next x positions on query that don't match)
# 2 As and 5 Cs, so nothing unusual (VAF 0.375 from this so nothing here)

sam_wt_hq = sam[mut_site=='C' & V2 %in% c(83, 99, 147, 163)] # 511
sam_wt_pos = data.table(table(sam_wt_hq[, V4]))
hist(sam_wt_pos[, N], xlim = c(1,8))

# okay what if I calculate VAF for each sample excluding those reads?
for (s in samples_names){
  vaf = dim(sam_hq[sample==s & !V1 %in% p_dup_id & mut_site=='A'])[1] / dim(sam_hq[sample==s & !V1 %in% p_dup_id])[1]
  print(paste('Sample:', s, 'VAF:', vaf))
  }

dim(sam_hq[sample %in% samples_normal_PD62341_clean & !V1 %in% p_dup_id & mut_site=='A'])[1] / dim(sam_hq[sample %in% samples_normal_PD62341_clean & !V1 %in% p_dup_id])[1]
# 0.31

dim(sam_hq[sample %in% samples_normal_PD63383_clean & !V1 %in% p_dup_id & mut_site=='A'])[1] / dim(sam_hq[sample %in% samples_normal_PD63383_clean & !V1 %in% p_dup_id])[1]
# 0.56 

###################################################################################################################################
# How many reads would have to be duplicates for the true VAF to be 0.5?

twins_filtered_dt[mut_ID == 'chr1_38827952_C_A', c('agg_normal_PD62341_clean_mtr', 
                                                   'agg_normal_PD62341_clean_dep',
                                                   'agg_normal_PD63383_clean_mtr', 
                                                   'agg_normal_PD63383_clean_dep'), with=F]

pd62341_mtr=56
pd62341_dep=199
pd63383_mtr=104
pd63383_dep=174

qbinom(p=c(0.05,0.95), size=pd62341_dep, prob=(pd62341_mtr/pd62341_dep))/pd62341_dep
qbinom(p=c(0.05,0.95), size=pd63383_dep, prob=(pd63383_mtr/pd63383_dep))/pd63383_dep

# say that due to PCR, there are 2 not-removed duplicates per sample, and all of those affect mutant reads
# would that explain the inflated VAF?
# there are 5 samples in each category so 2*5 = 10
qbinom(p=c(0.05,0.95), size=pd62341_dep-10, prob=((pd62341_mtr-10)/(pd62341_dep-10)))/(pd62341_dep-10)
qbinom(p=c(0.05,0.95), size=pd63383_dep-10, prob=((pd63383_mtr-10)/(pd63383_dep-10)))/(pd63383_dep-10)
# lower CI is still > 0.5

qbinom(p=c(0.05,0.95), size=pd62341_dep-15, prob=((pd62341_mtr-15)/(pd62341_dep-15)))/(pd62341_dep-15)
qbinom(p=c(0.05,0.95), size=pd63383_dep-15, prob=((pd63383_mtr-15)/(pd63383_dep-15)))/(pd63383_dep-15)

# imagine that 10 wt reads are thrown way (2 per sample)
qbinom(p=c(0.05,0.95), size=pd62341_dep+10, prob=((pd62341_mtr)/(pd62341_dep+10)))/(pd62341_dep+10)
qbinom(p=c(0.05,0.95), size=pd63383_dep+10, prob=((pd63383_mtr)/(pd63383_dep+10)))/(pd63383_dep+10)

qbinom(p=c(0.05,0.95), size=pd62341_dep+15, prob=((pd62341_mtr)/(pd62341_dep+15)))/(pd62341_dep+15)
qbinom(p=c(0.05,0.95), size=pd63383_dep+15, prob=((pd63383_mtr)/(pd63383_dep+15)))/(pd63383_dep+15)
# equally (un-reasonable) as the previous thing

###################################################################################################################################
# Can I just see the distribution of chr1_38825792 across samples?
twins_filtered_dt[mut_ID == 'chr1_38827952_C_A', c(samples_vaf), with=F]

# Is is the case that duplicated reads are more commonly A than C
table(sam[V2>512, mut_site]) # 0.4

###################################################################################################################################
# Can PCR duplicates specifically on the mutant strand in chr1_388 explain the inflated VAF?
qbinom(0.05, 174, 104/174)/174
# this is agg from 5 samples
# say each sample has 2 duplicate reads 
qbinom(0.05, 164, 94/164)/164
# okay this does actually get us to a reasonable number 
# but you would require > 2 duplicates per sample and your VAF is still really high  

qbinom(0.05, 199, 56/199)/199
qbinom(0.05, 189, 46/189)/189
# this would be 0.15 so still too high 
# again, you would have to have > 2 duplicates per sample for this to be viable
