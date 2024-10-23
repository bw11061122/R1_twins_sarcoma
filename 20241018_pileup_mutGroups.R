# Barbara Walkowiak bw18

# INPUT: merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# with filters 
# NOTE that for the moment, I am investigating if there is sth wrong in how the pileup was run
# I may alter the input df once everything is resolved

# NOTE ON ADDING COLUMNS WITH FILTERS 
# I want to add a column which says 0 if sample is fine and 1 if sample is wrong 

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
###################################################################################################################################
# Read the merged dataframe and list of mutations to include
setwd('/Users/bw18/Desktop/1SB/Analysis')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))
muts = read.csv('Data/mutations_include_20241018.txt')


######################################################################################################
# can we find mutations only present in one twin (look just at normal bc tumour is shared)
# turns out there are not really many of those that look real!
# we may want a different strategy for twin specific stuff ie which exact filters to use

twins_filtered_mtr[sum_normal_PD63383 == 0 & sum_normal_PD62341 == 4] # empty
twins_filtered_mtr[sum_normal_PD63383 == 0 & sum_normal_PD62341 == 3] # 4
# chr18_7097030_G_A a lot of insertions
# chr19_13472747_T_A a lot of insertions around
# chr20_54999342_C_A there are a few insertions but I'd think that it still looks okay

twins_filtered_mtr[sum_normal_PD63383 == 0 & sum_normal_PD62341 == 2] # 10
# chr11_24926209_A_G insertions so wouldn't necessarily believe it 
# chr1_12828564_A_G poor mapping quality 
# chr22_22531966_T_C insertions / poor mapping
# chr3_195982104_C_T poor mapping quality 
# chr3_62055057_C_G looks good!
# chr3_62055077_G_C looks good!
# chr5_22248315_A_G lots of insertions
# chr6_769976_C_T lots of insertions
# chr7_99888321_G_A poor mapping quality, looks tragic 
# chrX_77862341_C_A lots of deletions

twins_filtered_mtr[sum_normal_PD63383 == 6 & sum_normal_PD62341 == 0] # empty
twins_filtered_mtr[sum_normal_PD63383 == 5 & sum_normal_PD62341 == 0] # 2
# chr4_37131693_G_A lots of insertions
# chr5_23443714_T_G lots of insertions

twins_filtered_mtr[sum_normal_PD63383 == 4 & sum_normal_PD62341 == 0] # empty
twins_filtered_mtr[sum_normal_PD63383 == 3 & sum_normal_PD62341 == 0] # empty
# chr11_131162789_A_G deletions
# chr1_108488775_C_T very poor mapping
# chr2_54313444_T_G deletions
# chr9_6742174_T_A deletions + some reads poor mapping

# but then actually maybe what you should do is look for stuff called as germline in both twins
twins_germlinein1_mtr = twins_dt_filters[(f16_likelyGermline_PD62341 == 1 | f17_likelyGermline_PD63383 == 1) & sum_allfilters == 1,
                                         c(cols_info, cols_mtr, 'f16_likelyGermline_PD62341', 'f17_likelyGermline_PD63383'), with=FALSE] # 113 mutations called as germline in 1 twin
twins_germlinein1_vaf = twins_dt_filters[(f16_likelyGermline_PD62341 == 1 | f17_likelyGermline_PD63383 == 1) & sum_allfilters == 1,
                                         c('mut_ID', cols_info, cols_vaf, 'f16_likelyGermline_PD62341', 'f17_likelyGermline_PD63383'), with=FALSE] # 113 mutations called as germline in 1 twin
# looks bad unless specified below
# chr12_18097750_C_A looks kind of okay but quite low coverage in general 
# chr15_49480646_T_A looks ok, need to check how much i believe the germline assessment though
# chr15_55075744_C_G looks ok, need to check how much i believe the germline assessment though
# chr16_5479739_C_T looks very good
# chr17_33422229_C_A looks very good 
# chr17_80581632_C_T looks kind of okay, not sure about the germline tho
# chr18_12762524_C_A, chr18_12762527_C_A looks good, again not sure on the germline
# chr18_55426055_G_C looks good 
# chr1_157363251_C_T looks okay 
# chr1_16884197_C_G looks quite okay - would need to check
# chr1_5083107_A_C looks good
# chr1_66105222_T_C double check but maybe it's fine 
# chr1_76183938_A_G looks good
# chr20_171446_A_T looks ok
# chr2_55272895_A_G double check
# chr2_95662131_G_A double check
# chr3_50106043_C_T looks ok
# chr4_39717959_T_C looks ok
# chr5_61980224_A_G quite believable
# chr6_63596915_T_A mixed feelings

######################################################################################################
######################################################################################################

# classifying presence of mutation based on MTR
twins_mtr[,sum_normal := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_mtr[,sum_normal_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_mtr[,sum_normal_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]
twins_mtr[,sum_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_mtr[,sum_tumour_PD62341 := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_mtr[,sum_tumour_PD63383 := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_mtr[,sum_PD62341 := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_mtr[,sum_PD63383 := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]

# classifying presence of mutation based on VAF
samples_normal_VAF = paste(samples_normal, "VAF", sep="_")
samples_normal_PD62341_VAF = paste(samples_normal_PD62341, "VAF", sep="_")
samples_normal_PD63383_VAF = paste(samples_normal_PD63383, "VAF", sep="_")
samples_tumour_VAF = paste(samples_tumour, "VAF", sep="_")
samples_tumour_PD62341_VAF = paste(samples_tumour_PD62341, "VAF", sep="_")
samples_tumour_PD63383_VAF = paste(samples_tumour_PD63383, "VAF", sep="_")
samples_PD62341_VAF = paste(samples_PD62341, "VAF", sep="_")
samples_PD63383_VAF = paste(samples_PD63383, "VAF", sep="_")

twins_vaf[,sum_normal := rowSums(.SD>=4), .SDcols = samples_normal_VAF]
twins_vaf[,sum_normal_PD62341 := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_VAF]
twins_vaf[,sum_normal_PD63383 := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_VAF]
twins_vaf[,sum_tumour := rowSums(.SD>=4), .SDcols = samples_tumour_VAF]
twins_vaf[,sum_tumour_PD62341 := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_VAF]
twins_vaf[,sum_tumour_PD63383 := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_VAF]
twins_vaf[,sum_PD62341 := rowSums(.SD>=4), .SDcols = samples_PD62341_VAF]
twins_vaf[,sum_PD63383 := rowSums(.SD>=4), .SDcols = samples_PD63383_VAF]

# Filter dataframes with MTR, DEP and VAF
mut_exclude = twins_dt_filtered[sum_reqfilters==1, mut_ID] %>% unlist()
twins_mtr_filt = twins_mtr[!mut_ID %in% mut_exclude]
twins_dep_filt = twins_dep[!mut_ID %in% mut_exclude]
twins_vaf_filt = twins_vaf[!mut_ID %in% mut_exclude]

######################################################################################################
######################################################################################################
######################################################################################################
# CLASSES OF MUTATIONS OF INTEREST

# MUTATIONS PRIVATE TO 1 TWIN
# Presence in all samples ONLY from one twin
twins_mtr_filt[(sum_PD63383==0 & sum_PD62341==10)] # empty 
twins_mtr_filt[(sum_PD63383==8 & sum_PD62341==0)] # empty
twins_mtr_filt[(sum_PD63383==0 & sum_PD62341>=8)] # empty 
twins_mtr_filt[(sum_PD63383>=6 & sum_PD62341==0)] # empty
# but then the sums include tumour samples which we expect to be shared 

# Presence in all NORMAL samples ONLY from one twin
twins_mtr_filt[(sum_normal_PD63383==6 & sum_normal_PD62341==0)] # empty 
twins_mtr_filt[(sum_normal_PD63383==0 & sum_normal_PD62341==4)] # 2 

# Are these mutations also found in the tumour (from either or both twins?)
twins_mtr[(sum_normal_PD63383==0 & sum_normal_PD62341==4), sum_tumour] # 8 and 3 tumour samples 
twins_vaf[mut_ID == 'chr14_105458006_C_A'] # that's decent VAF but clearly not clonal
# "chr14_105458006_C_A": looks real!   
# "chr7_94130613_A_G": mapping okay, but lots of insertions / deletions (hard-to-seq region)

# "chr14_105458006_C_A"
# I CAN'T, THIS MUTATION IS IN A GENE CALLED MTA1 which is literally for metastasis-associated!
# on Jbrowse it looks like there could be a little bit of this in PD63383bb (skin) - MTR for this one is 3
# could that be due to contamination with tumour cells? or a seq / mapping artefact in this sample?
# mutation in MTA1 (intronic) - gene associated with metastasis 
# turns out expressed in the embryo: particularly nervous tissue https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4148129/ 

# Presence in SOME samples from one twin ONLY
twins_mtr_filt[(sum_PD63383==0 & sum_normal_PD62341>1)] # 19 (again, lots at low level - 1-4 reads)
# chr11_24926209_A_G - very limited confidence this is real 
# chr18_72837984_C_T - looks quite plausible
# chr20_54999342_C_A - looks quite plausible 

twins_mtr_filt[(sum_normal_PD63383>1 & sum_PD62341==0)] # 4 (again, lots at low level - 1-4 reads) 
# technically chr8_7421854_T_C looks vaguely plausible but not very confident 

twins_mtr_filt[(sum_normal_PD63383>1 & sum_normal_PD62341==0)] # 19 interestingly, in most cases, present in some tumour samples
# and actually more common in PD26341 tumour than PD63383 tumour
# which mutations of those look decent: literally none 

# Presence in all TUMOUR samples from one twin ONLY
twins_mtr_filt[(sum_tumour_PD63383==0 & sum_tumour_PD62341==6)] # 158 (can be present at low levels: MTR < 4)
twins_mtr_filt[(sum_tumour_PD63383==2 & sum_tumour_PD62341==0)] # empty

# Presence in MOST TUMOUR, and NO NORMAL samples (possible drivers?)
twins_mtr_filt[(sum_normal==0 & sum_tumour>=6)] # 14 mutations (but then you can imagine driver in the pre-cancer tissue too)
# GAL3ST1 - seems to be associated with several diseases (https://www.genecards.org/cgi-bin/carddisp.pl?gene=GAL3ST1#diseases)
# BID - seems to be involved in several adult cancers (https://www.genecards.org/cgi-bin/carddisp.pl?gene=BID#diseases)
# POLR3B - RNA Polymerase III subunit, associated with several developmental disorders (https://www.genecards.org/cgi-bin/carddisp.pl?gene=POLR3B&keywords=POLR3B#diseases)
# DAPK3 - there are disease associations for hematological cancer (https://www.genecards.org/cgi-bin/carddisp.pl?gene=DAPK3&keywords=DAPK3#diseases)
# KCND2 - K+ channel, mainly associated with dev disorders (cardiovascular) (https://www.genecards.org/cgi-bin/carddisp.pl?gene=KCND2&keywords=KCND2#diseases)

# Are there mutations only found in a tumour of one twin?
twins_mtr_filt[(sum_tumour_PD62341==0 & sum_tumour_PD63383==2)] # empty (NB also found in normal PD62341/PD63383 samples)
twins_mtr_filt[(sum_tumour_PD62341>=4 & sum_tumour_PD63383==0)] # 13,984 (NB also found in normal PD62341/PD63383 samples)

twins_mtr_filt[(sum_tumour_PD62341>=4 & sum_PD63383==0)] # 2
# only present in tumour samples from PD62341
# none of these look plausible to me 
# "chr8_106563310_A_G": lots of insertions (hard-to-seq region, string of As)
# "chr13_69469775_A_G": not convinced (lots of deletions in the region)

# Which tissues most often share mutations with the tumour?
twins_mtr_filt[(sum_tumour_PD62341>=2 & sum_tumour_PD62341>=1 & sum_normal_PD63383==0 & sum_normal_PD62341>0)] # 37 
# most often shared with 1 normal PD62341 

# which normal samples most commonly share mutations with tumour ones?
sum_PD62341_present = c()
sum_PD62341_absent = c()

for (name in samples_normal_PD62341_MTR){
  df = twins_mtr_filt[(sum_tumour_PD62341>=2 & sum_tumour_PD62341>=1 & sum_normal_PD63383==0 & sum_normal_PD62341>0), ..name]
  mut_present = sum(df>=4)
  mut_absent = sum(df==0)
  sum_PD62341_present = c(sum_PD62341_present, mut_present)
  sum_PD62341_absent = c(sum_PD62341_absent, mut_absent)
}

# also can do filtering out germline mutations (based on binomial test)
# first, select normal samples from MTR and DEP dataframes
twins_mtr_normal = twins_mtr_filt[, c('mut_ID', samples_normal_MTR), with=FALSE]

samples_normal_DEP = paste(samples_normal, "DEP", sep="_")
twins_dep_normal = twins_dep_filt[, c('mut_ID', samples_normal_DEP), with=FALSE]

# merge
twins_mtr_dep_normal = merge(twins_mtr_normal, twins_dep_normal)
twins_mtr_dep_normal[, sum_MTR := rowSums(.SD), .SDcols = 2:11]
twins_mtr_dep_normal[, sum_DEP := rowSums(.SD), .SDcols = 12:22]

# calculate the probability 
bin_test <- function(a, b, p = 0.5) {binom.test(a, b, 0.5, alternative=
                                                  c("less"), conf.level = 0.95)$p.value}
twins_mtr_dep_normal$p_val = mapply(bin_test, twins_mtr_dep_normal$sum_MTR, twins_mtr_dep_normal$sum_DEP)
twins_mtr_dep_normal = data.table(twins_mtr_dep_normal)
twins_mtr_dep_normal[, significance := as.factor(fcase( 
  p_val < 0.05/dim(twins_mtr_dep_normal)[1], 'significant',
  p_val >= 0.05/dim(twins_mtr_dep_normal)[1], 'not_significant'
))] 

paste('Number of mutations likely to not be germline (support from one-sided binomial):', dim(twins_mtr_dep_normal[significance=='significant'])[1])
# really looks like we have very, very few germline mutations in this dataset 
mut_likely_germline = twins_mtr_dep_normal[significance=='not_significant', mut_ID] %>% unlist()

# is it possible to get the likely VAF? or check if mutation is clonal in a subset of samples?

# how many mutations are present in normal samples?
twins_mtr_dep_normal[sum_MTR>=20] # 22,728 

# I think for each mutation, I'd like to know which samples it is present in 
# What are the most common combinations?

######################################################################################################
######################################################################################################
######################################################################################################


######################################################################################################
######################################################################################################
# ASCAT

# Read in ASCAT calls for tumour samples available
ascat_calls_PD62341ae = data.table(read.csv('Data/PD62341ae.ascat_ngs.summary.csv', header = F))
ascat_calls_PD62341ag = data.table(read.csv('Data/PD62341ag.ascat_ngs.summary.csv', header = F))
ascat_calls_PD62341aj = data.table(read.csv('Data/PD62341aj.ascat_ngs.summary.csv', header = F))
ascat_calls_PD62341ak = data.table(read.csv('Data/PD62341ak.ascat_ngs.summary.csv', header = F))
ascat_calls_PD62341am = data.table(read.csv('Data/PD62341am.ascat_ngs.summary.csv', header = F))
ascat_calls_PD62341ap = data.table(read.csv('Data/PD62341ap.ascat_ngs.summary.csv', header = F))
ascat_calls_PD63383ap = data.table(read.csv('Data/PD63383ap.ascat_ngs.summary.csv', header = F))
ascat_calls_PD63383aq = data.table(read.csv('Data/PD63383aq.ascat_ngs.summary.csv', header = F))

ascat_calls_PD62341ae[, sample := 'PD62341ae']
ascat_calls_PD62341ag[, sample := 'PD62341ag']
ascat_calls_PD62341aj[, sample := 'PD62341aj']
ascat_calls_PD62341ak[, sample := 'PD62341ak']
ascat_calls_PD62341am[, sample := 'PD62341am']
ascat_calls_PD62341ap[, sample := 'PD62341ap']
ascat_calls_PD63383ap[, sample := 'PD63383ap']
ascat_calls_PD63383aq[, sample := 'PD63383aq']

ascat_calls = data.table(rbind(ascat_calls_PD62341ae, ascat_calls_PD62341ag, ascat_calls_PD62341aj,
                               ascat_calls_PD62341ak, ascat_calls_PD62341am, ascat_calls_PD62341ap,
                               ascat_calls_PD63383ap, ascat_calls_PD63383aq))
setnames(ascat_calls, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8'), c('ID', 'chr', 'start', 'end', 'ploidy_ref', 'minorallele_ref', 'poidy_tumour', 'minorallele_tumour'))
ascat_calls[, mut_ID := paste(chr, start, end, sep='_')]
paste('Number of copy number changes identified:', length(ascat_calls[, mut_ID] %>% unlist() %>% unique()))

# sort by chromosome 
ascat_calls = ascat_calls[order(ascat_calls[,chr],ascat_calls[,start],ascat_calls[,end])]

# I find it quite likely that the whole thing would be lost really 
# Why is the first line called if it is the same in sample and tumour?

