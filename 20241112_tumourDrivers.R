# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-11-12
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# Searching for drivers: substitutions and indels  

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

# Import substitution data  
setwd('/Users/bw18/Desktop/1SB')
twins_dt_subs = data.table(read.csv('Data/twins_dt_filters_20241112.csv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt_subs), value = TRUE)
twins_dt_subs[, c(twins_PDv38is) := NULL]

# Filter to only include mutations retained for filtering 
muts = read.table('Data/mutations_include_20241106_1002.txt') %>% unlist()
paste('Number of mutations that passed required filters:', length(muts)) # 1002
twins_dt_filtered = twins_dt_subs[mut_ID %in% muts]

# Import indel data 
twins_dt_indels = data.table(read.table('Data/3434_indels_20241112.txt')) # import dataset with indels 
# to see how this file was generated, search README.txt under the data 12/11/2024

# Import rearrangement calls (see README.txt for the data location on the farm and file transfer script) 
rearrangements <- list.files("Data", full.names = TRUE, pattern = "removed_header.bedpe")
rearr_calls <- rbindlist(lapply(rearrangements, fread), fill = TRUE)

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

samples_mtr = grep("MTR", names(twins_dt_filters), value = TRUE)
samples_dep = grep("DEP", names(twins_dt_filters), value = TRUE)
samples_vaf = grep("VAF", names(twins_dt_filters), value = TRUE)

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
# SEARCH FOR PROTEIN CODING MUTATIONS IN ALL SUBSTITUTIONS PRESENT HERE 

# how many mutations do I find?
muts_coding_subs = sort(twins_dt_filters[Gene != '-' & Effect %in% c('ess_splice', 'missense', 'nc_variant', 'nonsense', 'start_lost', 'stop_lost'), mut_ID] %>% unlist())
paste('Number of protein-altering mutations:', length(muts_coding_subs)) # 1502 

# search across any mutations: look for genes and protein coding altering mutations (we are ignoring substitutions)
genes_coding_subs = sort(twins_dt_filters[Gene != '-' & Effect %in% c('ess_splice', 'missense', 'nc_variant', 'nonsense', 'start_lost', 'stop_lost'), Gene] %>% unique() %>% unlist())
paste('Number of genes with protein-altering mutations:', length(genes_coding_subs)) # 1230 

# search across mutations that pass quality filters (regardless of presence in germline)
genes_coding_qc_subs = twins_dt_filters[sum_req_filters_qual == 0 & Gene != '-' & Effect %in% c('ess_splice', 'missense', 'nonsense', 'start_lost', 'stop_lost'), Gene] %>% unique() %>% unlist()
paste('Number of genes with protein-altering mutations (mutations pass QC filter):', length(genes_coding_qc_subs)) # 1106

genes_missense_qc_subs = twins_dt_filters[sum_req_filters_qual == 0 & Gene != '-' & Effect %in% c('missense', 'nonsense'), Gene] %>% unique() %>% unlist()
paste('Number of genes with missense mutations (mutations pass QC filter):', length(genes_missense_qc_subs)) # 1087

# genes that are interesting in my dataset (based on my knowledge)
# APC2
# BCLC2L10, BCL2L12
# CDK12, CUL9
# DYRK2 
# KMT2A
# MYO genes (MYO3B, 5B, 5C)
# PARP 
# POLE
# RARA
# SMARCA1

# look at how those look like 
twins_dt_filters[Gene == 'SMARCA1', c('mut_ID', 'Gene', 'Effect', 'f7_likelyGermline_bothTwins'), with=FALSE] ## https://pmc.ncbi.nlm.nih.gov/articles/PMC6513346/?
twins_dt_filters[Gene == 'KMT2A', c('mut_ID', 'Gene', 'Effect', 'f7_likelyGermline_bothTwins'), with=FALSE]

######################################################################################################
# SEARCH FOR PROTEIN CODING MUTATIONS IN ALL INDELS PRESENT HERE 

setnames(twins_dt_indels, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10',  'V11'),
         c('Chrom', 'pos', 'rm1', 'Ref', 'Alt', 'rm2', 'Filter', 'Info', 'rm3', 'rm4', 'rm5'))

twins_dt_indels = twins_dt_indels[, c('Chrom', 'pos', 'Ref', 'Alt', 'Filter', 'Info'), with=FALSE]
twins_dt_indels[, Gene := tstrsplit(Info, 'VD=', fixed = T, keep=2)]
twins_dt_indels[, Gene := tstrsplit(Gene, '|', fixed = T, keep=1)]

# search across any mutations: look for genes and protein coding altering mutations (we are ignoring substitutions)
genes_coding_indels = sort(twins_dt_indels[, Gene] %>% unique() %>% unlist())
paste('Number of genes with protein-altering indels:', length(genes_coding_indels)) # 83 

# Looked at genes associated with indels: nothing clear that I can see from that 

# Genes hit by both a mutation and indel 
Reduce(intersect, list(genes_coding_qc_subs, genes_coding_indels))
# "CCDC63", "TSSK4", "CEP152", "PIF1", "ADNP2", "OR2T11", "OR5AC2", "VWDE"

# Are all those mutations germline?
table(twins_dt_filters[Gene %in% genes_coding_qc_subs, f7_likelyGermline_bothTwins]) # check number of mutations 
# 0: 297, 1: 22613

twins_dt_filters[Gene %in% genes_coding_qc_subs & f7_likelyGermline_bothTwins==0 & Effect == 'missense', c('mut_ID', 'Gene','Effect')]
# OR52I1, AHNAK, CCL4L2, AMY2A, PRAMEF11, GSTT4, TTN, MUC7, MUC12, DPP6, CYP11B2, OSR2
# checked all on Jbrowse, most don't look convincing 
# CCL4L2 mutation looks good on Jbrowse 
# MUC7 looks real and fun - have very clear mutations on two chromosomes (one on each)

twins_dt_filters[Gene %in% genes_coding_qc_subs & f7_likelyGermline_bothTwins==1 & Effect == 'missense', c('mut_ID', 'Gene','Effect')]

######################################################################################################
# GENES THAT COME UP IN GOOGLE SEARCH FOR GENES RELATED TO UNDIFFERENTIATED SARCOMAS
possible_genes = c('TP53', 'RB', 'ATRX1', 'PTEN', 'AKT', 'DKK1', 'PRDM10', 'CDKN2A', 'TRIO', 
                   'BCOR', 'CCNB3', 'DUX4', 'SMARCB1', 'SMARCA4', 'MTOR', 'MYCN', 'CRKL',
                   'PDGFRA', 'MLL', 'GSP2', 'PIK3CA', 'PIK3CG', 'TSC1', 'TSC2', 'TERT', 'FOXO1',
                   'PAX3A', 'PAX7A', 'SSX1', 'SSX2', 'ELP1', 'CHEK2', 'ERCC4', 'RAD51', 'NF1', 'NF2')

possible_genes %in% c(genes_coding_subs, genes_coding_indels)
# nothing comes up 

######################################################################################################
# SEARCHED THROUGH THE LIST OF GENES I RECEIVED FROM HENRY 

# file received from Henry Lee-Six, 12/11/2024
driver_genes_dt = data.table(read.csv('Data/HLS_fibromatoses_driver_list_with_fusions.csv', header=T))
driver_genes = driver_genes_dt[, gene] %>% unlist()

genes_muts_drivers = Reduce(intersect, list(c(genes_coding_subs, genes_coding_indels), driver_genes))
genes_muts_drivers = c(genes_muts_drivers, 'SMARCA1') # add SMARCA1 just in case 
paste('Number of mutated genes that could be plausible drivers:', length(genes_muts_drivers)) # 47

# for all of those genes, should inspect mutations on Jbrowse and see if VAF makes sense 
twins_dt_filters[Gene %in% genes_muts_drivers & Effect == 'nonsense', c('mut_ID', 'Gene', 'Effect', samples_vaf), with=FALSE] # no non-sense mutations 
twins_dt_filters[Gene %in% genes_muts_drivers & Effect == 'missense', c('mut_ID', 'Gene', 'Effect', samples_vaf), with=FALSE]
# looks okay unless indicated otherwise 
# "chr10_75029584_C_T"  
# "chr10_87364920_G_A" # poor mapping  
# "chr11_118471877_G_A" 
# "chr12_124344759_A_T" 
# "chr12_132676107_T_C" 
# "chr14_92809378_G_A" 
# "chr16_10901532_C_T"  
# "chr16_72797127_T_G"  
# "chr17_35844923_G_A" 
# "chr17_39526140_C_T"  
# "chr17_40331346_C_G"  
# "chr17_80393425_A_C" 
# "chr19_49665963_G_T"  
# "chr19_49673722_G_A"  
# "chr19_50402053_G_A"  
# "chr19_8885784_A_C"   
# "chr19_8902144_G_A" # poor mapping 
# "chr19_8915528_A_T" # poor mapping  
# "chr19_8915529_C_A" # poor mapping   
# "chr19_8964561_T_C"   
# "chr1_119922665_G_A"  
# "chr1_23559344_C_T" # lost in the tumour  
# "chr1_36468114_A_G" #   
# "chr1_47251573_C_T" # lost in the tumour 
# "chr1_77956689_T_C" # lost in the tumour  
# "chr20_20052468_G_T"  
# "chr22_19183443_C_G"  
# "chr22_35085668_G_T"  
# "chr2_108765265_A_T" # poor mapping  
# "chr2_140769245_G_C" 
# "chr2_140950243_G_C"  
# "chr2_47410107_A_G"   
# "chr3_195782474_C_T" # poor mapping   
# "chr3_47122093_G_A"   
# "chr4_1804404_T_C"   
# "chr4_80055981_T_C"  
# "chr5_150077432_C_A"  
# "chr5_180603325_C_T"  
# "chr6_117344068_T_C"  
# "chr6_138813286_C_A"  
# "chr7_2930035_C_G"    
# "chr7_2930064_G_T"   
# "chr7_92079299_A_G"   
# "chr8_112685531_G_A"  
# "chr8_32647871_G_A"   
# "chr8_41664889_G_C"   
# "chr8_70128479_T_C"   
# "chr9_115082732_T_C" 
# "chr9_121140773_A_G"  
# "chr9_133131037_G_A"  
# "chr9_95516658_C_G"   
# "chrX_129511839_G_A"  
# "chrX_41364334_G_A" 

######################################################################################################
# Get through rearrangement calls
setnames(rearr_calls, paste0('V', 1:46), c('chr1',	'start1','end1',	'chr2', 'start2',	'end2',	'id/name',	'brass_score',
  'strand1',	'strand2',	'sample',	'svclass',	'bkdist',	'assembly_score',	'readpair names',	'readpair count',	'bal_trans',	
  'inv',	'occL',	'occH',	'copynumber_flag',	'range_blat',	'Brass Notation',	'non-template',	'micro-homology	assembled', 'readnames	assembled',
  'read count',	'gene1',	'gene_id1', 'transcript_id1',	'strand1',	'end_phase1',	'region1',	'region_number1',	'total_region_count1',
  'first/last1',	'gene2',	'gene_id2',	'transcript_id2',	'strand2',	'phase2',	'region2',	'region_number2', 'total_region_count2', 'first/last2', 'fusion_flag'))

rearr_calls[gene1!='_', c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'sample', 'svclass', 'gene1', 'gene2'), with=FALSE]

rearr_calls[gene1=='MN1', c('sample', 'gene1', 'gene2'), with=FALSE] # why is there absolutely nothing in gene2??
# all tumour samples + skin + liver + heart (weird? - there are 4 read pair counts in the heart for this, cf 23 in clean tumour sample)

gene1 = rearr_calls[gene1!='_', c('gene1'), with=FALSE] %>% unlist() %>% unique()
gene2 = rearr_calls[gene1!='_', c('gene2'), with=FALSE] %>% unlist() %>% unique()
genes_rearr = c(gene1, gene2) %>% unique()
paste('Number of genes involved in rearrangements:', length(genes_rearr)) # 43

sum(genes_rearr %in% driver_genes) # 3 
# MN1  
# LRP1B # only PD62341b, 4 reads, with SCAF11, 200 Intron to something else fusion
# PTPRD # only PD62341ak, 9 reads (both gene1 and gene2 is PTPRB), 840 5UTR to 5UTR fusion, same gene. Same orientation
 
# this thing interprets MN1 as a deletion???? but we know from Nathan that there is a rearrangement 



