###################################################################################################################################
# SUPPLEMENTARY SCRIPT

# Script to check for the presence of possible driver mutations
# November - December 2024
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that pass QC filters (germline or somatic)
# 3 list of germline mutations 
# 4 list of driver genes involved in fibromatoses (from Henry Lee-Six)

# OUTPUT:
# 1 list of possible driver genes (if any found)

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

###################################################################################################################################
# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Import list of driver genes (from Henry Lee-Six, 12/11/2024)
driver_genes_dt = fread('Data/Drivers/HLS_fibromatoses_driver_list_with_fusions.csv', header=T)
driver_genes = driver_genes_dt[, gene] %>% unlist()

# Genes identified through my google / literature search (as associated with sarcomas)
driver_genes_search = c('TP53', 'RB', 'ATRX1', 'PTEN', 'AKT', 'DKK1', 'PRDM10', 'CDKN2A', 'TRIO', 
                                         'BCOR', 'CCNB3', 'DUX4', 'SMARCB1', 'SMARCA4', 'MTOR', 'MYCN', 'CRKL',
                                         'PDGFRA', 'MLL', 'GSP2', 'PIK3CA', 'PIK3CG', 'TSC1', 'TSC2', 'TERT', 'FOXO1',
                                         'PAX3A', 'PAX7A', 'SSX1', 'SSX2', 'ELP1', 'CHEK2', 'ERCC4', 'RAD51', 'NF1', 'NF2')
driver_genes = c(driver_genes, driver_genes_search)

# Import list of mutations that passed QC (germline or somatic)
muts_qc = read.table("Out/F1/F1_mutations_passedQuality_20241208.txt")
muts_qc = muts_qc$V1 %>% unlist()
paste('Number of mutations that passed QC:', length(muts_qc)) # 334,043

# Import list of likely germline mutations (based on exact binomial)
muts_germline = read.table("Out/F1/F1_mutations_putativeGermline_20241208.txt")
muts_germline = muts_germline$V1 %>% unlist()
paste('Number of likely germline mutations:', length(muts_germline)) # 332,974

# Import indel data 
twins_dt_indels = fread('Data/3434_indels_20241112.txt', sep = '\t') # import dataset with indels 
# to see how this file was generated, search README.txt under the data 12/11/2024

# Import rearrangement calls from BRASS (see README.txt for the data location on the farm and file transfer script) 
rearrangements = list.files("Data/Brass/", full.names = TRUE, pattern = "removed_header.bedpe")
rearr_calls = rbindlist(lapply(rearrangements, fread), fill = TRUE)

# Subset the SNV dt to only include mutations that passed basic QC
twins_dt_filters = twins_dt[mut_ID %in% muts_qc]

######################################################################################################
# PLOT SETTINGS

# Specify colors for plotting 
col_tumour = '#ad0505'
col_normal = '#07a7d0'
col_PD62341 = "#0ac368"
col_PD63383 = "#a249e8"

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
paste('Number of protein-altering mutations:', length(muts_coding_subs)) # 1354 

# search across any mutations: look for genes and protein coding altering mutations 
genes_coding_subs = sort(twins_dt_filters[Gene != '-' & Effect %in% c('ess_splice', 'missense', 'nc_variant', 'nonsense', 'start_lost', 'stop_lost'), Gene] %>% unique() %>% unlist())
paste('Number of genes with protein-altering mutations:', length(genes_coding_subs)) # 1128 

# search across mutations that pass quality filters (regardless of presence in germline)
genes_missense_subs = twins_dt_filters[Gene != '-' & Effect %in% c('missense', 'nonsense'), Gene] %>% unique() %>% unlist()
paste('Number of genes with missense mutations:', length(genes_missense_subs)) # 1094

# intersect the list of genes with the list of genes from Henry 
genes_possibleDriverSubs = Reduce(intersect, list(genes_missense_subs, driver_genes))
paste('Number of genes with mis-sense mutations present in the driver gene list:', length(genes_possibleDriverSubs)) # 42
# for POLE, POLD1, we would have likely seen signature mutations in the genome 

# look at how those look like 
# SMARCA1 (https://pmc.ncbi.nlm.nih.gov/articles/PMC6513346/)
twins_dt_filters[Gene == 'SMARCA1', c('mut_ID', 'Gene', 'Effect'), with=FALSE] # one mis-sense, no hotspots
# KMT2A (known leukemia gene)
twins_dt_filters[Gene == 'KMT2A', c('mut_ID', 'Gene', 'Effect'), with=FALSE] # one mis-sense, no hotspots

# are all these mutations classified as germline?
muts_in_possibleDrivers = twins_dt_filters[Gene %in% genes_possibleDriverSubs, mut_ID] %>% unlist()
sum(muts_in_possibleDrivers %in% muts_germline) / length(muts_in_possibleDrivers) # almost all classified as germline
twins_dt_filters[mut_ID %in% setdiff(muts_in_possibleDrivers, muts_germline), c('mut_ID', 'Gene', 'Effect', samples_normal_vaf), with=FALSE] # one mis-sense, no hotspots
# CARD11, intronic, NCOA2, 2 intronic mutations (side by side); none of those passed final QC (not included in the final dataset)

######################################################################################################
# SEARCH FOR PROTEIN CODING MUTATIONS IN INDELS

setnames(twins_dt_indels, c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10',  'V11'),
         c('Chrom', 'pos', 'rm1', 'Ref', 'Alt', 'rm2', 'Filter', 'Info', 'rm3', 'rm4', 'rm5'))

twins_dt_indels = twins_dt_indels[, c('Chrom', 'pos', 'Ref', 'Alt', 'Filter', 'Info'), with=FALSE]
twins_dt_indels[, Gene := tstrsplit(Info, 'VD=', fixed = T, keep=2)]
twins_dt_indels[, Gene := tstrsplit(Gene, '|', fixed = T, keep=1)]

# search across any mutations: look for genes and protein coding altering mutations (we are ignoring substitutions)
genes_coding_indels = sort(twins_dt_indels[, Gene] %>% unique() %>% unlist())
paste('Number of genes with protein-altering indels:', length(genes_coding_indels)) # 83 

# intersect with the driver gene list from Henry 
genes_possibleDriverIndels = Reduce(intersect, list(genes_coding_indels, driver_genes)) # 2 
# "ARID1B" "NFKBIE"

######################################################################################################
# SEARCH FOR REARRANGEMENTS AFFECTING PROTEIN-CODING GENES 
setnames(rearr_calls, paste0('V', 1:46), c('chr1',	'start1','end1',	'chr2', 'start2',	'end2',	'id/name',	'brass_score',
                                           'strand1',	'strand2',	'sample',	'svclass',	'bkdist',	'assembly_score',	'readpair names',	'readpair count',	'bal_trans',	
                                           'inv',	'occL',	'occH',	'copynumber_flag',	'range_blat',	'Brass Notation',	'non-template',	'micro-homology	assembled', 'readnames	assembled',
                                           'read count',	'gene1',	'gene_id1', 'transcript_id1',	'strand1',	'end_phase1',	'region1',	'region_number1',	'total_region_count1',
                                           'first/last1',	'gene2',	'gene_id2',	'transcript_id2',	'strand2',	'phase2',	'region2',	'region_number2', 'total_region_count2', 'first/last2', 'fusion_flag'))

rearr_calls[gene1!='_', c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'sample', 'svclass', 'gene1', 'gene2'), with=FALSE]

gene1 = rearr_calls[gene1!='_', c('gene1'), with=FALSE] %>% unlist() %>% unique()
gene2 = rearr_calls[gene1!='_', c('gene2'), with=FALSE] %>% unlist() %>% unique()
genes_rearr = c(gene1, gene2) %>% unique()
paste('Number of genes involved in rearrangements:', length(genes_rearr)) # 43

sum(genes_rearr %in% driver_genes) # 3 
genes_rearr[genes_rearr %in% driver_genes]
# MN1  
# LRP1B # only PD62341b, 4 reads, with SCAF11, 200 (Intron to something else fusion)
# PTPRD # only PD62341ak, 9 reads (both gene1 and gene2 is PTPRB), 840 (5UTR to 5UTR fusion, same gene, same orientation)

# Look up the MN1 gene (which we know is involved in the fusion with ZNF341; confirmed from clinical data + Nathan reconstructed rearrangement)
rearr_calls[gene1=='MN1', c('sample', 'gene1', 'gene2'), with=FALSE] 
# why is there absolutely nothing in gene2??
# all tumour samples + skin + liver + heart (weird? - there are 4 read pair counts in the heart for this, cf 23 in clean tumour sample)
# MN1 seems to be interpreted by BRASS as a deletion rather than a rearrangement, not sure why 



