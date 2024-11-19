###################################################################################################################################
# SCRIPT 4

# Script to analyse the pileup (high quality, run 15/10/2024)
# 2024-11-06
# Barbara Walkowiak bw18

# INPUT: 
# 1 merged pileup dataframes with mutation calls (from CaVEMan) for tumour and normal samples
# 2 list of mutations that passed required filters 

# OUTPUT:
# 1 lists of mutations in specific categories of interest for the phylogeny
# 2 plots for each category of mutations of interest 

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

# Read the merged dataframe 
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv')) # import high quality pileup

# Drop columns with PD38is_wgs (used as reference)
twins_PDv38is = grep("PDv38is", names(twins_dt), value = TRUE)
twins_dt[, c(twins_PDv38is) := NULL]

# Load QC-validated mutations (final set of mutations to be used)
muts_dt = data.table(read.csv('Data/20241114_599muts_QCJbrowse.csv', header = T))
muts = muts_dt[Jbrowse.quality == 'Y', mut_ID] %>% unlist()
paste('Number of mutations that passed QC:', length(muts)) # 258 

# Create a dataframe with mutations of interested retained
twins_filtered_dt = twins_dt[mut_ID %in% muts]

# Import dataframe with purity estimates
purity_dt = data.table(read.csv('Data/20241114_estimates_tumour_cont_27muts_median.csv'))

# Import list of driver genes (from Henry Lee-Six, 12/11/2024)
driver_genes_dt = data.table(read.csv('Data/HLS_fibromatoses_driver_list_with_fusions.csv', header=T))
driver_genes = driver_genes_dt[, gene] %>% unlist()

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
col_PD62341_spleen = "#148259"
col_PD63383_spleen = "#53088e"
col_PD63383_skin = '#d0bbe1'
col_bar = '#e87811'

######################################################################################################
# SAMPLES

# Create lists of possible samples of interest
samples_names = c("PD62341v", "PD62341q", "PD62341aa", "PD62341ad", "PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak", "PD63383bb", "PD62341ae", "PD62341ag", "PD62341aj", "PD62341ak", "PD62341am", "PD62341ap",  "PD63383ap", "PD63383aq", "PD62341b", "PD62341h", "PD62341n", "PD62341u")
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
# CALCULATE VAF FOR AGGREGATED SAMPLES

# Aggregated VAF in tumour samples 
twins_filtered_dt[, dep_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_dep]
twins_filtered_dt[, dep_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_dep]
twins_filtered_dt[, dep_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_dep]
twins_filtered_dt[, dep_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_dep]
twins_filtered_dt[, dep_all_normal := rowSums(.SD), .SDcols = samples_normal_dep]
twins_filtered_dt[, dep_all_tumour := rowSums(.SD), .SDcols = samples_tumour_dep]
twins_filtered_dt[, dep_all := rowSums(.SD), .SDcols = samples_dep]

twins_filtered_dt[, mtr_all_normal_PD62341 := rowSums(.SD), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, mtr_all_normal_PD63383 := rowSums(.SD), .SDcols = samples_normal_PD63383_mtr]
twins_filtered_dt[, mtr_all_tumour_PD62341 := rowSums(.SD), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, mtr_all_tumour_PD63383 := rowSums(.SD), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, mtr_all_normal := rowSums(.SD), .SDcols = samples_normal_mtr]
twins_filtered_dt[, mtr_all_tumour := rowSums(.SD), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, mtr_all := rowSums(.SD), .SDcols = samples_mtr]

twins_filtered_dt[, vaf_all_normal_PD62341 := mtr_all_normal_PD62341 / dep_all_normal_PD62341]
twins_filtered_dt[, vaf_all_normal_PD63383 := mtr_all_normal_PD63383 / dep_all_normal_PD63383]
twins_filtered_dt[, vaf_all_tumour_PD62341 := mtr_all_tumour_PD62341 / dep_all_tumour_PD62341]
twins_filtered_dt[, vaf_all_tumour_PD63383 := mtr_all_tumour_PD63383 / dep_all_tumour_PD63383]
twins_filtered_dt[, vaf_all_normal := mtr_all_normal / dep_all_normal]
twins_filtered_dt[, vaf_all_tumour := mtr_all_tumour / dep_all_tumour]
twins_filtered_dt[, vaf_all := mtr_all / dep_all]

# Add presence / absence based on MTR data
twins_filtered_dt[, sum_all_mtr := rowSums(.SD>=4), .SDcols = c(samples_tumour_mtr, samples_normal_mtr)]
twins_filtered_dt[, sum_tumour_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_mtr]
twins_filtered_dt[, sum_normal_mtr := rowSums(.SD>=4), .SDcols = samples_normal_mtr]
twins_filtered_dt[, sum_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_PD62341_mtr]
twins_filtered_dt[, sum_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_PD63383_mtr]
twins_filtered_dt[, sum_tumour_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD62341_mtr]
twins_filtered_dt[, sum_tumour_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_tumour_PD63383_mtr]
twins_filtered_dt[, sum_normal_PD62341_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_mtr]
twins_filtered_dt[, sum_normal_PD63383_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_mtr]

# Add presence / absence based on VAF data
twins_filtered_dt[, sum_all_vaf := rowSums(.SD>=0.1), .SDcols = c(samples_tumour_vaf, samples_normal_vaf)]
twins_filtered_dt[, sum_tumour_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_vaf]
twins_filtered_dt[, sum_normal_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_vaf]
twins_filtered_dt[, sum_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD62341_vaf]
twins_filtered_dt[, sum_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_PD63383_vaf]
twins_filtered_dt[, sum_tumour_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD62341_vaf]
twins_filtered_dt[, sum_tumour_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_tumour_PD63383_vaf]
twins_filtered_dt[, sum_normal_PD62341_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_vaf]
twins_filtered_dt[, sum_normal_PD63383_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_vaf]

twins_filtered_dt[, sum_all_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_all_mtr', 'sum_all_vaf')]
twins_filtered_dt[, sum_tumour_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_mtr', 'sum_tumour_vaf')]
twins_filtered_dt[, sum_normal_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_mtr', 'sum_normal_vaf')]
twins_filtered_dt[, sum_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_PD62341_mtr', 'sum_PD62341_vaf')]
twins_filtered_dt[, sum_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_PD63383_mtr', 'sum_PD63383_vaf')]
twins_filtered_dt[, sum_tumour_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_PD62341_mtr', 'sum_tumour_PD62341_vaf')]
twins_filtered_dt[, sum_tumour_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_tumour_PD63383_mtr', 'sum_tumour_PD63383_vaf')]
twins_filtered_dt[, sum_normal_PD62341_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD62341_mtr', 'sum_normal_PD62341_vaf')]
twins_filtered_dt[, sum_normal_PD63383_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD63383_mtr', 'sum_normal_PD63383_vaf')]

# Add a section for clean normal PD63383 samples (from script 2, evident that skin (PD63383bb) is contaminated)
samples_normal_PD63383_clean = c("PD63383w", "PD63383t", "PD63383u", "PD63383ae", "PD63383ak")
samples_normal_PD63383_clean_dep = paste0(samples_normal_PD63383_clean, '_DEP')
samples_normal_PD63383_clean_mtr = paste0(samples_normal_PD63383_clean, '_MTR')
samples_normal_PD63383_clean_vaf = paste0(samples_normal_PD63383_clean, '_VAF')

twins_filtered_dt[, sum_normal_PD63383_clean_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD63383_clean_mtr]
twins_filtered_dt[, sum_normal_PD63383_clean_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD63383_clean_vaf]
twins_filtered_dt[, sum_normal_PD63383_clean_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD63383_clean_mtr', 'sum_normal_PD63383_clean_vaf')]

twins_filtered_dt[, agg_normal_PD63383_clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD63383_clean_mtr]
twins_filtered_dt[, agg_normal_PD63383_clean_dep := rowSums(.SD), .SDcols = samples_normal_PD63383_clean_dep]
twins_filtered_dt[, agg_normal_PD63383_clean_vaf := agg_normal_PD63383_clean_mtr / agg_normal_PD63383_clean_dep]

# Add a section for clean normal PD62341 samples (exclude spleen which we know got transfusion)
samples_normal_PD62341_clean = c("PD62341n", "PD62341h", "PD62341ad", "PD62341aa", "PD62341q")
samples_normal_PD62341_clean_dep = paste0(samples_normal_PD62341_clean, '_DEP')
samples_normal_PD62341_clean_mtr = paste0(samples_normal_PD62341_clean, '_MTR')
samples_normal_PD62341_clean_vaf = paste0(samples_normal_PD62341_clean, '_VAF')

twins_filtered_dt[, sum_normal_PD62341_clean_mtr := rowSums(.SD>=4), .SDcols = samples_normal_PD62341_clean_mtr]
twins_filtered_dt[, sum_normal_PD62341_clean_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_PD62341_clean_vaf]
twins_filtered_dt[, sum_normal_PD62341_clean_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_PD62341_clean_mtr', 'sum_normal_PD62341_clean_vaf')]

twins_filtered_dt[, agg_normal_PD62341_clean_mtr := rowSums(.SD), .SDcols = samples_normal_PD62341_clean_mtr]
twins_filtered_dt[, agg_normal_PD62341_clean_dep := rowSums(.SD), .SDcols = samples_normal_PD62341_clean_dep]
twins_filtered_dt[, agg_normal_PD62341_clean_vaf := agg_normal_PD62341_clean_mtr / agg_normal_PD62341_clean_dep]

# Add a section to identify contaminated normal samples 
samples_normal_contaminated = c("PD62341aa", "PD62341h", "PD63383bb")
samples_normal_contaminated_dep = paste0(samples_normal_contaminated, '_DEP')
samples_normal_contaminated_mtr = paste0(samples_normal_contaminated, '_MTR')
samples_normal_contaminated_vaf = paste0(samples_normal_contaminated, '_VAF')

twins_filtered_dt[, sum_normal_contaminated_mtr := rowSums(.SD>=4), .SDcols = samples_normal_contaminated_mtr]
twins_filtered_dt[, sum_normal_contaminated_vaf := rowSums(.SD>=0.1), .SDcols = samples_normal_contaminated_vaf]
twins_filtered_dt[, sum_normal_contaminated_mtr_vaf := apply(.SD, 1, function(x) min(x)), .SDcols = c('sum_normal_contaminated_mtr', 'sum_normal_contaminated_vaf')]

twins_filtered_dt[, agg_normal_contaminated_mtr := rowSums(.SD), .SDcols = samples_normal_contaminated_mtr]
twins_filtered_dt[, agg_normal_contaminated_dep := rowSums(.SD), .SDcols = samples_normal_contaminated_dep]
twins_filtered_dt[, agg_normal_contaminated_vaf := agg_normal_contaminated_mtr / agg_normal_contaminated_dep]

######################################################################################################
# Analysis of QC-validated mutations 

# Plot the heatmap showing all mutations of interest 
mut_all_qc = as.matrix(twins_filtered_dt[, c(samples_vaf), with=FALSE])
rownames(mut_all_qc) = twins_filtered_dt[,1] %>% unlist()  
colnames(mut_all_qc) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()

col_annotation = data.frame('Tumour cell fraction' = c(0, 0.1, 0.2, 0, 0.9, 0.3, 0.6, 0.5, 0.8, 0.8, 0, 0, 0, 0, 0, 0.9, 0.9, 0.3, 0.7, 0.4, 0, 0.5), # fraction of tumour cells
                            Status = c(rep('normal', 4), rep('tumour', 6), rep('normal', 5), rep('tumour', 2), rep('normal', 1), c('tumour', 'normal', 'normal', 'tumour')),
                            Twin = c(rep('PD62341', 10), rep('PD63383', 8), rep('PD62341', 4)))
rownames(col_annotation) = colnames(mut_all_qc)
annotation_colors = list('Tumour.cell.fraction' = colorRampPalette(c('#f2e7e7', 'darkred'))(11),
                         Status = c(normal=col_normal, tumour=col_tumour),
                         Twin = c(PD62341=col_PD62341, PD63383=col_PD63383))
                        

# heatmap
pdf('Results/20241119_p4_heatmap_all_mutations_val258.pdf')
pheatmap(mut_all_qc,
         cellwidth=12, cellheight=0.4,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations, QC validated (258)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = F, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

pdf('Results/20241119_p4_heatmap_all_mutations_val258_large.pdf', height = 45)
pheatmap(mut_all_qc,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="All mutations, QC validated (259)", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=11, cexCol=2) 
dev.off()

######################################################################################################
# Trinucleotide context plots for mutations which are in the final 259 QC-validated set

# Function to identify the trinucleotide context of a given mutation, knowing its genomic cooordinates 
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

mybed_pass = twins_dt[mut_ID %in% muts, c('Chrom', 'Pos', 'Ref', 'Alt')]
trins_pass = get_trinucs(mybed_pass, BSgenome.Hsapiens.UCSC.hg38)
dt_pass = twins_dt[mut_ID %in% muts]
dt_pass$trins_pass=trins_pass

mut_sign_counts_pass = data.table(table(dt_pass[, trins_pass]))
setnames(mut_sign_counts_pass, c('V1', 'N'), c('trins', 'count'))
mut_sign_counts_pass[, mut_class := tstrsplit(trins, '.in', fixed=TRUE, keep=1)]
mut_sign_counts_pass[, mut_class := gsub("\\.", ">", mut_class)]
mut_sign_counts_pass[, context := tstrsplit(trins, '.', fixed=TRUE, keep=4)]

# aggregate by mutation class and context 
colors_sign = c('blue', 'black', 'red', 'grey', 'green', 'pink') # blue black red grey green pink 
ggplot(data=mut_sign_counts_pass, aes(x=context, y=count, fill=mut_class)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colors_sign)+
  facet_grid(~mut_class, scales = "free_x")+
  guides(fill="none")+ # remove legend
  ylim(c(0, 25))+
  labs(x = 'Context', y = 'Count', title = 'Final set of mutations (259)')+
  theme_classic(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme( strip.background = element_blank())+ 
  theme(panel.spacing = unit(0, "lines"))+
  theme(strip.text.x = element_text(size = 13))+
  geom_hline(yintercept = 0, colour="black", linewidth = 0.1)
ggsave('Results/20241119_p4_mut_trins_259muts_pass.pdf', width = 7.5, height = 3.5)

######################################################################################################
# Identify classes of mutations of interest 

# Mutations present in all normal samples
muts_all_normal = twins_filtered_dt[sum_normal_mtr_vaf>=9, mut_ID] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_all_normal, c('mut_ID', samples_vaf, 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf', 'vaf_all_tumour'), with=FALSE]
paste('Number of mutations present in all normal samples:', length(muts_all_normal)) # 5

# Mutations specific to PD62341 twin 
muts_PD62341_normal = twins_filtered_dt[sum_normal_PD63383_clean_mtr_vaf==0 & sum_normal_PD62341_mtr_vaf >= 4, mut_ID] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_PD62341_normal, c('mut_ID', 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf', 'vaf_all_tumour'), with=FALSE]
paste('Number of mutations specific to PD62341 twin:', length(muts_PD62341_normal)) # 7

# Mutations specific to PD63383 normal samples
muts_PD63383_normal = twins_filtered_dt[sum_normal_PD63383_mtr_vaf>=3 & sum_normal_PD62341_mtr_vaf <= 1, mut_ID] %>% unlist()
twins_filtered_dt[mut_ID %in% muts_PD63383_normal, c('mut_ID', 'agg_normal_PD62341_clean_vaf', 'agg_normal_PD63383_clean_vaf', 'vaf_all_tumour'), with=FALSE]
paste('Number of mutations specific to PD63383 twin:', length(muts_PD63383_normal)) # 11

# Mutations specific to the tumour 
muts_tumour = twins_filtered_dt[sum_tumour_mtr_vaf >=3 & sum_normal_contaminated_mtr_vaf <= 3 & sum_normal_mtr_vaf <= 3, mut_ID] %>% unlist()
paste('Number of mutations present in tumour samples but not normal samples:', length(muts_tumour)) # 87

# Mutations present in 2 PD62341 tumour samples 
muts_tumour_PD62341 = twins_filtered_dt[sum_tumour_PD62341_mtr_vaf == 2 & sum_normal_contaminated_mtr_vaf <= 1 & sum_tumour_mtr_vaf == 2, mut_ID] %>% unlist()
paste('Number of mutations specific to PD62341 tumour:', length(muts_tumour_PD62341)) # 11

# Mutations specific to PD63383 tumour samples
muts_tumour_PD63383 = twins_filtered_dt[sum_tumour_PD63383_mtr_vaf == 2 & sum_normal_contaminated_mtr_vaf <= 1 & sum_tumour_mtr_vaf == 2, mut_ID] %>% unlist()
paste('Number of mutations specific to PD63383 tumour:', length(muts_tumour_PD63383)) # 24

# Mutations present in 1 sample only
muts_one_sample_normal = twins_filtered_dt[sum_tumour_mtr_vaf==0 & sum_normal_mtr_vaf==1, mut_ID] %>% unlist() # 22
muts_one_sample_tumour = twins_filtered_dt[sum_tumour_mtr_vaf==1 & sum_normal_mtr_vaf==0, mut_ID] %>% unlist() # 88 
muts_one_sample = c(muts_one_sample_normal, muts_one_sample_tumour) 
paste('Number of mutations present in only 1 sample:', length(muts_one_sample)) # 110

######################################################################################################
# Plot distribution of VAF for mutations of interest 

twins_vaf_melt = melt(twins_filtered_dt[, c('mut_ID', samples_vaf), with=FALSE], id.vars = 'mut_ID')
twins_vaf_melt[, sample := tstrsplit(variable, '_VAF', fixed = TRUE, keep = 1)]
twins_vaf_melt[, sample := factor(sample, levels = c(
  'PD62341ad', 'PD62341aa', 'PD62341h', 'PD62341n', 'PD62341q', 'PD62341v',
  'PD63383ae', 'PD63383ak', 'PD63383bb', 'PD63383t', 'PD63383u', 'PD63383w',
  'PD62341ae', 'PD62341ag', 'PD62341aj', 'PD62341ak', 'PD62341am', 'PD62341ap', 'PD62341b', 'PD62341u',
  'PD63383ap', 'PD63383aq'))]
twins_vaf_melt[, twin := factor(fcase(sample %in% samples_PD62341, 'PD62341', sample %in% samples_PD63383, 'PD63383'))]
twins_vaf_melt[, status := factor(fcase(sample %in% samples_normal, 'normal', sample %in% samples_tumour, 'tumour'))]
twins_vaf_melt[, sample_type := factor(fcase(status == 'tumour', 'tumour', 
                                             status == 'normal' & twin == 'PD62341', 'PD62341 normal',
                                             status == 'normal' & twin == 'PD63383', 'PD63383 normal'))]

for (mut in muts_all_normal){
  
  ggplot(twins_vaf_melt[mut_ID==mut], aes(x = sample, y = value, col = sample_type))+
    geom_point(size=2.5)+
    theme_classic(base_size = 14)+
    labs(x = 'Sample', y = 'VAF', col = 'Sample')+
    ylim(c(0, 0.75))+
    ggtitle(glue('{mut}'))+
    scale_color_manual(values = c(col_PD62341, col_PD63383, col_tumour))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(glue('Results/20241114_p4_allNormal_dist_samples_{mut}.pdf'), height = 3, width = 6.5)
  
}


######################################################################################################
# Heatmap 
# create specific heatmaps for each class of mutation

# Mutations present in all normal samples
mut_all_normal = as.matrix(twins_filtered_dt[mut_ID %in% muts_all_normal, c(samples_vaf), with=FALSE])
rownames(mut_all_normal) = twins_filtered_dt[mut_ID %in% muts_all_normal,1] %>% unlist()  
colnames(mut_all_normal) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_allnormal4.pdf')
pheatmap(mut_all_normal,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations found in all normal samples", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Mutations specific to PD62341 twin 
mut_PD62341_normal = as.matrix(twins_filtered_dt[mut_ID %in% muts_PD62341_normal, c(samples_vaf), with=FALSE])
rownames(mut_PD62341_normal) = twins_filtered_dt[mut_ID %in% muts_PD62341_normal,1] %>% unlist()  
colnames(mut_PD62341_normal) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_allPD62341.pdf')
pheatmap(mut_PD62341_normal,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations specific to PD62341 twin", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Mutations specific to PD63383
mut_PD63383_normal = as.matrix(twins_filtered_dt[mut_ID %in% muts_PD63383_normal, c(samples_vaf), with=FALSE])
rownames(mut_PD63383_normal) = twins_filtered_dt[mut_ID %in% muts_PD63383_normal,1] %>% unlist()  
colnames(mut_PD63383_normal) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_allPD63383.pdf')
pheatmap(mut_PD63383_normal,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations specific to PD63383 twin", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()

# Mutations specific to the tumour 
mut_tumour = as.matrix(twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383), c(samples_vaf), with=FALSE])
rownames(mut_tumour) = twins_filtered_dt[mut_ID %in% c(muts_tumour, muts_tumour_PD62341, muts_tumour_PD63383),1] %>% unlist()  
colnames(mut_tumour) = tstrsplit(colnames(mut_all_qc), '_VAF', fixed=TRUE, keep=1) %>% unlist()
pdf('Results/20241119_p4_heatmap_muts_alltumour.pdf', height = 25)
pheatmap(mut_tumour,
         cellwidth=10, cellheight=10,
         annotation_col = col_annotation,
         annotation_colors = annotation_colors,
         main="Mutations specific to the tumour", 
         legend = T, 
         treeheight_row = 0,
         cluster_rows = T, cluster_cols = T, 
         show_rownames = T, show_colnames = T,
         fontsize=10, cexCol=2) 
dev.off()
