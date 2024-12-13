#! /usr/bin/env Rscript

####################################################################################################################################
# SCRIPT 8

# scRNA-seq genotyping 
# 2024-12-04
# Barbara Walkowiak 

###################################################################################################################################
# LIBRARIES
# the way to do this now is to use the conda environment that has been set up on the farm
# following the instructions from Behjati Lab GitHub

# Load libraries and necessary functions to run alleleIntegrator 
library(alleleIntegrator)
library(GenomicFeatures)
library(tidyverse)
library(data.table)
source('/lustre/scratch125/casm/team274sb/mt22/alleleIntegrator_BehjatiLab/R/alleleIntegrator_helperFunctions.R')
#source('/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/sc_utils.R')
#source('/lustre/scratch125/casm/team274sb/mt22/generalScripts/utils/misc.R')

# create output directory and set working directory to output 
outDir = '/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scRNA_genotyping/output'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}
setwd(outDir)

###################################################################################################################################
# Import mutations of interest (255 mutations used to reconstruct the phylogeny we consider informative)

setwd("/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/scRNA_genotyping")
phylo_mutations = fread('20241204_mutations255.bed', sep = '\t', header = TRUE) # import bed file with mutation calls 

# Convert to Granges object
loci = GRanges(phylo_mutations$seqnames,IRanges(phylo_mutations$start, phylo_mutations$end),
               varID = phylo_mutations$varID)
mcols(loci) = cbind(mcols(loci),phylo_mutations[match(loci$varID,phylo_mutations$varID)])

# Index of the sample
gosh_ID = c("GOSH100", "GOSH101") # GOSH IDs 
sanger_ID = c("PD62341", "PD63383") # Sanger / CASM IDs

tumourDNA = paste(sep="", '/nfs/cancer_ref01/nst_links/live/3434/',sanger_ID,'/',sanger_ID,'.sample.dupmarked.bam')
names(tumourDNA) = basename(dirname(tumourDNA))

###################################################################################################################################
# Run alleleCounter; aggregates base count at the specified positions 

# Get list of scRNAseq bam files of interest 
# GOSH IDs of interest: GOSH100, GOSH101
bams10X = list.files('/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_bams', full.names = T,recursive = T)
bams10X = bams10X[grepl('bam$', bams10X) & grepl('GOSH100|GOSH101', bams10X)]
print(paste('Bam files invesitgated:', bams10X))
namesvec = sapply(1:length(bams10X), function(x) strsplit(strsplit(bams10X[x], fixed=T, split="_")[[1]][3], fixed=T, split="/")[[1]][3])
names(bams10X) = basename(namesvec)
bamFiles = bams10X
names(bamFiles) = namesvec

# specify parameter values 
# path to the reference genome the same as provided in Nathan's script 
refGenome = '/nfs/srpipe_references/downloaded_from_10X/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
nParallel = 30 # parallelization? not 100% sure what this emans 
outputs = file.path(outDir,paste0(names(bamFiles),'_somaticVar_scRNA_alleleCounts.tsv')) # create the path to the output file
params = list(bams = bamFiles, refGenome = refGenome, tgtLoci = loci, outputs = outputs, nParallel = nParallel)

# call alleleCounter using parameters specified 
cnts = do.call(alleleCounter, params)

print('AlleleCounter completed!')
