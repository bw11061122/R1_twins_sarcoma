####################################################################################################################################
# SCRIPT 9

# scRNA-seq genotyping 
# 2024-12-04
# Barbara Walkowiak 

###################################################################################################################################
# LIBRARIES
# the way to do this now is to use the conda environment that has been set up on the farm
# following the instructions from Behjati Lab GitHub
# ssh <user>@farm22
# source /software/cellgen/team274/miniconda3/etc/profile.d/conda.sh
# which conda
# conda activate alleleIntegrator

# Load libraries and necessary functions to run alleleIntegrator 
library(alleleIntegrator)
library(GenomicFeatures)
library(tidyverse)
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
phylo_mutations = fread('20241204_mutations255.bed', sep = '\t') # import bed file with mutation calls 
phyloMut = phylo_mutations[,c('chrom', 'pos', 'pos'), with=FALSE]
phyloMut[, width := 1]
phyloMut[, strand := "*"]
phyloMut[, Ref := phylo_mutations[, Ref]]
phyloMut[, Alt := phylo_mutations[, Alt]]
setnames(phyloMut, c('chrom', 'pos', 'pos', 'width', 'strand', 'Ref', 'Alt'), c("seqnames","start","end","width","strand","Ref","Alt"))
phyloMut[, varID := paste(sep="_", seqnames, start, end)]
phyloMut[, seqnames := gsub("chr","", seqnames)]

# Convert to Granges object
loci = GRanges(somaticMut$seqnames,IRanges(somaticMut$start, somaticMut$end),
               varID = somaticMut$varID)
mcols(loci) = cbind(mcols(loci),somaticMut[match(loci$varID,somaticMut$varID),!colnames(somaticMut) %in% c('X')])

# Index of the sample
gosh_ID = c("GOSH100", "GOSH101") # GOSH IDs 
sanger_ID = c("PD63383", "PD62341") # Sanger / CASM IDs

tumourDNA = paste(sep="", '/nfs/cancer_ref01/nst_links/live/3434/',sanger_ID,'/',sanger_ID,'.sample.dupmarked.bam')
names(tumourDNA) = basename(dirname(tumourDNA))

###################################################################################################################################
# Run alleleCounter; aggregates base count at the specified positions 

# Get list of scRNAseq bam files of interest 
# GOSH IDs of interest: GOSH100, GOSH101
bams10X = list.files('/lustre/scratch126/casm/team274sb/project_folders/Sarcoma/sc_bams', full.names = T,recursive = T)
bams10X = bams10X[grepl('bam$', bams10X) & grepl('GOSH100|GOSH101', bams10X)]
print(paste('Bam files invesitgated:', bam10X))
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
