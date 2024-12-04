####################################################################################################################################
# SCRIPT 10

# Get bed file of 255 mutations used to reconstruct the final phylogeny 
# Barbara Walkowiak bw18
# 2024-12-04

# INPUT: 
# 1 list of 255 mutation IDs (coordinate, ref, alt) which are used to reconstruct the final phylogeny

# OUTPUT:
# bed file with coordinates of each mutation (used to do genotyping via alleleIntegrator)

####################################################################################################################################
# Read in the list of mutations

muts_dt = data.table(read.csv('Data/20241114_599muts_QCJbrowse.csv', header = T))
muts = data.table(muts_dt[Jbrowse.quality == 'Y', mut_ID])
setnames(muts, 'V1', 'mut_ID')

# create the desired bed file 
muts[, Chrom := tstrsplit(mut_ID, '_', fixed=T, keep=1)] # chromosome 
muts[, start := tstrsplit(mut_ID, '_', fixed=T, keep=2)] # position
muts[, start := as.numeric(start)]
muts[, end := as.numeric(start)]
muts[, width := 1]
muts[, strand := "*"]
muts[, Ref := tstrsplit(mut_ID, '_', fixed=T, keep=3)]
muts[, Alt := tstrsplit(mut_ID, '_', fixed=T, keep=4)]
muts[, varID := paste(sep="_", Chrom, start, end)]
muts[, Chrom := gsub("chr","", Chrom)]
muts = muts[, mut_ID := NULL]# drop mut_ID column
            
setnames(muts, c('Chrom', 'start', 'end', 'width', 'strand', 'Ref', 'Alt'), c("seqnames","start","end","width","strand","Ref","Alt"))
            
write.table(muts, 'Data/20241204_mutations255.bed', sep = '\t', quote=F, row.names=F)

library(GenomicFeatures)
loci = GRanges(muts$seqnames,IRanges(muts$start, muts$end),
               varID = muts$varID)
mcols(loci) = cbind(mcols(loci),muts[match(loci$varID,muts$varID)])


