
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))
mut_all = data.table(twins_dt[, mut_ID])
# add mutation category 
mut_all[, cat := factor(fcase(
  V1 %in% muts_tumour_spec, 'tumour mutations',
  !V1 %in% muts_tumour_spec, 'other'
))]

write.csv(mut_all, 'Data/20241204_all_mutation_calls.csv')