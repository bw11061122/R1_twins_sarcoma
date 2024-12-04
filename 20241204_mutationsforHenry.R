
setwd('/Users/bw18/Desktop/1SB')
twins_dt = data.table(read.csv('Data/pileup_merged_20241016.tsv'))
mut_all = data.table(twins_dt[, mut_ID])
# add mutation category 
mut_all[, cat := factor(fcase(
  V1 %in% muts_tumour_spec, 'tumour mutations',
  !V1 %in% muts_tumour_spec, 'other'
))]

write.csv(mut_all, 'Data/20241204_all_mutation_calls.csv')

# okay let's try doing this again
twins_dt_info[, in_final_255 := factor(fcase(
  mut_ID %in% muts, 1,
  !mut_ID %in% muts, 0
))]

twins_dt_info[, tumour_specific := factor(fcase(
  mut_ID %in% muts_tumour_spec, 1,
  !mut_ID %in% muts_tumour_spec, 0
))]

# save to a file
write.csv(twins_dt_info, 'Data/20241204_muts_with_qc.csv', quote=F)

# select mutation and filter columns only
twins_dt_info[, in_final := as.numeric(sum_req_filters==0 & in_final_255==1)]
