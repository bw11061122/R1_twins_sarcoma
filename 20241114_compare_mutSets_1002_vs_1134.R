# 2024-11-14 quick script to compare my mutation sets and track down why now I have more 

muts_1 = read.table('Data/mutations_include_20241106_1002.txt') %>% unlist()
muts_2 = read.table('Data/mutations_include_20241113_1134.txt') %>% unlist()

setdiff(muts_1, muts_2) # nothing was excluded 
setdiff(muts_2, muts_1) # 132

# mutations which are now excluded but were included in the previous dataframe 
twins_dt[mut_ID %in% setdiff(muts_1, muts_2), c('mut_ID', columns_req_filters), with=FALSE] # 0

# mutations which are now included but were excluded in the previous dataframe 
twins_dt[mut_ID %in% setdiff(muts_2, muts_1), c('mut_ID', samples_vaf), with=FALSE]

# load the previous dataframe with filters - what filter was failed previously?
dt = data.table(read.csv('Data/twins_dt_filters_20241112.csv'))
filter_cols = grep('^f', names(dt), value = TRUE)
dt[mut_ID %in% setdiff(muts_2, muts_1), c('mut_ID', filter_cols, 'sum_req_filters_qual'), with=FALSE]

# all of those mutations present in one sample only 
# I don't intend to remove them so will keep the extra 132 (won't change the phylogeny since present in 1 sample)


