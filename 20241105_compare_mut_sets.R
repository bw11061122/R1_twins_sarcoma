# Compare lists of mutations

muts_771 = read.table('Data/mutations_include_20241106_771.txt') %>% unlist()
muts_312 = read.table('Data/mutations_include_20241106_312.txt') %>% unlist()
muts_380 = read.table('Data/mutations_include_20241106_380.txt') %>% unlist()

setdiff(muts_312, muts_771)

# "chr16_22541122_G_A" 
# "chr1_111337001_A_G" 
# "chr1_16932907_G_A"  
# "chr7_152850206_G_A"

sum(muts_tumour_specific_69 %in% muts_312)
sum(muts_tumour_specific_18 %in% muts_312)
sum(muts_tumour_specific_69 %in% muts_771)
sum(muts_tumour_specific_18 %in% muts_771)
sum(muts_tumour_specific_69 %in% muts_771)
sum(muts_tumour_specific_18 %in% muts_771)
# mutations of interest present in all those sets

# another one sorry!

muts_1 = read.table('Data/mutations_include_20241106_1002.txt') %>% unlist()
muts_2 = read.table('Data/mutations_include_20241113_1134.txt') %>% unlist()

setdiff(muts_1, muts_2) # nothing was excluded 
setdiff(muts_2, muts_1) # 132

# mutations which are now excluded but were included in the previous dataframe 
twins_dt[mut_ID %in% setdiff(muts_1, muts_2), c('mut_ID', columns_req_filters), with=FALSE] # 0

# mutations which are now included but were excluded in the previous dataframe 
twins_dt[mut_ID %in% setdiff(muts_2, muts_1), c('mut_ID', samples_vaf), with=FALSE]


