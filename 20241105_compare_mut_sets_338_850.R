# Compare lists of mutations

muts_850 = read.table('Data/mutations_include_20241105_850.txt') %>% unlist()
muts_338 = read.table('Data/mutations_include_20241104_338.txt') %>% unlist()

setdiff(muts_338, muts_850)

# "chr10_49963606_G_A" # strand bias 
# "chr11_18922568_G_T" # strand bias 
# "chr13_18565005_G_A" # strand bias
# "chr14_73648334_G_A" # strand bias
# "chr16_22541122_G_A" # germline LQ
# "chr1_111337001_A_G" # germline LQ
# "chr1_145588334_T_C" # strand bias 
# "chr1_148153535_T_C" # strand bias
# "chr1_16932907_G_A" # germline LQ 
# "chr1_248574696_C_G" # strand bias
# "chr1_25404075_T_C" # strand bias 
# "chr1_83136026_A_T" # strand bias
# "chr1_83387756_A_G"  # strand bias
# "chr2_113525233_T_C" # strand bias
# "chr2_95747690_G_A"  # strand bias
# "chr5_795141_C_T" # strand bias  
# "chr6_167319952_G_A" # strand bias
# "chr6_93635834_A_G"  # strand bias
# "chr7_152850206_G_A" # germline LQ
# "chr7_57104457_T_A" # strand bias
# "chr8_86825503_C_T" # strand bias
# "chrX_47326908_G_A" # strand bias

length(setdiff(muts_850, muts_338))
setdiff(muts_850, muts_338)

sum(muts_tumour_specific_70 %in% muts_338)
sum(muts_tumour_specific_18 %in% muts_338)

muts_771 = read.table('Data/mutations_include_20241106_771.txt') %>% unlist()
muts_312 = read.table('Data/mutations_include_20241106_312.txt') %>% unlist()
sum(muts_tumour_specific_70 %in% muts_312)
sum(muts_tumour_specific_18 %in% muts_312)
sum(muts_tumour_specific_70 %in% muts_771)
sum(muts_tumour_specific_18 %in% muts_771)

muts_tumour_specific_70 %in% setdiff(muts_338, muts_312) %>% unlist()
muts_tumour_specific_70[28]


