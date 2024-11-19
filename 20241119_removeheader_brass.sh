# 2024-11-19
# BRASS rearrangement calls 
# remove headers from each file 

# cd ~/Desktop/1SB/Data

# remove the header 
for file in *.bedpe
do 
    base=`basename $file .bedpe`
    echo $base
    out=${base}_removed_header.bedpe 
    grep -v '#' $file > $out 
done 