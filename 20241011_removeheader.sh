#!/bin/bash

# Barbara Walkowiak bw18 / bw450
# written 2024/10/11, re-sed 2024/10/16

# determine what the header actually is and save to a file
grep '##' PD62341v_PD62341q_snp_vaf.tsv > 20241011_pileup_header.txt

# remove the header 
for file in *_vaf.tsv 
do 
    base=`basename $file .tsv`
    echo $base
    out=${base}_removed_header.tsv 
    grep -v '##' $file > $out 
done 
