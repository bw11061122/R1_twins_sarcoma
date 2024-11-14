module load telomerecat

for bam in /nfs/cancer_ref01/nst_links/live/3434/PD*/*telbam.bam

do
echo $bam
samp="$(echo $bam | cut -d'/' -f7)"
echo $samp
cp $bam .
telomerecat telbam2length ${samp}.telbam.bam --output ~/tel/${samp}_telomere_length.csv
rm ${samp}.telbam.bam

done
