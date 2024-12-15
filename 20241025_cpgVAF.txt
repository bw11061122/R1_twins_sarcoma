#########################################################################################
# Notes from 2024-10-25 from my general README file on how the low-quality pileup was run 

#########################################################################################
2024-10-25
# ok trying to understand how to actually do this 

ssh farm22
cd /lustre/scratch126/casm/team274sb/bw18
mkdir -p twins_sarcoma/subs/pileup_20241025
cp ../hl11/twins_sarcoma/subs/subs_unmatched/twins_unm_subs.bed twins_sarcoma/subs/twins_unm_subs.bed
cd twins_sarcoma/subs
source /software/CASM/etc/module.CGP
module load cgpVAFcommand
3 (choose caveman_c)
createVafCmd.pl -pid 3434 -o pileup_20241025 -g /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -hdr /lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/HiDepth_v3.bed.gz -bo 1 -b twins_unm_subs.bed -mq 5 -bq 5

edit run_bsub.sh to increase memory requirements (because there are 22 samples and failed on standard)

#!/usr/bin/env bash
set -e;
set -o pipefail;
JID=`uuidgen`
bsub -oo cgpvaf%I.log -G team274-grp -q normal -J "$JID[1-96]%5" -n 1 -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' -M2500 '/software/CASM/modules/installs/cgpVAFcommand/cgpVAFcommand-2.5.0/perl/bin/farm_idx_exec.pl cgpVafChr.cmd $LSB_JOBINDEX'
bsub -w "done($JID)||exit($JID)" -oo concat%I.log -q normal -G team274-grp -J 'batchjobs[1-4]%4' -n 1 -R 'select[mem>=1000] span[hosts=1] rusage[mem=1000]' -M1000 '/software/CASM/modules/installs/cgpVAFcommand/cgpVAFcommand-2.5.0/perl/bin/farm_idx_exec.pl cgpVafConcat.cmd $LSB_JOBINDEX'

bash run_bsub.sh

okay got this thing done!

cd /lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/pileup_20241025/output/PD62341v/snps
# get the tsv files onto your laptop 

rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/pileup_20241025/output/PD62341v/snp/*.tsv" /Users/bw18/Desktop/1SB/Data/
rsync -avhu bw18@farm22:"/lustre/scratch126/casm/team274sb/bw18/twins_sarcoma/subs/pileup_20241025/output/PD63383w/snp/*.tsv" /Users/bw18/Desktop/1SB/Data/

# on your local: cd Desktop/1SB/Data/pileup_20241025
bash ../../Scripts/R1/20241011_removeheader.sh 