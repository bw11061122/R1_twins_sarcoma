#!/bin/bash
#BSUB -o /nfs/users/nfs_b/bw18/logs/telomerecat.out
#BSUB -e /nfs/users/nfs_b/bw18/logs/telomerecat.err
#BSUB -n 4
#BSUB -M 1000
#BSUB -R 'select[mem>1000] rusage[mem=1000]'
#BSUB -q long
#BSUB -G team274-grp
 
SCRIPT_PATH=/nfs/users/nfs_b/bw18/scripts/20241113_telomere_length.sh
 
bash $SCRIPT_PATH

# to run:
# bsub < scripts/20241113_telomere_length_bsub.sh