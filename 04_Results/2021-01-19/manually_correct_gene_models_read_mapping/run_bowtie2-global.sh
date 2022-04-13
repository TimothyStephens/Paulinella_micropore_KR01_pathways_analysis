#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH=$PATH:"/home/timothy/programs/bowtie2-2.3.5.1-linux-x86_64"
export PATH=$PATH:"/home/timothy/programs/samtools-1.11/bin"
NCPUS=48
REF="Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt"


#### Start Script
while read LINE;
do
        R1=$(echo $LINE | awk '{print $3}')
        R2=$(echo $LINE | awk '{print $4}')
        OUT=$(echo $LINE | awk '{print $2}')
	run_cmd "bowtie2 --very-sensitive --no-unal -x $REF -1 $R1 -2 $R2 --threads $NCPUS 2>$OUT.bowtie2_global_mapping.stats | samtools sort -@ $NCPUS - > $OUT.global_mapping.coordsorted.bam"
done < samples.txt


