#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH=$PATH:"/home/timothy/programs/samtools-1.11/bin"
NCPUS=6
REF="Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt"


#### Start Script
for BAM in *.local_mapping.coordsorted.bam;
do
	run_cmd "samtools index -@ $NCPUS $BAM"
done

run_cmd "samtools merge -@ $NCPUS -r All_combined.local_mapping.coordsorted.bam CL*.local_mapping.coordsorted.bam HL*.local_mapping.coordsorted.bam"
run_cmd "samtools index -@ $NCPUS All_combined.local_mapping.coordsorted.bam"



