#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu


export PATH="$PATH:/home/timothy/programs/salmon-1.1.0_linux_x86_64/bin"
NCPUS=24

TRANSCRIPTS="Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt"
KMER=31
TYPE="puff"
INDEX="$TRANSCRIPTS.salmon_$TYPE.$KMER.idx"

#### Start Script
run_cmd "salmon index --transcripts $TRANSCRIPTS --index $INDEX --kmerLen $KMER --threads $NCPUS"


