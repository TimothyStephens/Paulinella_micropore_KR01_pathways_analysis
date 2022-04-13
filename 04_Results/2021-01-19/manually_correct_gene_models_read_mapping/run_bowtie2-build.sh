#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH=$PATH:"/home/timothy/programs/bowtie2-2.3.5.1-linux-x86_64"
REF="Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt"

#### Start Script
run_cmd "bowtie2-build $REF $REF"



