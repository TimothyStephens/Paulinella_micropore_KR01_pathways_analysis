#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu



##
## Build tree.
##
build_trees() {
	IQTREE="/home/timothy/programs/iqtree-1.6.12-Linux/bin"
	NCPUS=4
	FILE="${1}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.trimal1"
	NSEQS=$(grep -c '>' "$FILE")
	if [ $NSEQS -gt 3 ]; then
		run_cmd "$IQTREE/iqtree -s $FILE -nt AUTO -ntmax $NCPUS -m MFP -bb 2000 -nm 2000 -quiet"
	else
		echo "   - Too few seqs! (only $NSEQS seqs)"
	fi
}

export -f build_trees
parallel -j 24 build_trees :::: <(awk -F'\t' '$1~"MSTRG" {print $2}' SuppTable_S2.transcript_info.txt | sort | uniq)

echo ""; echo "Done running iqtree!"

