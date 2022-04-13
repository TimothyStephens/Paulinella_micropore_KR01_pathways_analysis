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
	NCPUS=6
	FILE="${1}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa"
	exec 1> "${FILE}.build_trees.log" 2>&1
	echo "## $FILE"
	NSEQS=$(grep -c '>' "$FILE")
	if [ $NSEQS -gt 3 ]; then
		run_cmd "$IQTREE/iqtree -s $FILE.mafft.aln -nt AUTO -ntmax $NCPUS -m LG+R7 -bb 2000 -nm 2000 -quiet"
	else
		echo "   - Too few seqs! (only $NSEQS seqs)"
	fi
}

export -f build_trees
parallel -j 20 build_trees :::: <(awk -F'\t' '$1~"MSTRG" {print $2}' SuppTable_S2.transcript_info.txt | sort | uniq)

echo ""; echo "Done running iqtree!"

