#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu



##
## Align and build tree.
##
build_trees() {
	MAFFT="/home/timothy/programs/mafft-7.453-with-extensions/bin"
	IQTREE="/home/timothy/programs/iqtree-1.6.12-Linux/bin"
	NCPUS=2
	FILE="${1}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa"
	exec 1> "${FILE}.build_trees.log" 2>&1
	echo "## $FILE"
	NSEQS=$(grep -c '>' "$FILE")
	if [ $NSEQS -gt 3 ]; then
		run_cmd "$MAFFT/mafft --thread $NCPUS --localpair --maxiterate 1000 $FILE > $FILE.mafft.aln 2> $FILE.mafft.log"
		run_cmd "$IQTREE/iqtree -s $FILE.mafft.aln -nt AUTO -ntmax $NCPUS -m MFP -bb 2000 -quiet"
	else
		echo "   - Too few seqs! (only $NSEQS seqs)"
	fi
}

export -f build_trees
parallel -j 36 build_trees :::: <(awk -F'\t' '$1~"MSTRG" {print $2}' SuppTable_S2.transcript_info.txt | sort | uniq)

echo ""; echo "Done running iqtree!"

