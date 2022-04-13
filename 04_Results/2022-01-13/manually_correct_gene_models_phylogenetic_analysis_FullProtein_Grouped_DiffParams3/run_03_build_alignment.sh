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
	NCPUS=8
	FILE="${1}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa"
	exec 1> "${FILE}.build_alignments.log" 2>&1
	echo "## $FILE"
	NSEQS=$(grep -c '>' "$FILE")
	if [ $NSEQS -gt 3 ]; then
		run_cmd "$MAFFT/mafft --thread $NCPUS --retree 2 --maxiterate 0 $FILE > $FILE.mafft.aln 2> $FILE.mafft.log"
	else
		echo "   - Too few seqs! (only $NSEQS seqs)"
	fi
}

export -f build_trees
parallel -j 10 build_trees :::: <(awk -F'\t' '$1~"MSTRG" {print $2}' SuppTable_S2.transcript_info.txt | sort | uniq)

echo ""; echo "Done running iqtree!"

