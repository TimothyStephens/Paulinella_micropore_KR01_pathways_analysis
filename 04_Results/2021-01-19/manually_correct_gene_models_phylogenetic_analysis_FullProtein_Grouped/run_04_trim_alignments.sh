#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu



##
## Trim alignment
##
TRIMAL="/home/timothy/programs/trimal-1.4.1/bin"

IDS="SuppTable_S2.transcript_info.txt"
while read KO; do
	FILE="${KO}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa"
	echo "## $FILE"
       	NSEQS=$(grep -c '>' "$FILE")
       	if [ $NSEQS -gt 3 ]; then
		run_cmd "$TRIMAL/trimal -automated1 -in $FILE.mafft.aln -out $FILE.mafft.aln.trimal1"
       	else
       	       	echo "   - Too few seqs! (only $NSEQS seqs)"
       	fi
done < <(awk -F'\t' '$1~"MSTRG" {print $2}' "$IDS" | sort | uniq)

echo ""; echo "Done trimming alignments!"

