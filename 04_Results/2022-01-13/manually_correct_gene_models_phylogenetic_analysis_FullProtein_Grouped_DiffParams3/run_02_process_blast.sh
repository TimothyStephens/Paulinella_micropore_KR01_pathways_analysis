#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"

#### Start Script

##
## Get just the seqs from the (large) DB so its quicker to extract hit regions later on.
##
REF="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
DB="/scratch/timothy/databases/db_from_Dana"
DB1="$DB/XX-BACTERIA-FINAL.faa"
DB2="$DB/XX-METAZOA-FINAL.faa"
DB3="$DB/XX-MMETSP-FINAL.faa"
DB4="$DB/XX-OTHER-FINAL.faa"

cat "$DB1" "$DB2" "$DB3" "$DB4" | ./grepf_fasta.py -s <(cut -f2 *".blastp_"{BACTERIA,METAZOA,MMETSP,OTHER}".outfmt6" | sort | uniq) \
 | seqkit fx2tab | sort | uniq | seqkit tab2fx > "BLAST_DB.all_subjects.faa"
cat $REF >> "BLAST_DB.all_subjects.faa"



##
## Combine BLASTP hits from each DB then down sample taxa. Also get coords of hits for later.
##
IDS="SuppTable_S2.transcript_info.txt"
while read KO; do
	echo "## $KO"
	
	rm -fr "$KO.blastp_ALL.outfmt6.parsed.filtered.faa"
	while read PEP; do
		# Filter blast output. Filter by e-value < 1e-10 (e-values upto 1000 was allowed during blast analysis).
		# Sort by query_id then bitscore (ignore cases where bitscore is the same [cant use e-value to break as we use different databases of different sizes])
		cat "$PEP.pep.faa.blastp_"{BACTERIA,METAZOA,MMETSP,OTHER}".outfmt6" | grep -v 'Eukaryota-Cercozoa-NA-Euglyphida-Paulinellidae-Paulinella' | awk -F'\t' '($11 < 1e-10)' | sort -k1,1 -k12,12nr > "$PEP.blastp_ALL.outfmt6"
		./down_sample_taxa.sh "$PEP.blastp_ALL.outfmt6"
		./grepf_column.py -c 2 -f "$PEP.blastp_ALL.outfmt6.parsed" -i "$PEP.blastp_ALL.outfmt6" -o "$PEP.blastp_ALL.outfmt6.parsed.filtered"
		
		# Get BLASTP subject + query seqs (need for building trees)
		./grepf_fasta.py -i "BLAST_DB.all_subjects.faa" -s <(awk -F'\t' '{print $1"\n"$2}' "$PEP.blastp_ALL.outfmt6.parsed.filtered" | sort | uniq) -o "$PEP.blastp_ALL.outfmt6.parsed.filtered.faa"
		
		cat "$PEP.blastp_ALL.outfmt6.parsed.filtered.faa" >> "$KO.blastp_ALL.outfmt6.parsed.filtered.faa"
	done < <(awk -F'\t' -vKO="${KO}" '$2==KO && $1~"MSTRG" {print $1}' "$IDS" | sort | uniq)
	
	cat "$KO.Other_Paulinella_seqs.pep.faa" >> "$KO.blastp_ALL.outfmt6.parsed.filtered.faa"
	
	cat "$KO.blastp_ALL.outfmt6.parsed.filtered.faa" | seqkit fx2tab | sort | uniq | seqkit tab2fx > "$KO.blastp_ALL.outfmt6.parsed.filtered.uniq.faa"
done < <(awk -F'\t' '$1~"MSTRG" {print $2}' "$IDS" | sort | uniq)


