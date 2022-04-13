#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu



#### Start Script
run_BLAST() {
	BLAST="/home/timothy/programs/ncbi-blast-2.10.1+/bin"
	NCPUS=1
	PEP="${1}"
	exec 1> "${PEP}.blast.log" 2>&1
	echo "## $PEP"
	
	DB="/scratch/timothy/databases/db_from_Dana/XX-BACTERIA-FINAL.faa"
	OUT="$PEP.blastp_BACTERIA.outfmt6"
	run_cmd "$BLAST/blastp -query $PEP -db $DB -out $OUT -num_threads $NCPUS -max_target_seqs 2000 -evalue 1000 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
	
	DB="/scratch/timothy/databases/db_from_Dana/XX-METAZOA-FINAL.faa"
	OUT="$PEP.blastp_METAZOA.outfmt6"
	run_cmd "$BLAST/blastp -query $PEP -db $DB -out $OUT -num_threads $NCPUS -max_target_seqs 2000 -evalue 1000 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
	
	DB="/scratch/timothy/databases/db_from_Dana/XX-MMETSP-FINAL.faa"
	OUT="$PEP.blastp_MMETSP.outfmt6"
	run_cmd "$BLAST/blastp -query $PEP -db $DB -out $OUT -num_threads $NCPUS -max_target_seqs 2000 -evalue 1000 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
	
	DB="/scratch/timothy/databases/db_from_Dana/XX-OTHER-FINAL.faa"
	OUT="$PEP.blastp_OTHER.outfmt6"
	run_cmd "$BLAST/blastp -query $PEP -db $DB -out $OUT -num_threads $NCPUS -max_target_seqs 2000 -evalue 1000 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""
}

export -f run_BLAST
parallel -j 48 run_BLAST :::: "files2run.txt"

echo ""; echo "Done running BLAST!"

