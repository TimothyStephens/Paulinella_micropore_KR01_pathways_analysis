#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate py27; set -eu

export PATH=$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin
export PATH=$PATH:/home/timothy/programs/hmmer-3.1b2/bin
TRANSDECODER="/home/timothy/programs/TransDecoder-v5.5.0"
BLAST_TOP_HIT="/home/timothy/scripts/blast_top_hits.py"

NCPUS=60

SWISSPROT="/scratch/timothy/databases/UniProt_2020_05/uniprot_sprot.fasta"
TREMBL="/scratch/timothy/databases/UniProt_2020_05/uniprot_trembl.fasta"
PFAM_A="/scratch/timothy/databases/PFAM_33.1/Pfam-A.hmm"

GENOME="Paulinella_micropora_KR01_nuclear.assembly.fasta"



#### Start Script
GTF="combined.manually_fixed.gtf"
OUT="${GTF%*.gtf}"

##
run_cmd "$TRANSDECODER/util/gtf_genome_to_cdna_fasta.pl $GTF $GENOME > $OUT.transcripts.fasta"
run_cmd "$TRANSDECODER/util/gtf_to_alignment_gff3.pl $GTF > $OUT.transcripts.gff3"

##
run_cmd "$TRANSDECODER/TransDecoder.LongOrfs -S -m 30 -t $OUT.transcripts.fasta > $OUT.transcripts.TransDecoder.LongOrfs.log 2>&1"

##
run_cmd "blastp -query $OUT.transcripts.fasta.transdecoder_dir/longest_orfs.pep  \
    -db $SWISSPROT  -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads $NCPUS \
    -out $OUT.transcripts.fasta.transdecoder_dir/longest_orfs.pep.blastp_SwissProt.outfmt6"

run_cmd "hmmscan --cpu $NCPUS \
	--domtblout $OUT.transcripts.fasta.transdecoder_dir/longest_orfs.pep.pfam.domtblout \
	$PFAM_A \
	$OUT.transcripts.fasta.transdecoder_dir/longest_orfs.pep \
	> $OUT.transcripts.fasta.transdecoder_dir/longest_orfs.pep.hmmscan.log 2>&1"

##
run_cmd "$TRANSDECODER/TransDecoder.Predict --single_best_only \
	-t $OUT.transcripts.fasta \
	--retain_blastp_hits $OUT.transcripts.fasta.transdecoder_dir/longest_orfs.pep.blastp_SwissProt.outfmt6 \
	--retain_pfam_hits $OUT.transcripts.fasta.transdecoder_dir/longest_orfs.pep.pfam.domtblout \
	> $OUT.transcripts.TransDecoder.Predict.log 2>&1"

##
run_cmd "$TRANSDECODER/util/cdna_alignment_orf_to_genome_orf.pl \
	$OUT.transcripts.fasta.transdecoder.gff3 \
	$OUT.transcripts.gff3 \
	$OUT.transcripts.fasta > $OUT.transcripts.fasta.transdecoder.genome.gff3"

## 
run_cmd "$TRANSDECODER/util/force_start_codon_refinement.pl --gff3_file $OUT.transcripts.fasta.transdecoder.gff3 --transcripts $OUT.transcripts.fasta > $OUT.transcripts.fasta.transdecoder.forceStart.gff3"
run_cmd "$TRANSDECODER/util/gff3_file_to_proteins.pl --gff3 $OUT.transcripts.fasta.transdecoder.forceStart.gff3 --fasta $OUT.transcripts.fasta --seqType prot > $OUT.transcripts.fasta.transdecoder.forceStart.pep"
run_cmd "$TRANSDECODER/util/gff3_file_to_proteins.pl --gff3 $OUT.transcripts.fasta.transdecoder.forceStart.gff3 --fasta $OUT.transcripts.fasta --seqType CDS  > $OUT.transcripts.fasta.transdecoder.forceStart.cds"
run_cmd "$TRANSDECODER/util/cdna_alignment_orf_to_genome_orf.pl $OUT.transcripts.fasta.transdecoder.forceStart.gff3 $OUT.transcripts.gff3 $OUT.transcripts.fasta > $OUT.transcripts.fasta.transdecoder.forceStart.genome.gff3"



