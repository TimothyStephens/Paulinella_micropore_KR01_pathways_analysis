# Visualize manually corrected genes using IGV

Link to files with info about the manually corrected genes. This info will be visualized using `IGV`.

## Setup analysis directory

Link data files with manually corrected genes.

```bash
F="manually_correct_gene_models"
ln -s ../$F/combined.manually_fixed.transcripts.fasta
ln -s ../$F/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.genome.gff3
ln -s ../$F/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.gff3
ln -s ../$F/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep

F="manually_correct_gene_models_blastp_UniProt_KEGG_seqs"
ln -s ../$F/combined.manually_fixed.blastp_KEGG_UniProt_tophit.CDS_coords.bed
ln -s ../$F/combined.manually_fixed.blastp_KEGG_UniProt_tophit_fullAlnInfo.PEP_coords.bed

F="manually_correct_gene_models_read_mapping"
ln -s ../$F/combined.manually_fixed.transcripts.fasta.All_combined.global_mapping.coordsorted.bam
ln -s ../$F/combined.manually_fixed.transcripts.fasta.All_combined.global_mapping.coordsorted.bam.bai
ln -s ../$F/combined.manually_fixed.transcripts.fasta.All_combined.local_mapping.coordsorted.bam
ln -s ../$F/combined.manually_fixed.transcripts.fasta.All_combined.local_mapping.coordsorted.bam.bai

```

Link to custom script for:

-  Parsing text files
- Extracting coords along CDS where exon/intron boundaries are
- Find coords of alternative start codons

```bash
ln -s ../../../02_Scripts/add_value_to_table.py
ln -s ../../../02_Scripts/splice_site_positions_along_CDS.py
ln -s ../../../02_Scripts/find_alternate_start_sites_in_a_protein.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"
```

## Reformat data ready for visualization

```bash
FASTA="combined.manually_fixed.transcripts.fasta"
GFF3="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.gff3"
PEP="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
PFAM="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep.pfam_scan.out"
```

Get offset of ORF in transcript.

```bash
awk -F'\t' '$3=="CDS" {print $9"\t"($4-1)"\t"$1}' $GFF3 | sed -e 's/.*Parent=//' > $GFF3.ORF_offset.txt
```

Adjust domains from `pfam_scan.pl` so they are correctly positioned along the original transcript (i.e. need to add length of 5-prime UTR to coords to adjust for ORF position).

Run `pfam_scan.pl`.

```bash
./run_pfamscan.sh
```

Check log file `run_pfamscan.sh.log.*` for errors. **None found.**

```bash
awk '$1!~"^#" && $1!="" && $13<0.001 {print $1"\t"($2-1)"\t"$3"\t"$6"___"$7}' $PFAM > $PFAM.bed
./add_value_to_table.py -i $PFAM.bed -a $GFF3.ORF_offset.txt \
  | awk -F'\t' '{print $6"\t"(($2*3)+$5)"\t"($3*3)+$5"\t"$4}' \
  > $FASTA.pfam_scan.CDS_coords.ORF_offset.bed
```

Adjust KEGG UniProt BLASTp top hit so they are correctly positioned along the original transcript (i.e. need to add length of 5-prime UTR to coords to adjust for ORF position).

```bash
./add_value_to_table.py \
    -i combined.manually_fixed.blastp_KEGG_UniProt_tophit.CDS_coords.bed \
    -a $GFF3.ORF_offset.txt \
  | awk -F'\t' '{print $6"\t"($2+$5)"\t"($3+$5)"\t"$4}' \
  > $FASTA.blastp_KEGG_UniProt_tophit.CDS_coords.ORF_offset.bed
```

Adjust possible start codons so they are correctly positioned along the original transcript (i.e. need to add length of 5-prime UTR to coords to adjust for ORF position).

```bash
./find_alternate_start_sites_in_a_protein.py -i $PEP \
  | awk -F'\t' '{print $1"\t"$2"\t"$2+1"\t"$4}' \
  > $PEP.start_codons.bed
./add_value_to_table.py -i $PEP.start_codons.bed -a $GFF3.ORF_offset.txt \
  | awk -F'\t' '{print $6"\t"(($2*3)+$5)"\t"($3*3)+$5"\t"$4}' \
  > $FASTA.start_codons.CDS_coords.ORF_offset.bed
```

Get the position along the transcript where the intron/exon boundaries are (don't need to adjust).

```bash
./splice_site_positions_along_CDS.py --use_ftype exon \
    -i combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.genome.gff3 \
  | sed -e 's/\.p[0-9]*\t/\t/' \
  > $FASTA.splice_sites.bed
```

Get the 5-prime UTR regions >20bp and BLASTx against RefSeq Complete

```bash
seqkit subseq \
    --bed <(cat combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.gff3.ORF_offset.txt \
  | awk '$2>=20{print $1"\t0\t"$2}' \
  | sed -e 's/\.p[0-9]*\t/\t/') \
  combined.manually_fixed.transcripts.fasta \
  > combined.manually_fixed.transcripts.5prime_UTR.min20bp.fasta

blastx \
  -query combined.manually_fixed.transcripts.5prime_UTR.min20bp.fasta \
  -db /scratch/timothy/databases/refseq/complete_204/complete.all.protein.faa.gz \
  -out combined.manually_fixed.transcripts.5prime_UTR.min20bp.fasta.blastx_complete.outfmt6 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -num_threads 36
```

Filter, merge and convert BLASTx results into bed format.

```bash
awk -F'\t' '$11<1e-5 { if($8>$7) {print $1"\t"$7-1"\t"$8"\t"$2} else {print $1"\t"$8-1"\t"$7"\t"$2} }' \
    combined.manually_fixed.transcripts.5prime_UTR.min20bp.fasta.blastx_complete.outfmt6 \
  | bedtools sort | bedtools merge -d -10 -c 4 -o first \
  | sed -e 's/_[^\t]*\t/\t/' \
  > combined.manually_fixed.transcripts.5prime_UTR.min20bp.fasta.blastx_complete.evalue-5.merged.bed
```



```bash
seqkit subseq \
    --bed <(cat combined.manually_fixed.blastp_KEGG_UniProt_tophit_fullAlnInfo.PEP_coords.bed \
              | awk '$2>=20{print $1"\t0\t"$2}') \
    combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep \
  > combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.5prime_terminus.min20aa.pep

blastp \
  -query combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.5prime_terminus.min20aa.pep \
  -db /scratch/timothy/databases/refseq/complete_204/complete.all.protein.faa.gz \
  -out combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.5prime_terminus.min20aa.pep.blastp_complete.outfmt6 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -num_threads 36

```

Filter, merge and convert BLASTp results into bed format.

```bash
awk -F'\t' '$11<1e-5 {print $1"\t"$7-1"\t"$8"\t"$2}' combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.5prime_terminus.min20aa.pep.blastp_complete.outfmt6 \
  | bedtools sort | bedtools merge -d -10 -c 4 -o first \
  | sed -e 's/_[^\t]*\t/\t/' \
  > combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.5prime_terminus.min20aa.pep.blastp_complete.evalue-5.merged.bed

./add_value_to_table.py \
    -i combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.5prime_terminus.min20aa.pep.blastp_complete.evalue-5.merged.bed \
    -a $GFF3.ORF_offset.txt \
  | awk -F'\t' '{print $6"\t"(($2*3)+$5)"\t"($3*3)+$5"\t"$4}' \
  > combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.5prime_terminus.min20aa.pep.blastp_complete.evalue-5.merged.ORF_offset.bed
```

### Download results from server

```bash
F="combined.manually_fixed.transcripts.fasta"
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/manually_correct_gene_models_transcript_visualization/"
rsync -zarv --prune-empty-dirs --relative --include="*/" \
 --include="5prime_terminus.blastp.bed" \
 --include="combined.manually_fixed.transcripts.5prime_UTR.min20bp.fasta.blastx_complete.evalue-5.merged.bed " \
 --include="*.coordsorted.bam*" \
 --include="$F" \
 --include="$F.fai" \
 --include="$F.blastp_KEGG_UniProt_tophit.CDS_coords.ORF_offset.bed" \
 --include="$F.pfam_scan.CDS_coords.ORF_offset.bed" \
 --include="$F.splice_sites.bed" \
 --include="$F.start_codons.CDS_coords.ORF_offset.bed" \
 --include="$F.transdecoder.forceStart.5prime_terminus.min20aa.pep.blastp_complete.evalue-5.merged.ORF_offset.bed" \
 --include="$F.transdecoder.forceStart.gff3" \
 --exclude="*" ${WD}/./ \
 . --dry-run
```
