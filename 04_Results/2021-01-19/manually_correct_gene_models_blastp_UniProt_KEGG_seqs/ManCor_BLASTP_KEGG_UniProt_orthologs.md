# BLASTP search the manually corrected proteins against UniProt KEGG orthologs

BLAST manually corrected KR01 proteins against KEGG proteins from UniProt. This will show us how our KR01 KEGG proteins compare to other proteins with the same KO numbers, i.e., are they the same length, too short (false positives?), or have long 5- or 3-prime termini? Do they even still have hits to other KEGG proteins after manual correction? (can happen when the initial gene was severely mis-predicted)  .

## Setup analysis directory

Link to manually corrected proteins.

```bash
ln -s ../manually_correct_gene_models/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep 
```

Link to the UniProt KEGG Orthologs that we have already retrieved for analysis using the "initial" proteins.

```bash
for F in ../initial_sequences_blastp_UniProt_KEGG_seqs/*.uniprot_id_list.cleaned.faa*;
do
	ln -s $F
done
```

Link to custom scripts for getting BLAST top hits and parsing text files.

```bash
ln -s ../../../02_Scripts/blast_top_hits.py
ln -s ../../../02_Scripts/add_value_to_table.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
```

## BLAST against UniProt

Split the KR01 manually corrected proteins into separate files (one seq per file). Also clean description from `fasta` header and make sequence into single line (from multiline) `fasta`.

```bash
PEP="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
seqkit fx2tab "$PEP" | awk -F'\t' '{split($1,a," "); print ">"a[1]"\n"$2 > a[1]".pep.faa"}'
```

BLASTP KR01_pep against UniProt_KEGG_pep

For each KEGG<->KR01_pep pair run BLASTP and take just the top hit.

```bash
blastp_KEGG_UniProt_seqs() {
        K=$(echo $1 | awk '{print $1}')
        B=$(echo $1 | awk '{print $2}')
        Q=$(echo $B.pep.faa)

        # Print all info to log file
        exec 1> "${B}_vs_${K}.blastp.log" 2>&1

        echo "Running BLASTp $B -vs- $K"

        blastp -query "$Q" -db "$K.uniprot_id_list.cleaned.faa" \
            -out "${Q%*.pep.faa}_vs_${K}.blastp.outmft6" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        ./blast_top_hits.py -i "${Q%*.pep.faa}_vs_${K}.blastp.outmft6" \
          | awk -F'\t' -vK="$K" '{print $1"\t"($7-1)*3"\t"$8*3"\t"K":"$9"-"$10"/"$14}' \
          > "${Q%*.pep.faa}_vs_${K}.blastp.outmft6.tophit.CDS_coords.bed"
}

export -f blastp_KEGG_UniProt_seqs
parallel -j 10 blastp_KEGG_UniProt_seqs :::: <(cut -f1,3 SuppTables.transcript_info.txt | awk '$2~"MSTRG"')
```

Concatenate all CDS bed features together.

```bash
cat MSTRG.*.blastp.outmft6.tophit.CDS_coords.bed \
  > combined.manually_fixed.blastp_KEGG_UniProt_tophit.CDS_coords.bed
```

## Reformat top hit results

Get query (transcript) AND subject (KEGG) info in a bed-like format (0-based coords; protein coords for both files).

- Filter blast hits (e-vlaue < 1e-3 [weak cut-off])

qseqid [tab] qstart [tab] qend [tab] qlen [tab] KEGG_id [tab] sstart [tab] send [tab] slen

```bash
while read pair; 
do
        K=$(echo $pair | awk '{print $1}');
        B=$(echo $pair | awk '{print $2}');
        Q=$(echo $B.pep.faa);
        cat "${Q%*.pep.faa}_vs_${K}.blastp.outmft6" \
          | awk -F'\t' '$11<1e-5' \
          | ./blast_top_hits.py \
          | awk -F'\t' -vK="$K" '{print $1"\t"($7-1)"\t"$8"\t"$13"\t"K"\t"($9-1)"\t"$10"\t"$14}'
done < <(cut -f1,3 SuppTables.transcript_info.txt | awk '$2~"MSTRG"' | sort | uniq) \
  > combined.manually_fixed.blastp_KEGG_UniProt_tophit_fullAlnInfo.PEP_coords.bed
```

Add BLASTP to KEGG seqs info to table in same format as the Master Excel spread sheet.

```bash
PEP="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
BED="combined.manually_fixed.blastp_KEGG_UniProt_tophit_fullAlnInfo.PEP_coords.bed"
#               |--------------new_transcript_ORF-------|       |--------KEGG_top_hit-------------------|
# Gene_id [tab] qseqid [tab] qstart [tab] qend [tab] qlen [tab] sseqid [tab] sstart [tab] send [tab] slen [tab] new_transcript_ORF_type
cat SuppTables.transcript_info.txt \
 | ./add_value_to_table.py -c 3 -a <(awk '{print $1"\t"$0}' "$BED") \
 | ./add_value_to_table.py -c 3 -a <(awk '$1~">" {print $1"\t"$4}' "$PEP" | sed -e 's/>//' -e 's/type://') \
 > SuppTables.transcript_info.txt.blastKEGG_info.txt
```

`SuppTables.transcript_info.txt` is a file that lists the "initial" and "new" protein IDs in the order that they are in in my Master Excel table. The order in this file helps transfer the results into the spread sheet. 

## Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/manually_correct_gene_models_blastp_UniProt_KEGG_seqs"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="SuppTables.transcript_info.txt*" --include="combined.manually_fixed.blastp_KEGG_UniProt_tophit*" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```
