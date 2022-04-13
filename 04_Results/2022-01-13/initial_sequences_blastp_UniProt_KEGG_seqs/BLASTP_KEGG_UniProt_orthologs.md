# BLASTP KR01 proteins against UniProt KEGG orthologs

BLAST KR01 proteins against KEGG proteins from UniProt. This will show us how our KR01 KEGG proteins compare to other proteins with the same KO numbers, i.e., are they the same length, too short (false positives?), or have long 5- or 3-prime termini?  

## Setup analysis directory

Link files with KR01 gene IDs and KO numbers + extra info we collected for each protein.

```bash
ln -s ../initial_sequences_info/Histidine_subset.txt
ln -s ../initial_sequences_info/Histidine_subset.txt.extra_info_added.txt
```

Link to KR01 proteins.

```bash
ln -s /scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa
```

Link to a combined UniProt protein file (SwissProt + TrEMBL).

```bash
ln -s ../../../../0016_Paulinella_micropora_KR01_BioPhysics/03_Analysis/2021-04-03/BioPhysics/uniprot.faa
ln -s ../../../../0016_Paulinella_micropora_KR01_BioPhysics/03_Analysis/2021-04-03/BioPhysics/uniprot.faa.fai
```

Link to custom script for filtering BLAST results.

```bash
ln -s ../../../02_Scripts/blast_top_hits.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/ncbi-blast-2.10.1+/bin"
export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
# Needs seqkit: /home/timothy/programs/seqkit_v0.15.0/seqkit
# Also needs gnu parallel
```

### Get list of target IDs

Get unique list of KEGG IDs that are in our target pathways.

```bash
cat Histidine_subset.txt.extra_info_added.txt \
  | awk -F'\t' '{print $1}' | sort | uniq \
  > KEGG_ids.txt
```

Get unique list of KR01 gene IDs that are in our target pathways.

```bash
cat Histidine_subset.txt.extra_info_added.txt \
  | awk -F'\t' '{print $1"\t"$2}' | sort | uniq \
  > KEGG_and_GENE_id_pairs.txt
```

### Extract UniProt sequences annotated with the KO numbers that we are interested in

For each KEGG pathway;

- Get associated UniProt seq IDs
- Get UniProt protein seqs from UniProt seq IDs

```bash
get_KEGG_UniProt_seqs() {
        # Print all info to log file
        exec 1> "${1}.get_KEGG_UniProt_seqs.log" 2>&1
        
        K="$1"
        SEQKIT="/home/timothy/programs/seqkit_v0.15.0/seqkit"
        
        echo "Getting IDs for $K"
        
        wget i-q https://www.kegg.jp/kegg-bin/uniprot_list?ko=$K -O $K.uniprot_id_list.html
        cat $K.uniprot_id_list.html \
          | awk '$0~"http://www.uniprot.org/uniprot/" {T=$0; getline; print T$1}' \
          | sed -e "s@ <td><a href='http://www.uniprot.org/uniprot/@@" -e "s@'>.*</a> </td><td>@\t@" -e 's@</td>@@' \
          | awk '{ if ($2~$1) {print "tr|"$1"|"$2} else {print "sp|"$1"|"$2} }' \
          > $K.uniprot_id_list.cleaned.txt
        
        "$SEQKIT" faidx \
            --infile-list $K.uniprot_id_list.cleaned.txt \
            --out-file $K.uniprot_id_list.cleaned.faa uniprot.faa
}

export -f get_KEGG_UniProt_seqs
parallel -j 4 get_KEGG_UniProt_seqs :::: KEGG_ids.txt
```

Make BLAST DB of UniProt KEGG seqs

```bash
for F in *.uniprot_id_list.cleaned.faa
do
        makeblastdb -dbtype prot -in $F
done > makeblastdb.log 2>&1
```

### Compare KR01 proteins against UniProt KEGG seqs

Get KR01 pep seqs into individual files.

```bash
PEP="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa"
for F in `cut -f2 KEGG_and_GENE_id_pairs.txt | awk '$1~"^g"'`
do
        samtools faidx "$PEP" "Paulinella_micropora_KR01_nuclear___$F" > "$F.pep.faa"
done
```

BLASTP KR01_genes against UniProt_KEGG_seqs.

- For each KEGG<->KR01_nuc_pep pair run BLASTP and take just the top hit.
- Convert pep coods to bed formatted CDS coords: (pep_start-1)\*3 and pep_stop\*3
- Also print subject top hit pep coords + length (useful for computing subject coverage later).

```bash
blastp_KEGG_UniProt_seqs() {
        K=$(echo $1 | awk '{print $1}')
        B=$(echo $1 | awk '{print $2}')

        # Print all info to log file
        exec 1> "${B}_vs_${K}.blastp.log" 2>&1

        echo "Running BLASTp $B -vs- $K"

        blastp -query "$B.pep.faa" -db "$K.uniprot_id_list.cleaned.faa" \
          -out "${B}_vs_${K}.blastp.outmft6" \
          -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        ./blast_top_hits.py -i "${B}_vs_${K}.blastp.outmft6" \
          | awk -F'\t' -vK="$K" '{print $1"\t"($7-1)*3"\t"$8*3"\t"K":"$9"-"$10"/"$14}' \
          > "${B}_vs_${K}.blastp.outmft6.tophit.CDS_coords.bed"
}

export -f blastp_KEGG_UniProt_seqs
parallel -j 10 blastp_KEGG_UniProt_seqs :::: <(awk '$2~"^g"' KEGG_and_GENE_id_pairs.txt)
```

Cat all CDS bed features together.

```bash
cat *.blastp.outmft6.tophit.CDS_coords.bed > KEGG_and_GENE_id_pairs.txt.blastp_KEGG_UniProt_tophit.CDS_coords.bed
```

## Results

### Download results from server

```bash
scp timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2022-01-13/initial_sequences_blastp_UniProt_KEGG_seqs/KEGG_and_GENE_id_pairs.txt.blastp_KEGG_UniProt_tophit.CDS_coords.bed .
```
