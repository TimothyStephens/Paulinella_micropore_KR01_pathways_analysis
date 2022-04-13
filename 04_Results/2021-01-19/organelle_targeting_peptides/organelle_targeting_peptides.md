# Predict organelle targeting peptides in protein sequences

## Setup analysis directory

Link to manually corrected + KR01 + *P. Chromatophora* (in Dana DB) + *P. Chromatophora* (my naming format) + *P. Ovalis* protein files.

```bash
ln -s ../manually_correct_gene_models/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep

ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa

ln -s ../../../../../genome_data/Paulinella_chromatophora_CCAC0185/nuclear_genome/CCAC0185-TGI-SCAFF.transdecoder.pep.fasta.gz

ln -s ../../../../../genome_data/Paulinella_chromatophora_CCAC0185/nuclear_genome/databases/Paulinella_chromatophora_CCAC0185_nuclear.pep.faa

ln -s ../../../../../genome_data/Paulinella_ovalis_6_SingleCells/nuclear_genome/databases/Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep.faa
```

Link to crTP HMM that we created when we were analyzing the "initial" proteins.

```bash
ln -s ../initial_sequences_plot_read_mapping/KR01_crTP_StartingSet_crTPonly_ManualCorr.aln.fasta.hmm
```

Link to custom scripts for text parsing.

```bash
ln -s ../../../02_Scripts/add_value_to_table.py 
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/targetp-2.0/bin"
export PATH="$PATH:/home/timothy/programs/hmmer-3.1b2/bin"
HMM="KR01_crTP_StartingSet_crTPonly_ManualCorr.aln.fasta.hmm"
```

## Reformat names of *P. chromatophora* sequences to match what is in Dana's DB

Need to format the names of the *P. chromatophora* sequences to match what they are in the database + trees; need to replace `|` with `_` and add taxonomic info to sequence header. 

```bash
zcat CCAC0185-TGI-SCAFF.transdecoder.pep.fasta.gz | awk '{print $1}' | sed -e 's/>/>Eukaryota-Cercozoa-NA-Euglyphida-Paulinellidae-Paulinella_chromatophora_/' -e 's/|/_/' > CCAC0185-TGI-SCAFF.transdecoder.pep.fasta
```

## Mitochondrial Targeting Peptides (mtTP)

Assess if a sequence has a mitochondrial targeting peptide using `targets`.

**Manually corrected proteins**

```bash
PEP="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
targetp -fasta "$PEP" -gff3 -prefix "$PEP.targetp"
```

**KR01**

```bash
PEP="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa"
targetp -fasta "$PEP" -gff3 -prefix "$PEP.targetp"
```

***P. Chromatophora*** (in Dana DB)

```bash
PEP="CCAC0185-TGI-SCAFF.transdecoder.pep.fasta"
targetp -fasta "$PEP" -gff3 -prefix "$PEP.targetp"
```

***P. Chromatophora*** (my naming format)

```bash
PEP="Paulinella_chromatophora_CCAC0185_nuclear.pep.faa"
targetp -fasta "$PEP" -gff3 -prefix "$PEP.targetp"
```

***P. Ovalis***

```bash
PEP="Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep.faa"
targetp -fasta "$PEP" -gff3 -prefix "$PEP.targetp"
```

## Chromatophore Targeting Peptide (crTP)

Assess if a sequence has a chromatophore targeting peptide using `hmmer` and the KR01 crTP HMM.

Output bed formatted file has the following columns:

> seqid \t hit_start \t hit_stop \t hmm_start \t hmm_stop \t hmm_hit_len

**Manually corrected proteins**

```bash
PEP="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
hmmsearch --domtblout $PEP.crTP.domtblout $HMM $PEP
awk '$1!~"#" && $12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17"\t"$17-($16-1)}' \
  $PEP.crTP.domtblout > $PEP.crTP.domtblout.bed
```

**KR01**

```bash
PEP="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa"
hmmsearch --domtblout $PEP.crTP.domtblout $HMM $PEP
awk '$1!~"#" && $12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17"\t"$17-($16-1)}' \
  $PEP.crTP.domtblout > $PEP.crTP.domtblout.bed
```

***P. Chromatophora*** (in Dana DB)

```bash
PEP="CCAC0185-TGI-SCAFF.transdecoder.pep.fasta"
hmmsearch --domtblout $PEP.crTP.domtblout $HMM $PEP
awk '$1!~"#" && $12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17"\t"$17-($16-1)}' \
  $PEP.crTP.domtblout > $PEP.crTP.domtblout.bed
```

***P. Chromatophora*** (my naming format)

```bash
PEP="Paulinella_chromatophora_CCAC0185_nuclear.pep.faa"
hmmsearch --domtblout $PEP.crTP.domtblout $HMM $PEP
awk '$1!~"#" && $12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17"\t"$17-($16-1)}' \
  $PEP.crTP.domtblout > $PEP.crTP.domtblout.bed
```

***P. Ovalis***

```bash
PEP="Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep.faa"
hmmsearch --domtblout $PEP.crTP.domtblout $HMM $PEP
awk '$1!~"#" && $12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17"\t"$17-($16-1)}' \
  $PEP.crTP.domtblout > $PEP.crTP.domtblout.bed
```

## Combine crTP and mtTP predictions into a single table

Combine crTP and mtTP results into a single table listing the status of each sequence in the different datasets.

Output bed formatted file has the following columns:

> seqid \t has_crTP \t has_mtTP
>
> has_crTP: either "No" or the length of the crTP
>
> has_mtTP: either "No" or "Yes"

**Manually corrected proteins**

```bash
PEP="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
awk '$1~"^>" {print $1}' $PEP \
  | sed -e 's/>//' \
  | ./add_value_to_table.py -d "No" -a <(awk '!seen[$1]++' $PEP.crTP.domtblout.bed | cut -f1,6) \
  | ./add_value_to_table.py -d "No" -a <(awk -F'\t' '$1!~"^#" && $2=="mTP" {print $1"\tYes"}' $PEP.targetp_summary.targetp2) \
  | awk -F'\t' 'BEGIN{print "seqid\tcrTP\tmtTP"} {print}' \
  > $PEP.organelle_TP.txt
```

**KR01**

```bash
PEP="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa"
awk '$1~"^>" {print $1}' $PEP \
  | sed -e 's/>//' \
  | ./add_value_to_table.py -d "No" -a <(awk '!seen[$1]++' $PEP.crTP.domtblout.bed | cut -f1,6) \
  | ./add_value_to_table.py -d "No" -a <(awk -F'\t' '$1!~"^#" && $2=="mTP" {print $1"\tYes"}' $PEP.targetp_summary.targetp2) \
  | awk -F'\t' 'BEGIN{print "seqid\tcrTP\tmtTP"} {print}' \
  > $PEP.organelle_TP.txt
```

***P. Chromatophora*** (in Dana DB; If multiple crTP hits, just just the first one)

```bash
PEP="CCAC0185-TGI-SCAFF.transdecoder.pep.fasta"
awk '$1~"^>" {print $1}' $PEP \
  | sed -e 's/>//' \
  | ./add_value_to_table.py -d "No" -a <(awk '!seen[$1]++' $PEP.crTP.domtblout.bed | cut -f1,6) \
  | ./add_value_to_table.py -d "No" -a <(awk -F'\t' '$1!~"^#" && $2=="mTP" {print $1"\tYes"}' $PEP.targetp_summary.targetp2) \
  | awk -F'\t' 'BEGIN{print "seqid\tcrTP\tmtTP"} {print}' \
  > $PEP.organelle_TP.txt
```

***P. Chromatophora*** (my naming format)

```bash
PEP="Paulinella_chromatophora_CCAC0185_nuclear.pep.faa"
awk '$1~"^>" {print $1}' $PEP \
  | sed -e 's/>//' \
  | ./add_value_to_table.py -d "No" -a <(awk '!seen[$1]++' $PEP.crTP.domtblout.bed | cut -f1,6) \
  | ./add_value_to_table.py -d "No" -a <(awk -F'\t' '$1!~"^#" && $2=="mTP" {print $1"\tYes"}' $PEP.targetp_summary.targetp2) \
  | awk -F'\t' 'BEGIN{print "seqid\tcrTP\tmtTP"} {print}' \
  > $PEP.organelle_TP.txt
```

***P. Ovalis***

```bash
PEP="Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep.faa"
awk '$1~"^>" {print $1}' $PEP \
  | sed -e 's/>//' \
  | ./add_value_to_table.py -d "No" -a <(awk '!seen[$1]++' $PEP.crTP.domtblout.bed | cut -f1,6) \
  | ./add_value_to_table.py -d "No" -a <(awk -F'\t' '$1!~"^#" && $2=="mTP" {print $1"\tYes"}' $PEP.targetp_summary.targetp2) \
  | awk -F'\t' 'BEGIN{print "seqid\tcrTP\tmtTP"} {print}' \
  > $PEP.organelle_TP.txt
```

## Print results in the order of sequences in the master Excel document

```bash
cat GeneID.Host.txt | ./add_value_to_table.py -d "No" -a <(awk -F'\t' '$1!~"^#" && $2=="mTP" {print $1"\tmtTP"}' Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep.faa.targetp_summary.targetp2)
```

```bash
awk '{ if($1==""){print "."} else{print $1}}' GeneID.Chrom.txt | ./add_value_to_table.py -a <(awk '!seen[$1]++' CCAC0185-TGI-SCAFF.transdecoder.pep.fasta.crTP.domtblout.bed | sed -e 's/Eukaryota-Cercozoa-NA-Euglyphida-Paulinellidae-Paulinella_chromatophora_//' | cut -f1,6)
```

```bash
awk '{ if($1==""){print "."} else{print $1}}' GeneID.Chrom.txt | ./add_value_to_table.py -a <(awk -F'\t' '$1!~"#" {print $1"\t"$2}' CCAC0185-TGI-SCAFF.transdecoder.pep.fasta.targetp_summary.targetp2 | sed -e 's/Eukaryota-Cercozoa-NA-Euglyphida-Paulinellidae-Paulinella_chromatophora_//')
```

