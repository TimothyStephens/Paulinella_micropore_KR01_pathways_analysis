# Phylogenetic analysis of crTP sequences in photosynthetic *Paulinella*

## Setup analysis directory

Link KR01 + *P. chromatophora* predicted proteins - Will use to re-predict crTP locations since we only know which ones had a crTP not where it was along the protein. 

```bash
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa
ln -s ../../../../../genome_data/Paulinella_chromatophora_CCAC0185/nuclear_genome/databases/Paulinella_chromatophora_CCAC0185_nuclear.pep.faa
```

Link to crTP predictions for manually corrected proteins.

```bash
ln -s ../manually_correct_gene_models_crTP/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep.crTP.domtblout
ln -s ../manually_correct_gene_models_crTP/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep
```

Link to crTP HMM that we created when we were analyzing the "initial" proteins.

```bash
ln -s ../initial_sequences_plot_read_mapping/KR01_crTP_StartingSet_crTPonly_ManualCorr.aln.fasta.hmm
```

Custom scripts.

```bash
ln -s ../../../02_Scripts/add_value_to_table.py
ln -s ../../../02_Scripts/grepf_fasta.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/hmmer-3.1b2/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
HMM="KR01_crTP_StartingSet_crTPonly_ManualCorr.aln.fasta.hmm"
```

## Search for crTP in KR01 + *P. chromatophora*

Use HMM we built earlier to search for crTP in KR01 + *P. chromatophora* gene models. Make results into a `bed`-like formatted file.

>  seqid \t start \t stop \t hmm_start \t hmm_stop

```bash
P="Paulinella_chromatophora_CCAC0185_nuclear.pep.faa"
hmmsearch --domtblout $P.crTP.domtblout $HMM $P
grep -v '^#' $P.crTP.domtblout | awk '$12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17}' > $P.crTP.domtblout.bed
seqkit subseq --bed $P.crTP.domtblout.bed $P > $P.crTP_regions.faa
```

Rename sequences that are annotated with KO numbers that we are interested in.

```bash
awk '{print $1}' $P.crTP_regions.faa | sed -e 's/:.*//' > $P.crTP_regions.renamed.faa
awk -F'\t' '$3!="" {print $3"\t"$1}' KO_to_PepID.txt | sort | uniq | sed -e 's/-size[^_]*_m\./-m/' | awk -F'\t' '{print "sed -i -e \"s/"$1"/"$1"_"$2"/\" Paulinella_chromatophora_CCAC0185_nuclear.pep.faa.crTP_regions.renamed.faa"}' | bash
```

KR01

```bash
P="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa"
hmmsearch --domtblout $P.crTP.domtblout $HMM $P
grep -v '^#' $P.crTP.domtblout | awk '$12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17}' > $P.crTP.domtblout.bed
seqkit subseq --bed $P.crTP.domtblout.bed $P > $P.crTP_regions.faa
```

Rename sequences that are annotated with KO numbers that we are interested in.

```bash
awk '{print $1}' $P.crTP_regions.faa | sed -e 's/:.*//' > $P.crTP_regions.renamed.faa
awk '{print $1}' $P.crTP_regions.faa | sed -e 's/:.*//' | seqkit fx2tab | grep -v -f <(awk -F'\t' '{print "___"$2"_"}' ../../../01_Data/2021-01-19/*.txt | sort | uniq) | seqkit tab2fx > $P.crTP_regions.renamed.faa
```

Manually corrected proteins.

```bash
P="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
grep -v '^#' $P.crTP.domtblout | awk '$12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17}' > $P.crTP.domtblout.bed
seqkit subseq --bed $P.crTP.domtblout.bed $P > $P.crTP_regions.faa
awk '{print $1}' $P.crTP_regions.faa | sed -e 's/:.*//' > $P.crTP_regions.renamed.faa
```

## crTP phylogenetic analysis

Combine crTP protein sequences from KR01 (- initial proteins) + *P. chromatophora* + KR01 manually corrected proteins.  

```bash
cat *.crTP_regions.renamed.faa > Paulinella_crTP.pep.faa
```



