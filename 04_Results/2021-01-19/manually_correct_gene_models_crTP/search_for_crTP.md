# Search for crTP in manually corrected proteins

Search the new manually corrected gene models for crTP sequences (since we have reprodicted the proteins and some might have lost or gained a crTP). Use the KR01 crTP alignment that was used by Eva Nowack for the 2020 KR01 genome paper. 

## Setup analysis directory

Link to manually corrected proteins.

```bash
ln -s ../manually_correct_gene_models/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep
```

Link to crTP HMM that we created when we were analyzing the "initial" proteins.

```bash
ln -s ../initial_sequences_plot_read_mapping/KR01_crTP_StartingSet_crTPonly_ManualCorr.aln.fasta.hmm
```

Link to custom script for parsing text files.

```bash
ln -s ../../../02_Scripts/add_value_to_table.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/hmmer-3.1b2/bin"
```

## Search for crTP

Use HMM we built earlier to search for crTP in the manually corrected gene models. Keep coords for the pep sequences (i.e., don't transform to CDS coords); pep coords are used in Master Excel table for other functional annotations. Simple print to screen and manually add to Master Excel table; simple since there are only a few proteins with crTP.

```bash
P="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
## Search for crTP in new rebuilt KEGG proteins.
hmmsearch --domtblout $P.crTP.domtblout $A.fasta.hmm $P

## Print crTP hits with e-value <1e-5 and >50% HMM coverage
grep -v '^#' $P.crTP.domtblout | awk '$12<1e-5 && ($17-$16)+1 > $6*0.5 {print $1"\t"$18-1"\t"$19}'

## Print crTP info in the order from master Excel table (e-value <1e-5)
## Dont filter by coverage (less stringent but still might be informative)
awk -F'\t' '{print $3}' ../blastp_UniProt_KEGG_seqs/SuppTables.transcript_info.txt \
	| ./add_value_to_table.py -a <(grep -v '^#' $P.crTP.domtblout | awk '$12<1e-5 {print $1"\t"$18-1"\t"$19"\t"$16-1"\t"$17}')
```

Header/column order of output file.

```bash
#      1                 2         3           4               5         6       7       8     9    10  11     12       13       14    15     16   17    18     19    20    21   22       23
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
```

