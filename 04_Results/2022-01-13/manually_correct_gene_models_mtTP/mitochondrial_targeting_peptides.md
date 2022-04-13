# Predict mitochondrial targeting peptides in protein sequences

## Setup analysis directory

Link to manually corrected proteins.

```bash
ln -s ../manually_correct_gene_models/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep
```

Link to custom scripts for text parsing.

```bash
ln -s ../../../02_Scripts/add_value_to_table.py 
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/targetp-2.0/bin"
```

## Mitochondrial Targeting Peptides (mtTP)

Assess if a sequence has a mitochondrial targeting peptide using `targetp`.

**Manually corrected proteins**

```bash
PEP="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep"
targetp -fasta "$PEP" -gff3 -prefix "$PEP.targetp"
```

If a protein has "mTP" in the 2nd column then it has a mtTP. **Neither of the proteins have a mtTP**

