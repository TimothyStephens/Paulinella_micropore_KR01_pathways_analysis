# Phylogenetic analysis of manually corrected proteins

Perform phylogenetic analysis using the manually corrected proteins. This analysis is to see if the proteins are of host, endosymbiont, or foreign (HGT) origin. 

## Setup analysis directory

Link manually corrected proteins.

```bash
ln -s ../manually_correct_gene_models/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep 
```

Link to custom scripts for text parsing and for filtering `fasta` files.

```bash
ln -s ../../../02_Scripts/grepf_column.py
ln -s ../../../02_Scripts/grepf_fasta.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
```

## Prepare proteins for analysis

Get each manually verified protein sequence into its own `fasta` file.

```bash
awk '{print $1}' combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep \
  | seqkit fx2tab | awk '{print ">"$1"\n"$2 > $1".pep.faa"}'

ls -1 --color=none MSTRG.* > files2run.txt
```

## 01 - BLASTP against taxonomically diverse database

For each manually corrected protein, compare it against a taxonomically diverse database using BLASTP. Database is split into four parts (Bacteria, Metazoa, MMETSP, and Other), BLAST against each separately using very lose e-value cutoff (will apply a proper filter in later steps). Use a loose filter now as small proteins often have weak e-values even for good hits, so keep all hits now so later we can play with the filtering (if needed).

```bash
./run_01_blastp_DanaDB.sh
```

See if there are any obvious error in the log files from BLAST.

```bash
grep -i 'error\|warn\|fault\|seg' *.blast.log
cat *.blast.log | grep -v 'CMD:' | grep -v '#'
```

Nothing was returned from the log files so it looks like they all completed OK. 

## 02 - Process BLAST results

For each manually corrected protein, combine the BLAST results from each database, filter the hits, down sample the subject sequences using taxonomic information, and extract the proteins from the final set of selected sequences (to be used for phylogenetic analysis). 

The `run_02_process_blast.sh` script:

- Extracts all subject sequences (with hits) across all query proteins + add manually corrected proteins
    - This reduces the number of sequences we need to handle in subsequent steps, so it speeds up the analysis if we only have to parse the large database files once at the beginning. 
    - NOTE: these is some redundancy across the four database files (mostly between the MMETSP and OTHER databases I think) so the script will remove proteins with identical names AND sequences if they do occur.
- For each manually corrected (query) protein:
    1. Concatenate the BLAST results from the different databases.
    2. Filter *e*-value < **1e-5**
    3. Sort hits by query_id then by bitscore
    4. Down sample subject sequences (with `down_sample_taxa.sh`) using taxonomic information about the different sequences. The script will:
        1. Calculate the number of taxonomic 'Classes' present in the hit subjects.
        2. Calculate the max number of 'Genera' per 'Class' to keep (`numGen=$((180/numClass))`)
        3. If the number of 'Genera' is < 25: retrieve 2 subject sequences per 'Genera'; else retrieve 1 subject sequence per 'Genera'
    5. Get the hits for the down sampled subject sequences
    6. Get the protein sequences for the down sampled subject sequences + the current manually corrected (query) protein

```bash
./run_02_process_blast.sh
```

Check log file `run_02_process_blast.sh.log.*` for errors.

A few samples had no hits that passed the *e*-value filter. This caused some 'division by 0' errors in the log file but that OK. Some of these sequences were short or had long regions of simple repeats, which is likely why they had no hits. Can probably safely ignore these sequences for now. 

## 03 - Align proteins and build trees

For each manually corrected protein, align the proteins that we just collected (and down sampled), and build a phylogenetic tree. 

The `run_03_build_trees.sh` script:

- Runs `mafft` on the KR01 protein + down sampled BLASTP hits. 
- Runs `iqtree` on the alignment to build a maximum likelihood tree with substitution model picked automatically and bootstrap support added via 2000 ultrafast bootstraps. 
  - If an alignment had <= 3 sequences it couldn't/didn't have a tree built using `iqtree`. 

```bash
./run_03_build_trees.sh
```

Check log file `run_03_build_trees.sh.log.*` for errors. **None found.**

Check `*.build_trees.log` for errors. **None found.**

```bash
cat *.build_trees.log | grep -v '#\|CMD'
```

Check `mafft` and `iqtree` log files for errors. 

```bash
grep -i 'error\|warn' *.mafft.log
grep -i 'error\|warn' *.mafft.aln.log
```

No errors found for `mafft`. Lots of warnings were found for `iqtree` but most are unavoidable or just warning the user about features of the data, so are probably not critical. 

## Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/manually_correct_gene_models_phylogenetic_analysis_FullProtein"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="*.sh*" --include="*.log" --include="*.aln" --include="*.contree" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

