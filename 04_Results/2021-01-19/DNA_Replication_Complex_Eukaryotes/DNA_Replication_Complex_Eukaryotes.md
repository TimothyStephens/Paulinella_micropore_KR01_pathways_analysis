# Phylogenetic analysis of genes annotated as part of the Eukaryotic DNA Replication pathway

## Setup analysis directory

Link to KO annotations for the KR01 nuclear genome.

```bash
ln -s /scratch/timothy/Paulinella_micropora_KR01/GeneAnnotation/Paulinella_micropora_KR01_Nuclear/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa.kofamscan_results.mapper.txt
ln -s /scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa
```

Link custom scripts for text parsing.

```bash
ln -s ../../../02_Scripts/grepf_column.py
ln -s ../../../02_Scripts/add_value_to_table.py
ln -s ../../../02_Scripts/groupby.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
```

## Get KO numbers that are part of the eukaryotic pathway

Get KO numbers associated with the DNA Replication pathway (ko03030).

```bash
grep -A 60 'ko03030' /scratch/timothy/databases/KEGG/ko00001.keg | awk 'NR>1{L=""; for(i=3; i<=NF; i++){L=L" "$i}; print $2"\t"L }' > ko03030.KO_numbers.txt 
cp ko03030.KO_numbers.txt ko03030.KO_numbers.Eukayotic.txt
```

Manually sort through`ko03030.KO_numbers.Eukayotic.txt` and remove IDs that are not part of the Eukaryotic pathway. Use https://www.genome.jp/pathway/ko03030 to identify which pathways each ID is associated with. 

KO's selected:

```less
K02320	POLA1; DNA polymerase alpha subunit A [EC:2.7.7.7]
K02321	POLA2; DNA polymerase alpha subunit B
K02684	PRI1; DNA primase small subunit [EC:2.7.7.102]
K02685	PRI2; DNA primase large subunit
K02327	POLD1; DNA polymerase delta subunit 1 [EC:2.7.7.7]
K02328	POLD2; DNA polymerase delta subunit 2
K03504	POLD3; DNA polymerase delta subunit 3
K03505	POLD4; DNA polymerase delta subunit 4
K02324	POLE; DNA polymerase epsilon subunit 1 [EC:2.7.7.7]
K02325	POLE2; DNA polymerase epsilon subunit 2 [EC:2.7.7.7]
K02326	POLE3; DNA polymerase epsilon subunit 3 [EC:2.7.7.7]
K03506	POLE4; DNA polymerase epsilon subunit 4 [EC:2.7.7.7]
K02540	MCM2; DNA replication licensing factor MCM2 [EC:3.6.4.12]
K02541	MCM3; DNA replication licensing factor MCM3 [EC:3.6.4.12]
K02212	MCM4, CDC54; DNA replication licensing factor MCM4 [EC:3.6.4.12]
K02209	MCM5, CDC46; DNA replication licensing factor MCM5 [EC:3.6.4.12]
K02542	MCM6; DNA replication licensing factor MCM6 [EC:3.6.4.12]
K02210	MCM7, CDC47; DNA replication licensing factor MCM7 [EC:3.6.4.12]
K07466	RFA1, RPA1, rpa; replication factor A1
K10739	RFA2, RPA2; replication factor A2
K10741	RPA4; replication factor A4
K10740	RPA3; replication factor A3
K04802	PCNA; proliferating cell nuclear antigen
K10754	RFC1; replication factor C subunit 1
K10755	RFC2_4; replication factor C subunit 2/4
K10756	RFC3_5; replication factor C subunit 3/5
K03469	rnhA, RNASEH1; ribonuclease HI [EC:3.1.26.4]
K10743	RNASEH2A; ribonuclease H2 subunit A [EC:3.1.26.4]
K10744	RNASEH2B; ribonuclease H2 subunit B
K10745	RNASEH2C; ribonuclease H2 subunit C
K10742	DNA2; DNA replication ATP-dependent helicase Dna2 [EC:3.6.4.12]
K04799	FEN1, RAD2; flap endonuclease-1 [EC:3.-.-.-]
K10747	LIG1; DNA ligase 1 [EC:6.5.1.1 6.5.1.6 6.5.1.7]
```

Get KR01 proteins that are annotated with the KO numbers that we just collected.

```bash
F="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa.kofamscan_results.mapper.txt"
./grepf_column.py -c 2 -f <(cut -f1 ko03030.KO_numbers.Eukayotic.txt) -i <(awk -F'\t' '$2!=""' $F) > KR01.ko03030_Eukayotic.txt
```

Print KO numbers without annotated KR01 proteins.

```bash
./grepf_column.py -v -f <(cut -f2 KR01.ko03030_Eukayotic.txt | sort | uniq) -i ko03030.KO_numbers.Eukayotic.txt
```

>K03506	POLE4; DNA polymerase epsilon subunit 4 [EC:2.7.7.7]
>
>K02210	MCM7, CDC47; DNA replication licensing factor MCM7 [EC:3.6.4.12]
>
>K10741	RPA4; replication factor A4

## Separate proteins into different files and BLAST against database

Get each KR01 protein into a different file ready for BLASTP.

```bash
FAA="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa"
while read SEQ;
do
	seqkit faidx $FAA $SEQ > $SEQ.pep.faa;
done < <(cut -f1 KR01.ko03030_Eukayotic.txt)
```

Create list of protein files to run.

```bash
while read SEQ; do echo "$SEQ.pep.faa"; done < <(cut -f1 KR01.ko03030_Eukayotic.txt) > files2run.txt
```

## 01 - BLASTP against taxonomically diverse database

For each KR01 protein, compare it against a taxonomically diverse database using BLASTP. Database is split into four parts (Bacteria, Metazoa, MMETSP, and Other), BLAST against each separately using very lose e-value cutoff (will apply a proper filter in later steps). Use a loose filter now as small proteins often have weak e-values even for good hits, so keep all hits now so later we can play with the filtering (if needed).

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

For each KR01 protein, combine the BLAST results from each database, filter the hits, down sample the subject sequences using taxonomic information, combine the results for each protein annotated to the same KO number, and extract the proteins from the final set of selected sequences (to be used for phylogenetic analysis). 

The `run_02_process_blast.sh` script:

- Extracts all subject sequences (with hits) across all query proteins + add KR01 proteins
  - This reduces the number of sequences we need to handle in subsequent steps, so it speeds up the analysis if we only have to parse the large database files once at the beginning. 
  - NOTE: these is some redundancy across the four database files (mostly between the MMETSP and OTHER databases I think) so the script will remove proteins with identical names AND sequences if they do occur.
- For each KO number:
  - For each protein annotated with that KO number:
    1. Concatenate the BLAST results from the different databases
    2. Filter *e*-value < **1e-10**
    3. Sort hits by query_id then by bitscore
    4. Down sample subject sequences (with `down_sample_taxa.sh`) using taxonomic information about the different sequences. The script will:
       1. Calculate the number of taxonomic 'Classes' present in the hit subjects.
       2. Calculate the max number of 'Genera' per 'Class' to keep (`numGen=$((180/numClass))`)
       3. If the number of 'Genera' is < 25: retrieve 2 subject sequences per 'Genera'; else retrieve 1 subject sequence per 'Genera'
    5. Get the hits for the down sampled subject sequences
    6. Get the protein sequences for the down sampled subject sequences + the current KR01 (query) protein
    7. Add (concatenate) the down sampled proteins to a file with other proteins from the other KR01 proteins annotated with the same KO number. 
  - Remove redundant sequences from the combined `fasta` file with down sampled sequences from each of the manually corrected (query) proteins. 

```bash
./run_02_process_blast.sh
```

Check log file `run_02_process_blast.sh.log.*` for errors. **None found.**

## 03 - Align proteins and build trees

For each set of KR01 proteins, align the proteins that we just collected (and down sampled), and build a phylogenetic tree. 

The `run_03_build_trees.sh` script:

- Runs `mafft` on the KR01 proteins + down sampled BLASTP hits. 
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

## 04 - Trim alignments

Trim `mafft` alignments using `trimal`. We will use these trimmed/cleaned alignments to rebuild trees for each KO number; should be quicker to run, as the alignments will be smaller and the alignment will hopefully be more informative, giving us similar or better results. These trees should give use similar results to the full alignment trees, if they do not then this could tell us something about how robust the phylogenetic signal is for the sequences.

The `run_04_trim_alignments.sh` script:

- Runs `trimal` using the `-automated1` flag: *Use a heuristic selection of the automatic method based on similarity statistics. (see User Guide). (Optimized for Maximum Likelihood phylogenetic tree reconstruction).*
- Won't run trimming on alignment with <= 3 sequences as there is no point as be can't use them to build trees.
- NOTE:
  - `trimal` will remove any sequences from the alignment that after trimming are composed only of gaps. 
  - `trimal` doesn't like ambiguous amino acids either (such as 'Z'). Will print a bunch of errors but doesn't seem to stop `trimal` from producing a trimmed alignment. The ambiguous nucleotide can still be preserved in the final alignment as well, so it doesn't just ignore/remove that column in the alignment. 

```bash
./run_04_trim_alignments.sh
```

Check file `run_04_trim_alignments.sh.log.*` for unexpected errors. **None found that we didn't expect.** Mostly just warning about removing gap filled sequences after trimming or ambiguous amino acids.

## 05 - Build trees using trimmed alignments

Build maximum likelihood phylogenetic trees using the `trimal` trimmed alignments. 

The `run_05_build_trees.sh` script:

- Runs `iqtree` on the alignment to build a maximum likelihood tree with substitution model picked automatically and bootstrap support added via 2000 ultrafast bootstraps. 
  - If an alignment had <= 3 sequences it couldn't/didn't have a tree built using `iqtree`. 

```bash
./run_05_build_trees.sh
```

Check log file `run_05_build_trees.sh.log.*` for errors. **None found.**

Check `*.build_trees.log` for errors. **None found.** Only one KO number (K03651) had <= 3 sequences and had to be ignored during tree construction.

```bash
cat *.faa.mafft.aln.trimal1.build_trees.log | grep -v '#\|CMD'
```

Check `iqtree` log files for errors. 

```bash
grep -i 'error\|warn' *.faa.mafft.aln.trimal1.log
```

No errors found for `mafft`. Lots of warnings were found for `iqtree` but most are unavoidable or just warning the user about features of the data, so are probably not critical.

## 06 - Prepare for visualization

Prepare information (color and annotations) for each alignment ready for visualization using `TreeViewer`. 

The `run_06_prepare_for_visualization.sh` script:

- Runs `AliStat` to generate diagnostic/QC statistics and plots for each alignment. 
  - The "Cr" statistic from `AliStat` can show us which sequences don't fit well in the alignment (e.g., have very few residues that overlap with other sequences in the alignment) .
  - The "Cc" histogram can also help identify if there are pairs of sequences in the alignment that don't overlap in the alignment (e.g., one sequence aligns to the start of the alignment and the second aligns to the end, the two sequences don't actual overlap in the alignment and so their position in the resulting tree will be inferred based on their similarity to sequences that overlap both). 
- Generate cumulative Cc, and Cc and Cr histograms using the `R` scripts produced by `AliStat`
- Create an annotation file listing each sequence in the alignment. Annotate each sequence with a color (based on which broad domain of life it is from, e.g., Bacteria, Eukaryotes, etc.) and the `AliStat` Cr score (how well that sequence fits/overlaps with other sequences in the alignment).

```bash
./run_06_prepare_for_visualization.sh
```

Check file `run_06_prepare_for_visualization.sh.log.*` for errors. **None found.** K03651 could be run because is had <= 3 sequences in its alignment.

## Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/manually_correct_gene_models_phylogenetic_analysis_FullProtein"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="*.sh*" --include="*.log" --include="*.aln" --include="*.contree" --exclude="*" \
  ${WD}/./ \
   . --dry-run
```

