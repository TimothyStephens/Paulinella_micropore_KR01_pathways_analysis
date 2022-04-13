# Phylogenetic analysis of manually corrected proteins grouped by KO number

Perform phylogenetic analysis using the manually corrected proteins but grouping those sequences annotated with the same KO numbers into the same tree. This analysis is to see if the proteins are of host, endosymbiont, or foreign (HGT) origin, and to see if different proteins annotated with the same KO number have different evolutionary histories. 

Assumes that analysis is `manually_correct_gene_models_phylogenetic_analysis_FullProtein` has been completed.

## Setup analysis directory

Link manually corrected proteins.

```bash
ln -s ../manually_correct_gene_models/combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.pep 
```

Link to custom scripts for text parsing and for filtering `fasta` files.

```bash
ln -s ../../../02_Scripts/grepf_column.py
ln -s ../../../02_Scripts/grepf_fasta.py
ln -s ../../../02_Scripts/add_value_to_table.py
```

Link to BLAST results from previous phylogenetic analysis where we kept the proteins separate (so we don't have to redo the slow BLAST step).

```bash
for F in ../manually_correct_gene_models_phylogenetic_analysis_FullProtein/MSTRG.*.pep.faa.blastp_{BACTERIA,METAZOA,MMETSP,OTHER}.outfmt6; 
do 
	ln -s $F
done
```

Setup bash environment.

```bash
conda activate py27
```

## Link proteins to KO numbers

We need a file that links the manually corrected proteins with their annotated KO numbers. This file in our case is called `SuppTable_S2.transcript_info.txt`, and was created by copying the information for the Master Excel spreadsheet. So this file not only can be used to link protein_id<->KO_number but can also be used to format results in the same order as the Excel document, so it makes it easier to copy-paste results into the master spreadsheet. 

## 02 - Process BLAST results

Don't have to do step 01 as we already have the BLAST results from when we built trees separately for each manually corrected protein; results in `manually_correct_gene_models_phylogenetic_analysis_FullProtein`. 

For each manually corrected protein, combine the BLAST results from each database, filter the hits, down sample the subject sequences using taxonomic information, and extract the proteins from the final set of selected sequences (to be used for phylogenetic analysis). 

The `run_02_process_blast.sh` script:

- Extracts all subject sequences (with hits) across all query proteins + add manually corrected proteins
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
        6. Get the protein sequences for the down sampled subject sequences + the current manually corrected (query) protein
        7. Add (concatenate) the down sampled proteins to a file with other proteins from the other manually corrected proteins annotated with the same KO number. 
    - Remove redundant sequences from the combined `fasta` file with down sampled sequences from each of the manually corrected (query) proteins. 

```bash
./run_02_process_blast.sh
```

Check log file `run_02_process_blast.sh.log.*` for errors.

A one query protein (MSTRG.15268.1.p1) had no hits that passed the *e*-value filter; most hits had weak (>0.001) e-values. This caused a 'division by 0' errors in the log file but that OK. The rest of the proteins annotated with the KO number had hits.

> \>MSTRG.15268.1.p1
> MMELVDCVRTAALFEYGSQHPLILALPLRACTCPCSCTPIPPTPHPASPNKQPPAFYPISPQKTPCDHISVISAGKSGQFPANLLAGPPAPRRPGWYCCRCVGMVGVRETLWWLWVLRKYGTTHVVIADGGECGHAPDDVGTGGALHVAKQLQKELGVELLVLPPLLYLPDEDSYKASHNIAPGERCEQSCCSGLVSKQIRSVLARTFGGQKTTFVTAQLKNKEQKKEQGPCKQQNRPKKKGRPHVEEKMEEKKKKKKKEKKEEKNEEKKETKREKEKREKMDKKEKKEKKEEEKSEQLEPDLVVLKKASSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSPFSSPSATNISAVHTTADRGPYLSTLKAWTFPHVVSVLQKHSYPTSLGVSVHASSSEAELSYFADHASGLAAISSNSNGNGGRRIVTVGGVWGDRLPNRPHRTPYSVGQVLERYGFEALWLVGLFATGYGIFALGRFGTRWVRLRLTVLGRTYYAKYWPGYSVDVDWLTGVVHSHWRTLRAVSPQVVTDTLHNGSQALADVLRAVWGGTLSQISQLAASTLPALLASAAVAASSSSPSDADAAASALSSPPLLLRPPSPPSPPSPPSSSSSSSSSSSSSPCSSSSSSAPPSFSPSFSSSSSSSSSSSSSCCSSSSSSSSSSSSSSSSSSSSSSS*

MSTRG.15268.1.p1 only had 18% UniProt KEGG Ortholog coverage, other protein annotated with this KO number has 98.5% coverage, so this protein is probably a false positive or highly diverged.

## 03 - Align proteins and build trees

For each manually corrected protein, align the proteins that we just collected (and down sampled), and build a phylogenetic tree. 

The `run_03_build_trees.sh` script:

- Runs `mafft` on the KR01 proteins + down sampled BLASTP hits. 
- Runs `iqtree` on the alignment to build a maximum likelihood tree with substitution model picked automatically and bootstrap support added via 2000 ultrafast bootstraps. 
  - If an alignment had <= 3 sequences it couldn't/didn't have a tree built using `iqtree`. 

```bash
./run_03_build_trees.sh
```

Check log file `run_03_build_trees.sh.log.*` for errors. **None found.**

Check `*.build_trees.log` for errors. **None found.** Only one KO number (K03651) had <= 3 sequences and had to be ignored during tree construction.

```bash
cat *.faa.build_trees.log | grep -v '#\|CMD'
```

Check `mafft` and `iqtree` log files for errors. 

```bash
grep -i 'error\|warn' *.faa.mafft.log
grep -i 'error\|warn' *.faa.mafft.aln.log
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

## Clean '/' characters from sequence names

Some sequences in the taxonomically broad database have `/` characters in their names. `iqtree` replaces these characters with `_` , we need to do the same to the alignment files and annotation files otherwise we can't link these sequences across the input files for when we visualize them. 

```bash
grep -l '/' *.mafft.aln *.mafft.aln.annots.txt *.mafft.aln.trimal1 *.mafft.aln.trimal1.annots.txt  | sort
```

The affected datasets are: K00876, K01520, K01933, K02335, K10807, K10808.

```bash
for F in K00876 K01520 K01933 K02335 K10807 K10808; 
do
	sed -i -e 's@/@_@g' $F.*.mafft.aln $F.*.mafft.aln.annots.txt $F.*.mafft.aln.trimal1 $F.*.mafft.aln.trimal1.annots.txt
done
```

## Download results from server

Download results files + alignments and trees for visualization using `TreeViewer`.

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/manually_correct_gene_models_phylogenetic_analysis_FullProtein_Grouped"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="SuppTable_S2.transcript_info.txt*" --include="*.sh*" --include="*.log" --include="*.aln" --include="*.trimal1" --include="*.contree" --include="*.annots.txt" --include="*.pdf" --include="*.AliStat.Summary.txt" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

## Visualize trees using `TreeViewer`

Run `TreeViewer` on **FULL** alignments. Ignore K03651 which had too few sequences for phylogenetic analysis.

`TreeViewer` needs a few extra template file to work. They are:

- TreeView.AliStat_Cr_Legend.svg
  - `SVG` image to use as the legend showing the color range used to plot the `AliStat` Cr values for each sequence in the tree.
- TreeView.Legend.Histogram_Cc.txt
  - `Markdown` text used to help insert the `AliStat` Cc histogram into the tree.
- TreeView.Legend.Histogram_Cr.txt
  - `Markdown` text used to help insert the `AliStat` Cr histogram into the tree.
- TreeView.Legend.txt
  - Legend detailing the colors used for the different taxa + the highlighted nodes + Cr value color scale.
- TreeView.Node_Shape.txt
  - Code to control the placement of a black circle on nodes in the tree with bootstrap support >= 95%
- TreeView.Node_states_AliStat.txt
  - Code to setup the color gradient used then plotting the `AliStat` Cr values next to the tips in the tree.
- TreeView.Plot_alignment.txt
  - Code to set the color of the alignment the same as the tips in the tree (i.e., color the alignment by taxonomic group)
- TreeView.plot.txt
  - Main template with the `TreeView` commands to run to build the trees.

```bash
while read KO; do
ALN="${KO}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln"
TREE="${ALN}.contree"

if [ ! -f "$ALN" ]; then continue; fi #Ignore missing K03651

sed -e "s/:TREE:/${TREE}/" \
    -e "s/:ANNOT:/${ALN}.annots.txt/" \
    -e "s/:ALN:/${ALN}/" \
    -e "s/:ALISTAT_HISTOGRAM_CR_FILE:/${ALN}.AliStat.Histogram_Cr.pdf/" \
    -e "s/:ALISTAT_HISTOGRAM_CC_FILE:/${ALN}.AliStat.Histogram_Cc.pdf/" \
    -e '/^#.*$/d' \
    TreeView.plot.txt  | TreeViewerCommandLine 1> "${TREE}.plot.txt.log" 2>&1
done < <(awk -F'\t' '$1~"MSTRG" {print $2}' SuppTable_S2.transcript_info.txt | sort | uniq)
```

Run `TreeViewer` on **TRIMMED** alignments. Ignore K03651 which had too few sequences for phylogenetic analysis.

```bash
while read KO; do
ALN="${KO}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.trimal1"
TREE="${ALN}.contree"

if [ ! -f "$ALN" ]; then continue; fi #Ignore missing K03651

sed -e "s/:TREE:/${TREE}/" \
    -e "s/:ANNOT:/${ALN}.annots.txt/" \
    -e "s/:ALN:/${ALN}/" \
    -e "s/:ALISTAT_HISTOGRAM_CR_FILE:/${ALN}.AliStat.Histogram_Cr.pdf/" \
    -e "s/:ALISTAT_HISTOGRAM_CC_FILE:/${ALN}.AliStat.Histogram_Cc.pdf/" \
    -e '/^#.*$/d' \
    TreeView.plot.txt  | TreeViewerCommandLine 1> "${TREE}.plot.txt.log" 2>&1
done < <(awk -F'\t' '$1~"MSTRG" {print $2}' SuppTable_S2.transcript_info.txt | sort | uniq)
```

Check for any errors in the `TreeViewer` log files.

```bash
grep 'Unknown command' *.plot.txt.log
grep 'The specified file does not exist' *.plot.txt.log
```

Merge the individual `PDF` files into a single multi-page `PDF` so its easier to read.

Merge **FULL** alignment results.

```bash
./merge_PDFs.py -i $(awk -F'\t' '$1~"MSTRG" && $2!="K03651" {print $2".blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.contree.pdf"}' SuppTable_S2.transcript_info.txt | awk '!seen[$1]++ {printf "%s ", $1}') -o Trees.merged.pdf --debug
```

Merge **TRIMMED** alignment results.

```bash
./merge_PDFs.py -i $(awk -F'\t' '$1~"MSTRG" && $2!="K03651" {print $2".blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.trimal1.contree.pdf"}' SuppTable_S2.transcript_info.txt | awk '!seen[$1]++ {printf "%s ", $1}') -o Trees.trimal1.merged.pdf --debug
```

