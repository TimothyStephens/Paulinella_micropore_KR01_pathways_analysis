# Perform phylogenetic analysis of manually corrected genes

Perform phylogenetic analysis of manually corrected genes with orthologs from other *Paulinella* species added; *Paulinella* orthologs are selected using results from KAAS annotation of all predicted proteins from each species. 

## Setup analysis directory

Link to *P. chromatophora* and *P. ovalis* predicted proteins and KAAS annotations.

```bash
ln -s ../../../../../genome_data/Paulinella_chromatophora_CCAC0185/nuclear_genome/databases/Paulinella_chromatophora_CCAC0185_nuclear.pep.faa
ln -s ../../../../../genome_data/Paulinella_chromatophora_CCAC0185/nuclear_genome/databases/Paulinella_chromatophora_CCAC0185_nuclear.pep.KAAS_Combined.ko
ln -s ../../../../../genome_data/Paulinella_ovalis_6_SingleCells/nuclear_genome/databases/Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep.faa
ln -s ../../../../../genome_data/Paulinella_ovalis_6_SingleCells/nuclear_genome/databases/Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep.KAAS_Combined.ko
```

Link to custom scripts for parsing `fasta` and text files.

```bash
ln -s ../../../02_Scripts/grepf_fasta.py
ln -s ../../../02_Scripts/grepf_column.py
```

Setup bash environment.

```bash
conda activate py27
```

## 01 - Link to BLAST results

Link to BLAST results that we generated previously.

```bash
awk -F'\t' '$1~"^MSTRG" {print $1}' SuppTable_S2.transcript_info.txt \
  | while read F; do 
  ln -s ../manually_correct_gene_models_phylogenetic_analysis_FullProtein/$F.pep.faa.blastp_BACTERIA.outfmt6;
  ln -s ../manually_correct_gene_models_phylogenetic_analysis_FullProtein/$F.pep.faa.blastp_METAZOA.outfmt6;
  ln -s ../manually_correct_gene_models_phylogenetic_analysis_FullProtein/$F.pep.faa.blastp_MMETSP.outfmt6;
  ln -s ../manually_correct_gene_models_phylogenetic_analysis_FullProtein/$F.pep.faa.blastp_OTHER.outfmt6;
done
```

## 02 - Process BLAST results

Extract *P. chromatophora* and *P. ovalis* proteins annotated with target KO numbers.

```bash
PC="Paulinella_chromatophora_CCAC0185_nuclear.pep"
PO="Paulinella_ovalis_nuclear_SingleCell_ALL.assembly.fasta.augustus.pep"

awk -F'\t' '$1~"^MSTRG" {print $2}' SuppTable_S2.transcript_info.txt \
  | while read KO; 
  do
    ./grepf_fasta.py -i "$PC.faa" -f <(awk -vKO="$KO" '$2==KO {print $1}' "$PC.KAAS_Combined.ko" | sort | uniq) \
      > "$KO.Other_Paulinella_seqs.pep.faa"
    ./grepf_fasta.py -i "$PO.faa" -f <(awk -vKO="$KO" '$2==KO {print $1}' "$PO.KAAS_Combined.ko" | sort | uniq) \
      >> "$KO.Other_Paulinella_seqs.pep.faa"
  done
```

Run `run_02_process_blast.sh` script to filter (*e*-value < 1e-10) and down sample the BLAST results for all the KR01 proteins associated with each KO number.

```bash
./run_02_process_blast.sh
```

## 03 - Align proteins

Align proteins in each KO set using `mafft`

```bash
./run_03_build_alignment.sh
```

## 04 - Build trees - Full alignments

Build phylogenetic trees using `iqtree`; 

```bash
./run_04_build_trees.sh
```

## 05 - Trim alignment

Trim full protein alignments (from `mafft`) using `trimal`

```bash
./run_05_trim_alignments.sh
```

## 06 - Build trees - Trimmed alignments

Build phylogenetic trees using `iqtree`; let `iqtree` pick the best evolutionary model for each trees.

```bash
./run_06_build_trees.sh
```

## 07 - Visualize trees

Prepare phylogenetic trees (both full and trimmed) for visualization using `TreeViewer`; runs `AliStat` to generate metrics about each alignment (can help diagnose problems downstream).

```bash
./run_07_prepare_for_visualization.sh
```

**Clean '/' characters from sequence names**

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

## Download results

Download results files + alignments and trees for visualization using `TreeViewer`.

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2022-01-13/manually_correct_gene_models_phylogenetic_analysis_FullProtein_Grouped_DiffParams3"
rsync -zarvL --prune-empty-dirs --relative \
 --include="*/" --include="SuppTable_S2.transcript_info.txt*" --include="*.sh*" --include="*.log" --include="*.aln" --include="*.trimal1" --include="*.contree" --include="*.annots.txt" --include="*.pdf" --include="*.AliStat.Summary.txt" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

## Visualize trees using `TreeViewer`

Run `TreeViewer` on **FULL** and **TRIMMED** alignments. Ignore K03651 which had too few sequences for phylogenetic analysis.

`TreeViewer` needs a few extra template file to work. They are:

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

sed -e "s/:TREE:/${TREE}/" \
    -e "s/:ANNOT:/${ALN}.annots.txt/" \
    -e "s/:ALN:/${ALN}/" \
    -e '/^#.*$/d' \
    TreeView.plot.txt  | TreeViewerCommandLine 1> "${TREE}.plot.txt.log" 2>&1
done < <(awk -F'\t' '$1~"MSTRG" {print $2}' SuppTable_S2.transcript_info.txt | sort | uniq)
```

Run `TreeViewer` on **TRIMMED** alignments. Ignore K03651 which had too few sequences for phylogenetic analysis.

```bash
while read KO; do
ALN="${KO}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.trimal1"
TREE="${ALN}.contree"

sed -e "s/:TREE:/${TREE}/" \
    -e "s/:ANNOT:/${ALN}.annots.txt/" \
    -e "s/:ALN:/${ALN}/" \
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
./merge_PDFs.py -i $(awk -F'\t' '$1~"MSTRG" {print $2".blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.contree.pdf"}' SuppTable_S2.transcript_info.txt | awk '!seen[$1]++ {printf "%s ", $1}') -o Trees.merged.pdf --debug
```

Merge **TRIMMED** alignment results.

```bash
./merge_PDFs.py -i $(awk -F'\t' '$1~"MSTRG" {print $2".blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.trimal1.contree.pdf"}' SuppTable_S2.transcript_info.txt | awk '!seen[$1]++ {printf "%s ", $1}') -o Trees.trimal1.merged.pdf --debug
```

