# Map reads against KR01 + manually corrected CDS using `salmon`  

Align RNA-seq reads from control + high light samples against KR01 CDS (with initial sequences removed) + manually corrected CDs added. Can use these results to infer the diurnal pattern of the genes annotated to the pathways of interest. 

## Setup analysis directory

Link to fast files with KR01 CDS (without initial genes) + KR01 transcripts.

```bash
ln -s ../manually_correct_gene_models_read_mapping/Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt
```

Link to `samples.txt` file that lists the RNA-seq datasets that we want to align.

```bash
ln -s ../manually_correct_gene_models_read_mapping/samples.txt
```

Link to custom scripts for post processing of read count matrices. 

```bash
ln -s ../../../02_Scripts/add_value_to_table.py
ln -s ../../../02_Scripts/grepf_column.py
ln -s ../../../02_Scripts/normalize_matrix.R 
ln -s ../../../02_Scripts/plot_time_course_expression_data.R
ln -s ../../../02_Scripts/melt_table.py
```

Setup bash environment.

```bash
conda activate py27
```

## Salmon index

Index fasta file with KR01 CDS + manually corrected CDS

```bash
./run_salmon_index.sh 
```

Check `run_salmon_index.sh.log.*` for errors. **None found**.

## Salmon quant

Align reads from each sample against the CDS file using `salmon quant`

```bash
./run_salmon_quant.sh
```

Check `run_salmon_quant.sh.log.*` for errors. **None found**.

## Combine and post-process count results

Merge results from each sample into one big counts table.

```bash
export PATH="$PATH:/home/timothy/programs/salmon-1.1.0_linux_x86_64/bin"

## Merge salmon results.
DIR_LIST=$(awk -F'\t' '{printf " %s", $2}' samples.txt)
F="Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt"

salmon quantmerge --column numreads --output "$F.salmon.allSamples.numreads.matrix" --quants $DIR_LIST
salmon quantmerge --column tpm --output "$F.salmon.allSamples.tpm.matrix" --quants $DIR_LIST
```

Get read mapping rate for each sample. Can be used downstream to highlight problematic samples. 

```bash
echo -e "sample_id\tSalmon_mapping_rate" > stats.txt
for D in $(awk -F'\t' '{printf " %s", $2}' samples.txt); do echo -e "$D\t"$(grep 'Mapping rate' $D/logs/salmon_quant.log | awk '{print $8}'); done >> stats.txt
```

Compress and remove each salmon sample directory

```bash
for D in $(awk -F'\t' '{printf " %s", $2}' samples.txt); do echo "$D"; tar -zcf "$D.tar.gz" "$D" && rm -r "$D"; done
```

## Plot

Rename and Normalize counts using median-of-ratios

```bash
F="Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt"

# Have to remove '-' in names as R does not like them. Change to '_'
cat $F.salmon.allSamples.numreads.matrix \
 | awk -F'\t' 'BEGIN{OFS="\t"} { if(NR==1) {for (i=2; i<=NF; i++) {split($i,a,"-"); $i=a[1]"_"a[2]"_"a[3]}; print $0 } else {print $0} }' \
 > $F.salmon.allSamples.numreads.matrix.renamed

conda activate r4_env

# median of ratios
Rscript normalize_matrix.R $F.salmon.allSamples.numreads.matrix.renamed
# [1] "max SizeFactor: 1.39548182369577"
# [1] "min SizeFactor: 0.600952578222596"
```

List genes to plot in `genes2plot.txt`. The enzyme reactions 1.17.4.2, 2.7.4.6, and 3.1.3.5 are repent in both the Purine and Pyrimidine pathways, so remove the second entry for these IDs to make things simpler.

Example text file:

| EC Number | GeneID           |
| --------- | ---------------- |
| 6.5.1.2   | MSTRG.15911.1.p1 |
| 6.5.1.2   | MSTRG.8876.1.p1  |
| 6.5.1.2   | MSTRG.7660.1.p1  |
| 2.7.7.7   | MSTRG.17738.1.p1 |
| 2.7.7.7   | MSTRG.27742.1.p1 |

Plot expression patterns over time points using just the "Control" (CL) samples using **median-of-ratios** normalized data.

```bash
## Split CL and HL into diff files and rename sample names (make TP numbers)
cut -f 1,5-25 $F.salmon.allSamples.numreads.matrix.renamed.normalized_counts.txt \
| awk -F'\t' 'BEGIN{OFS="\t"} { if(NR==1) {for (i=2; i<=NF; i++) {split($i,a,"_"); $i=substr(a[2], 1, length(a[2])-1)}; print $0 } else {print $0} }' \
> $F.salmon.allSamples.numreads.matrix.renamed.normalized_counts.txt.CL

## Melt CL results
python melt_table.py \
  -i $F.salmon.allSamples.numreads.matrix.renamed.normalized_counts.txt.CL \
  -o $F.salmon.allSamples.numreads.matrix.renamed.normalized_counts.txt.CL.melted

## Grep just the genes associated with the EC numbers that we are after.
## Also reformat table so its EC \t Time \t Count \t GeneID (format required so plotting script displays all genes from same EC in a plot; also allows for same gene to be associated with multiple EC's)
echo -e "Name\tTime\tExpression\tCondition" > genes2plot.txt.norm_count_matrix.melted
while read line; 
do 
  EC=$(echo "$line" | cut -f1); 
  GENE=$(echo "$line" | cut -f2); 
  ./grepf_column.py -f <(echo "$GENE" | sed -e 's/\.p.*//') \
      -i "$F.salmon.allSamples.numreads.matrix.renamed.normalized_counts.txt.CL.melted" \
    | awk -F'\t' -vEC="$EC" '{print EC"\t"$2"\t"$3"\t"$1}' \
    >> genes2plot.txt.norm_count_matrix.melted
done < genes2plot.txt

## Plot expression patterns - both single plot per file and 4x4 combined
conda activate r4_env
Rscript plot_time_course_expression_data.R \
  genes2plot.txt.norm_count_matrix.melted \
  genes2plot.txt.norm_count_matrix.melted.plots.pdf \
  <(cut -f1 genes2plot.txt | uniq)
```

Plot expression patterns over time points using just the "Control" (CL) samples using **Transcripts Per Million (TPM)** normalized data.

```bash
## Split CL and HL into diff files and rename sample names (make TP numbers)
cut -f 1,5-25 $F.salmon.allSamples.tpm.matrix \
| awk -F'\t' 'BEGIN{OFS="\t"} { if(NR==1) {for (i=2; i<=NF; i++) {split($i,a,"-"); $i=substr(a[2], 1, length(a[2])-1)}; print $0 } else {print $0} }' \
> $F.salmon.allSamples.tpm.matrix.CL

## Melt CL results
python melt_table.py \
  -i $F.salmon.allSamples.tpm.matrix.CL \
  -o $F.salmon.allSamples.tpm.matrix.CL.melted

## Grep just the genes associated with the EC numbers that we are after.
## Also reformat table so its EC \t Time \t Count \t GeneID (format required so plotting script displays all genes from same EC in a plot; also allows for same gene to be associated with multiple EC's)
echo -e "Name\tTime\tExpression\tCondition" > genes2plot.txt.tpm_matrix.melted
while read line; 
do 
  EC=$(echo "$line" | cut -f1); 
  GENE=$(echo "$line" | cut -f2); 
  ./grepf_column.py -f <(echo "$GENE" | sed -e 's/\.p.*//') \
      -i "$F.salmon.allSamples.tpm.matrix.CL.melted" \
    | awk -F'\t' -vEC="$EC" '{print EC"\t"$2"\t"$3"\t"$1}' \
    >> genes2plot.txt.tpm_matrix.melted
done < genes2plot.txt

## Plot expression patterns - both single plot per file and 4x4 combined
conda activate r4_env
Rscript plot_time_course_expression_data.R \
  genes2plot.txt.tpm_matrix.melted \
  genes2plot.txt.tpm_matrix.melted.plots.pdf \
  <(cut -f1 genes2plot.txt | uniq)
```

## Just the DNA Replication proteins

Plot just the DNA Replication genes that are nuclear-encoded. Groups them by KO number rather than EC number as different complex subunits have the same EC number and we would want to keep them separate. 

```bash
## Grep just the genes associated with the KO numbers that we are after.
## Also reformat table so its EC \t Time \t Count \t GeneID (format required so plotting script displays all genes from same EC in a plot; also allows for same gene to be associated with multiple EC's)
echo -e "Name\tTime\tExpression\tCondition" > genes2plot_DNArep.txt.tpm_matrix.melted
while read line; 
do 
  EC=$(echo "$line" | cut -f1); 
  GENE=$(echo "$line" | cut -f2); 
  ./grepf_column.py -f <(echo "$GENE" | sed -e 's/\.p.*//') \
      -i "$F.salmon.allSamples.tpm.matrix.CL.melted" \
    | awk -F'\t' -vEC="$EC" '{print EC"\t"$2"\t"$3"\t"$1}' \
    >> genes2plot_DNArep.txt.tpm_matrix.melted
done < genes2plot_DNArep.txt

## Plot expression patterns - both single plot per file and 4x4 combined
conda activate r4_env
Rscript plot_time_course_expression_data.R \
  genes2plot_DNArep.txt.tpm_matrix.melted \
  genes2plot_DNArep.txt.tpm_matrix.melted.plots.pdf \
  <(cut -f1 genes2plot_DNArep.txt | uniq)
```





## Results

### Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/manually_correct_gene_models_Salmon_mapping"
rsync -zarvL --prune-empty-dirs --relative \
 --include="*/" --include="*.sh*" --include="*.py" --include="*.R" --include="*.pdf" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```
