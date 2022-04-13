# Collect extra info about identified proteins

Collect info about KR01 KEGG proteins. Info such as length, genome coords, and number of exons can help with interpreting results from downstream analysis. 

## Setup analysis directory

Link to files that list KR01 gene IDs and the KO numbers that they are annotated to.

```bash
ln -s ../../../01_Data/2022-01-13/Histidine_subset.txt
```

Setup bash environment.

```bash
conda activate py27
```

## Collect and prepare extra info from different data sources

KEGG ID description/info.

```bash
KEGG_INFO="/scratch/timothy/databases/KEGG/ko00001.keg"
awk '$1=="D"' $KEGG_INFO | sed -e 's/D      //' -e 's/  /\t/' > info.KEGG_descriptions
```

Protein length. Extract from `*.fai` files for Nuclear and Plastid/Chromatophore proteins. 

```bash
NUCLEAR_PEP_FAI="/scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa.fai"
PLASTID_PEP_FAI="/scratch/timothy/genome_data/Paulinella_micropora_KR01/plastid_genome/databases/Paulinella_micropora_KR01_plastid.pep.faa.fai"
cat "$PLASTID_PEP_FAI" "$NUCLEAR_PEP_FAI" | cut -f1,2 | sed -e 's/.*___//' -e 's/^lcl_//' > info.protein_lengths
```

No. exons; only need nuclear counts, plastid genes are all single exon. Extract from a file I prepared previously; this file was created using the KR01 predicted genes `gff3` file.

```bash
NUCLEAR_EXON_COUNTS="/scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3.exon_counts"
sed -e 's/.*___//' "$NUCLEAR_EXON_COUNTS" > info.nuclear_exon_counts
```

Gene reference sequence coords. Need to parse the plastid protein IDs first and reformat them so they match that of the `gff3` file. 

Need to convert protein ID

From: 

> Paulinella_micropora_KR01_plastid___lcl_KX897545.1_prot_APP87810.1_1
>
> Paulinella_micropora_KR01_plastid___lcl_KX897545.1_prot_APP87811.1_2
>
> Paulinella_micropora_KR01_plastid___lcl_KX897545.1_prot_APP87812.1_3
>
> Paulinella_micropora_KR01_plastid___lcl_KX897545.1_prot_APP87813.1_4
>
> ...

To: 

>KX897545.1_prot_APP87810.1_1 <tab> APP87810.1
>
>KX897545.1_prot_APP87811.1_2 <tab> APP87811.1
>
>KX897545.1_prot_APP87812.1_3 <tab> APP87812.1
>
>KX897545.1_prot_APP87813.1_4 <tab> APP87813.1
>
>...

```bash
PLASTID_PEP_FAI="/scratch/timothy/genome_data/Paulinella_micropora_KR01/plastid_genome/databases/Paulinella_micropora_KR01_plastid.pep.faa.fai"

awk -F'\t' '{print $1"\t"$1}' "$PLASTID_PEP_FAI" \
  | sed -e 's/[^\t]*_prot_//' -e 's/_[^\t]*//' -e 's/\t[^\t]*___/\t/' -e 's/lcl_//' \
  | awk '{print $2"\t"$1}' \
  > info.plastid_gff3_seq_names_2_protein_seq_names


NUCLEAR_GFF3="/scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3"
PLASTID_GFF3="/scratch/timothy/genome_data/Paulinella_micropora_KR01/plastid_genome/Paulinella_micropora_KR01_plastid.gff3"

cat \
 <(cat "$NUCLEAR_GFF3" | awk -F'\t' '$3=="mRNA" {print $1"\t"$4"\t"$5"\t"$9}' \
    | sed -e 's/^[^\t]*___//' -e 's/ID=[^;]*___//' -e 's/;.*//' \
    | awk -F'\t' '{print $4"\t"$1"\t"$2"\t"$3}') \
 <(cat "$PLASTID_GFF3" | awk -F'\t' '$3=="CDS"  {print $1"\t"$4"\t"$5"\t"$9}' \
    | sed -e 's/[^\t]*;Name=//' -e 's/;.*//' \
    | ~/scripts/add_value_to_table.py -c 4 \
        -a <(awk '{print $2"\t"$1}' info.plastid_gff3_seq_names_2_protein_seq_names) \
    | awk -F'\t' '{print $5"\t"$1"\t"$2"\t"$3}') \
 > info.ref_genome_coords
```

Add all the info that we just collected to each of the pathway text files (with KEGG_ID <tab> gene_id).

>1. KEGG_ID
>2. gene_id
>3. KEGG description/info
>4. Protein length
>5. No. exons (if not in '-a list' [genes from plastid] assume 1 exon)
>6. Ref seq genome seq_is
>7. Ref seq genome start
>8. Ref seq genome stop

```bash
function add_info() {
        IDS="$1"
        cat "$IDS" \
         | ~/scripts/add_value_to_table.py -c 1 -a info.KEGG_descriptions 2> /dev/null \
         | ~/scripts/add_value_to_table.py -c 2 -a info.protein_lengths \
         | ~/scripts/add_value_to_table.py -c 2 -d 1 -a info.nuclear_exon_counts \
         | ~/scripts/add_value_to_table.py -c 2 -a info.ref_genome_coords \
         > "$IDS.extra_info_added.txt"
}

add_info "Histidine_subset.txt"
```

## Results

### Download results from server

```bash
scp timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2022-01-13/initial_sequences_info/*.extra_info_added.txt .
```

