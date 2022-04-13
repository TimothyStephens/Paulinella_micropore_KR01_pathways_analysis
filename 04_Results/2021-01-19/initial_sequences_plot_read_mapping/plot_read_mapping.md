# Plot read mapping along target gene CDS

Visualize read mapping along target gene CDS. Use to identify CDS with putative mis-assembled/mis-predicted regions.

Assumes analysis in `initial_sequences_info/` and `initial_sequences_blastp_UniProt_KEGG_seqs/` has been completed.

## Setup analysis directory

Link files with KR01 gene IDs and KO numbers + extra info we collected for each protein.

```bash
ln -s ../initial_sequences_info/DNA_Replication.txt.extra_info_added.txt
ln -s ../initial_sequences_info/Purine.txt.extra_info_added.txt
ln -s ../initial_sequences_info/Pyrimidine.txt.extra_info_added.txt
```

Link to BAM file with RNA-seq reads locally aligned against the KR01 CDS and the BAM file with RNA-seq reads aligned (splice-aware) against the KR01 genome.

```bash
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/All_combined.HISAT2_RNAseq_mapping.coordSorted.bam
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/All_combined.HISAT2_RNAseq_mapping.coordSorted.bam.bai
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/PredGenes_CDS_RNAseq_bowtie2_local/All_combined.local_mapping.coordsorted.bam
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/PredGenes_CDS_RNAseq_bowtie2_local/All_combined.local_mapping.coordsorted.bam.bai
```

Link to the KR01 (1) genome assembly, (2) the CDS sequences, and (3) the gene models gff3 file. 

```bash
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.assembly.fasta
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.fai
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3

ln -s Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fa
```

Link to KR01 protein sequences (used for crTP analysis) + crTP alignment that we will use as the basis for our search. 

```bash
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/crTP_containing_proteins/KR01_crTP_StartingSet_crTPonly_ManualCorr.aln
```

Setup bash environment.

```bash
conda activate py27

export HHLIB="/home/timothy/programs/hhsuite-3.3.0-SSE2-Linux"
export PERL5LIB="$HHLIB/scripts"
export PATH="$PATH:$HHLIB/bin:$HHLIB/scripts"
export PATH="$PATH:/home/timothy/programs/hmmer-3.1b2/bin"

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
```

## Extract info for plotting

### Extract genome coords of target genes

Used for defining the bounds of the read mapping plots that we will generate.

```bash
function get_coords() {
        INFO="$1"
        python ~/scripts/add_value_to_table.py \
         -i <(awk -F'\t' '$2~"^g"{print "Paulinella_micropora_KR01_nuclear___"$2}' $INFO.txt.extra_info_added.txt) \
         -a <(awk '{print $1"\t0\t"$2"\t"$1}' Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna.fai) \
         > $INFO.CDS_coords.bed

        awk -F'\t' '$2~"^g"{print "Paulinella_micropora_KR01_nuclear___"$6"\t"$7-1"\t"$8"\tPaulinella_micropora_KR01_nuclear___"$2}' $INFO.txt.extra_info_added.txt > $INFO.GENOME_coords.bed
}

get_coords "DNA_Replication"
get_coords "Purine"
get_coords "Pyrimidine"
```

### Get crTP coords from HMM search

Run `hmmsearch` using the crTP model built for *P. micropora* KR01, given to us by Eva Nowack.

```bash
A=KR01_crTP_StartingSet_crTPonly_ManualCorr.aln
## Build a hmmer *.hmm file from the manualy curated alignment. 
reformat.pl clu fas $A $A.fasta
hmmbuild $A.fasta.hmm $A.fasta

P="Paulinella_micropora_KR01_nuclear.augustus_190918v4.pep.faa"
## Search for crTP in proteins using HMM
hmmsearch --domtblout $P.domtblout $A.fasta.hmm $P
```

Extract coords for crTP from a HMM search I completed previously. Filter >30% HMM length coverage. 

Need to multiple coords by x3 to get CDS coords not protein coords (assumes that partial CDS don't have any nucleotides from partial codons at the ends of the sequences).

```bash
grep -v '^#' "$P.domtblout" \
  | awk '$12<1e-5 && ($17-$16)+1 > $6*0.3 {print $1"\t"($18-1)*3"\t"$19*3"\tcrTP"}' \
  > crTP.filtered_30pct_HMM_len.CDS_coords.bed
```

### Combine bed features to plot (crTP locations + KEGG BLASTP top hit)

```bash
cat crTP.filtered_30pct_HMM_len.CDS_coords.bed ../initial_sequences_blastp_UniProt_KEGG_seqs/KEGG_and_GENE_id_pairs.txt.blastp_KEGG_UniProt_tophit.CDS_coords.bed \
  | bedtools sort > Combined.CDS_coords.bed
```

### Run script to plot read mapping

The `Gviz_plot_read_mapping_CDS.R` and `Gviz_plot_read_mapping_GENOME.R` R scripts use the `Gvis` package to plot the mapped reads along the target sequences (either CDS or a region of the genome).

```bash
conda activate r4_env

F="DNA_Replication.CDS_coords.bed"
Rscript Gviz_plot_read_mapping_CDS.R "$F" "$F.read_mapping_plots.pdf" \
  > "$F.read_mapping_plots.log" 2>&1
F="Purine.CDS_coords.bed"
Rscript Gviz_plot_read_mapping_CDS.R "$F" "$F.read_mapping_plots.pdf" \
  > "$F.read_mapping_plots.log" 2>&1
F="Pyrimidine.CDS_coords.bed"
Rscript Gviz_plot_read_mapping_CDS.R "$F" "$F.read_mapping_plots.pdf" \
  > "$F.read_mapping_plots.log" 2>&1

F="DNA_Replication.GENOME_coords.bed"
Rscript Gviz_plot_read_mapping_GENOME.R "$F" "$F.read_mapping_plots.pdf" \
  > "$F.read_mapping_plots.log" 2>&1
F="Purine.GENOME_coords.bed"
Rscript Gviz_plot_read_mapping_GENOME.R "$F" "$F.read_mapping_plots.pdf" \
  > "$F.read_mapping_plots.log" 2>&1
F="Pyrimidine.GENOME_coords.bed"
Rscript Gviz_plot_read_mapping_GENOME.R "$F" "$F.read_mapping_plots.pdf" \
  > "$F.read_mapping_plots.log" 2>&1
```

## Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/initial_sequences_plot_read_mapping/"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="*.R" --include="*.log" --include="*.pdf" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

