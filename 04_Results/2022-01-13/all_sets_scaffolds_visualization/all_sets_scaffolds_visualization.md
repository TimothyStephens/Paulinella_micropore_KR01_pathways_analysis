# Subset data for visualization

Get just the features/reads from each file that are associated with the target scaffolds (i.e., the scaffolds with our target genes). Getting just the features on the scaffolds of interest means that we don't have to keep handling the full datasets which (in the case of the aligned reads) can be quite large.

## Setup analysis directory

Link to custom scripts for parsing text files.

```bash
ln -s ../../../02_Scripts/grepf_column.py
ln -s ../../../02_Scripts/add_value_to_table.py
ln -s ../../../02_Scripts/groupby.py
```

Link to KR01 genome assembly, gene models `gff3` file, StringTie2 constructed transcript models, HISAT2 aligned RNA-seq reads, and Minimap2 aligned PacBio reads. 

```bash
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/Paulinella_micropora_KR01_nuclear.assembly.fasta
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/Paulinella_micropora_KR01_nuclear.assembly.fasta.fai
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/Paulinella_micropora_KR01_nuclear.stringtie2.estAbd.gtf
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/All_combined.HISAT2_RNAseq_mapping.coordSorted.bam
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/All_combined.HISAT2_RNAseq_mapping.coordSorted.bam.bai 
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_DNAseq_PacBio_RSII_minimap2/SRR10230249.coordsorted.bam
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_DNAseq_PacBio_RSII_minimap2/SRR10230249.coordsorted.bam.bai
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
```

## Subset data

**Genome scaffolds**

```bash
seqkit faidx --infile-list scaffolds_to_extract.txt \
  Paulinella_micropora_KR01_nuclear.assembly.fasta \
  > extracted_scaffolds.fasta
seqkit faidx extracted_scaffolds.fasta
```

**Gene models**

```bash
./grepf_column.py \
  -i Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3 \
  -f scaffolds_to_extract.txt \
  -o extracted_geneModels.gff3
```

**StringTie2 transcripts**

```bash
./grepf_column.py \
  -i Paulinella_micropora_KR01_nuclear.stringtie2.estAbd.gtf \
  -f scaffolds_to_extract.txt \
  -o extracted_stringtie2Transcripts.gtf
```

**Hisat2 aligned RNA-seq**

```bash
samtools view -b -L <(awk -F'\t' '{print $1"\t0\t"$2}' extracted_scaffolds.fasta.fai) \
  All_combined.HISAT2_RNAseq_mapping.coordSorted.bam \
  > extracted_Hisat2RNAseqReads.bam
samtools index extracted_Hisat2RNAseqReads.bam
```

## Download subset-ed data

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2022-01-13/all_sets_scaffolds_visualization"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="extracted_*" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

## Combine manually corrected gene models

Combine the gene models created for each set of manually corrected genes.

```bash
GFF="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.genome.gff3"
cat ../manually_correct_gene_models/$GFF \
    ../../2021-01-19/manually_correct_gene_models/$GFF \
    > allsets_combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.genome.gff3
```

## Visualize genes using IGV

Visualize each datasets along the scaffolds which have manually corrected gene models using IGV. 

IGV was run using KR01:

1. Reference genome
2. Gene models (`gff3` file; where the initial gene were extracted from)
3. StringTie2 transcripts (`gff3` file; with estimated abundances added to each feature)
4. `BAM` file of RNA-seq reads aligned against the genome using HISAT2 (same aligned reads used by StringTie2)
4. `BAM` file of PacBio reads aligned against the genome using Minimap2
4. Manually corrected gene models (`gff3` file)

IGV setting are saved in `igv_session.xml`

