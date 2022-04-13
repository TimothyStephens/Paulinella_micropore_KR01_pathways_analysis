# Manually check and correct KR01 gene models of interest if needed

Use IGV to manually rebuild the KR01 genes using transcripts constructed by StringTie2. This will correct any mis-predicted regions as the StringTie2 transcripts are based on RNA-seq evidence. ORFs will be re-predicted in the new transcripts even if a gene exactly matches the StringTie2 transcript; this will act as a double check for the ORFs structure + will add UTRs to each gene. These new gene models will be used in all downstream analysis.

## Setup analysis directory

Link files with KR01 gene IDs and KO numbers.

```bash
ln -s ../../../01_Data/2022-01-13/Histidine_subset.txt
```

Link to KR01 genome assembly, gene models `gff3` file, and StringTie2 constructed transcript models. 

```bash
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/Paulinella_micropora_KR01_nuclear.assembly.fasta
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/Paulinella_micropora_KR01_nuclear.assembly.fasta.fai
ln -s ../../../../../genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/Paulinella_micropora_KR01_nuclear.stringtie2.estAbd.gtf
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/All_combined.HISAT2_RNAseq_mapping.coordSorted.bam
ln -s ../../../../../Paulinella_micropora_KR01/Read_Mapping/Genome_RNAseq_HISAT2_StringTie2_NoRefGff/All_combined.HISAT2_RNAseq_mapping.coordSorted.bam.bai 
```

Link to custom scripts for parsing text files.

```bash
ln -s ../../../02_Scripts/grepf_column.py
ln -s ../../../02_Scripts/add_value_to_table.py
ln -s ../../../02_Scripts/groupby.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/bedtools-2.29.2/bin"
export PATH="$PATH:/home/timothy/programs/seqkit_v0.15.0"
export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
```

## Extract just features associated with our target genes

- Extract gene model and StringTie2 `gff3` features that overlap with out target genes.

- Will visualize this information (but using the original files) using IGV so I can sort through it manually and correct the gene models where necessary.  
- The StringTie2 features that we extract here will be used as a basis/scaffold for manually correcting the gene models. It is easier to mix/update StringTie2 records then build a whole set of gene model `gff` features by hand. 
- Also, even if the original gene model is correct we are still interested in the full transcript, including UTRs.

Get bed formatted mRNA features from KR01 predicted genes.

```bash
PRED_GENES="Paulinella_micropora_KR01_nuclear.augustus_190918v4.gff3"
awk -F'\t' '$3=="mRNA" {print $1"\t"$4-1"\t"$5"\t"$9}' "$PRED_GENES" \
        | sed -e 's/ID=Paulinella_micropora_KR01_nuclear___//' -e 's/;Parent=.*//' \
        > "$PRED_GENES.mRNA_features.bed"
```

Get mRNA features for the KR01 genes that are annotated to each pathway. 

Overlap these features with transcripts constructed with StringTie2.

```bash
STRINGTIE_GTF="Paulinella_micropora_KR01_nuclear.stringtie2.estAbd.gtf"

## Get just "transcript" features from gtf file -> bed format
bedtools sort -i "$STRINGTIE_GTF" | awk -F'\t' '$3=="transcript"' > "$STRINGTIE_GTF.transcript_features.bed"

function get_overlap() {
        INFO="$1"
        ./grepf_column.py -c 4 \
            -i <(bedtools sort -i "$PRED_GENES.mRNA_features.bed") \
            -f <(awk -F'\t' '$2~"^g" {print $2}' "$INFO" | sort | uniq) \
            -o "$INFO.mRNA_features.bed"
        bedtools intersect -wa -wb -a "$INFO.mRNA_features.bed" \
            -b "$STRINGTIE_GTF.transcript_features.bed" \
          | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$9"\t"$13}' \
          | sed -e 's/gene_id.*transcript_id "//' -e 's/";.*//' \
          > "$INFO.overlapping_transcript_features.bed"

        # Check every KR01 gene have overlapping transcript features.
        cat "$INFO.overlapping_transcript_features.bed" | cut -f4 | sort | uniq | wc -l
        cat $INFO.mRNA_features.bed | cut -f4 | sort | uniq | wc -l

        # Group overlapping transcripts and add them to the "info" list (gets it all 
        # in the correct order so we can join [cut/paste] it to our exisiting excel table)
        ./add_value_to_table.py -d NA \
            -i <(cut -f2 "$INFO") \
            -a <(cut -f4,7 "$INFO.overlapping_transcript_features.bed" \
                   | ./groupby.py --delim_groups ',') \
            -o "$INFO.overlapping_transcript_features.ordered.txt"
}

get_overlap "Histidine_subset.txt"
```

Extract StringTie2 `gtf` features that overlap with out target genes. We edit the gene features in these files manually, as required, based off what I find after visulaization of each gene using IGV.

```bash
INFO="Histidine_subset.txt"
grep -F -f <(awk -F'\t' '{print "\""$7"\""}' "$INFO.overlapping_transcript_features.bed") \
    "$STRINGTIE_GTF" > "$INFO.overlapping_transcript_features.gtf"
```

## Subset data for visualization

Get just the features/reads from each file that are associated with the target scaffolds (i.e., the scaffolds with our target genes). Getting just the features on the scaffolds of interest means that we don't have to keep handling the full datasets which (in the case of the aligned reads) can be quite large.

```bash
ln -s ../initial_sequences_info/Histidine_subset.txt.extra_info_added.txt
```



```bash
awk -F'\t' '$2~"^g" {print "Paulinella_micropora_KR01_nuclear___"$6}' \
  Histidine_subset.txt.extra_info_added.txt \
  > scaffolds_to_extract.txt 
```

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
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2022-01-13/manually_correct_gene_models"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="extracted_*" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

## Visualize genes using IGV

Visualize each "initial" target gene using IGV. Pick the StringTie2 gene models that best represents the initial gene. We will repredict the ORFs for each just to double check. Should not affect most genes but will add info about UTRs to our predictions which might be helpful down the line.

IGV was run using KR01:

1. Reference genome
2. Gene models (`gff3` file; where the initial gene were extracted from)
3. StringTie2 transcripts (`gff3` file; with estimated abundances added to each feature)
4. `BAM` file of RNA-seq reads aligned against the genome using HISAT2 (same aligned reads used by StringTie2)

IGV setting are saved in `igv_session.xml`

**Assumptions/Processes**

For each initial protein/gene:

- A StringTie2 transcript was chosen that most closely matched the intron/exon structure of the initial gene (often the two exactly matched; the StringTie2 transcript would obviously have extra regions at the 5- and 3-prime UTRs).
- If the initial gene had obvious mis-predicted regions or slightly altered intron/exon boundaries (based on manually checking the gene model against the aligned reads) then the StringTie2 transcripts with the highest estimated abundance was chosen as the "correct" gene model for the initial gene. Assuming that the StringTie2 transcript also appeared to be free from mis-predictions and was congruent with the aligned reads. 
- If the best StringTie2 transcript model had terminal exons that had much lower support then the rest of the gene then they were removed. This often occurs at the 5-prime end when a gene has an SL sequence attached; one base of the SL with align to the genome and form a new exon that is either 1bp in length of only a few bp with very strange read coverage patterns. 

New manually corrected gene models are in `*.overlapping_transcript_features.manually_fixed.gtf` files.

## Extract ORFs from manually checked transcripts

After manually checking transcripts using IGV (fixing problems when they were identified) the next step is to predict the most likely ORF in each transcript. 

- One ORF per transcript.

- Only predict ORFs on the same strand as the transcript (ORFs on opposite strands cant usually be propagated [unless from a single exon transcript] to the genome so are normally not helpful).

- Have to use a custom script to force 5_prime_partial genes to next available start codons (assume most of these genes are full length and that start codon prediction is a little inaccurate).
- Use BLASTP and HMMSCAN searches against SwissProt and PfamA (respectivaly) to guide ORF prediction.

```bash
cat *.overlapping_transcript_features.manually_fixed.gtf | sort -k1,1 -k4,4n -k3,3r | uniq > combined.manually_fixed.gtf

./run_Transdecoder.sh
```

## Download manually corrected genes

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2022-01-13/manually_correct_gene_models"
rsync -zarv --prune-empty-dirs --relative \
 --include="*/" --include="*.sh*" --include="*.pl" --include="*.overlapping_transcript_features.manually_fixed.gtf" --include="combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.*" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

Can redo IGV visualization using the `combined.manually_fixed.transcripts.fasta.transdecoder.forceStart.genome.gff3` file. This way we can double check that we constructed our new new gene models correctly.
