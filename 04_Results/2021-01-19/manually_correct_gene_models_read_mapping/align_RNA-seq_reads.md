# Align RNA-seq reads against KR01 CDS with the manually corrected sequences added

Align RNA-seq reads against all KR01 CDS + manually corrected transcripts (not CDS); the "initial" CDS were removed from the set before mapping reads. Align against all sequences so we don't inadvertently force the aligner to suboptimally map reads to the wrong transcripts; by giving the aligner all CDS we are sure the reads that align to the transcripts are in the optimal position and are at the best positions possible.

## Setup analysis directory

Link to files listing KR01 target genes.

```bash
ln -s ../../../01_Data/2021-01-19/DNA_Replication.txt
ln -s ../../../01_Data/2021-01-19/Purine.txt
ln -s ../../../01_Data/2021-01-19/Pyrimidine.txt
```

Link manually corrected transcripts (full mRNA not just CDS) + original KR01 CDS.

```bash
ln -s ../manually_correct_gene_models/combined.manually_fixed.transcripts.fasta
ln -s /scratch/timothy/genome_data/Paulinella_micropora_KR01/nuclear_genome/databases/Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna
```

Link to samples file which lists the pairs of read files for each KR01 sample that we want to align against our database.

```bash
ln -s /scratch/timothy/Paulinella_micropora_KR01/Data_from_PRJNA568118/samples/samples.txt
```

Link custom scripts for parsing text files + search aligned reads for SL sequences.

```bash
ln -s ../../../02_Scripts/grepf_fasta.py
ln -s ../../../02_Scripts/mapped_reads_with_SL_sequences.py
```

Setup bash environment.

```bash
conda activate py27

export PATH="$PATH:/home/timothy/programs/samtools-1.11/bin"
```

## Setup CDS/transcript database for mapping

Remove KEGG genes from KR01 CDS and add in the new manually corrected transcripts.

```bash
./grepf_fasta.py -v \
    -f <(cat DNA_Replication.txt Purine.txt Pyrimidine.txt | awk -F'\t' '$2~"^g" {print "Paulinella_micropora_KR01_nuclear___"$2}' | sort | uniq) \
    -i Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.fna \
    -o Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.non_target.fna

cat combined.manually_fixed.transcripts.fasta \
    Paulinella_micropora_KR01_nuclear.augustus_190918v4.cds.non_target.fna \
  > Paulinella_micropora_KR01_nuclear.corrected_genes_v1.txt
```

Build Bowtie2 index/database.

```bash
./run_bowtie2-build.sh
```

Get `BED` features that covered the whole sequence (will be used later to extract reads that map to our manually corrected target transcripts).

```bash
samtools faidx combined.manually_fixed.transcripts.fasta
awk -F'\t' '{print $1"\t0\t"$2}' combined.manually_fixed.transcripts.fasta.fai > combined.manually_fixed.transcripts.fasta.whole_gene.bed
```

## Align reads using Bowtie2

Align paired-reads from each sample listed in `samples.txt` against the new sequence database we just created. Use both local and global read mapping.

```bash
./run_bowtie2.sh
./run_bowtie2-global.sh
```

Index each coord-sorted aligned `BAM` file and them merge them together into a big combined `BAM` file. Do this for both local and global results. 

```bash
./run_samtools_merge.sh
./run_samtools_merge-global.sh
```

## Extract reads aligned to target transcripts

```bash
PREFIX="combined.manually_fixed.transcripts.fasta"
```

Extract just reads that map to corrected gene models (local mapping)

```bash
samtools view -@ 6 -b \
    -L $PREFIX.whole_gene.bed \
    -o $PREFIX.All_combined.local_mapping.coordsorted.bam \
    All_combined.local_mapping.coordsorted.bam
samtools index $PREFIX.All_combined.local_mapping.coordsorted.bam
```

Extract just reads that map to corrected gene models (global mapping)

```bash
samtools view -@ 6 -b \
    -L $PREFIX.whole_gene.bed \
    -o $PREFIX.All_combined.global_mapping.coordsorted.bam \
    All_combined.global_mapping.coordsorted.bam
samtools index $PREFIX.All_combined.global_mapping.coordsorted.bam
```

## Find SL addition sites along target transcripts

```bash
PREFIX="combined.manually_fixed.transcripts.fasta"
```

Use a custom script to search through the aligned reads for SL addition sites (positions where the read data support the addition of the *Paulinella* spliced leader [SL] sequence into the transcript).

```bash
./mapped_reads_with_SL_sequences.py \
    -s CCGGCTTTTCTG \
    -i <(cut -f1 "$PREFIX.whole_gene.bed") \
    -b "$PREFIX.All_combined.local_mapping.coordsorted.bam" \
    -o "$PREFIX.SL_from_mapped_reads.txt" \
    --info "$PREFIX.SL_from_mapped_reads.SL_read_info.txt.gz" \
  2> "$PREFIX.SL_from_mapped_reads.log"
```

## Download results from server

```bash
WD="timothy@coral.rutgers.edu:/scratch/timothy/projects/0011_Paulinella_micropora_KR01_KEGG_pathway_analysis/03_Analysis/2021-01-19/manually_correct_gene_models_read_mapping"
rsync -zarvL --prune-empty-dirs --relative \
 --include="*/" --include="*.sh*" --include="*.stats" --include="samples.txt" --include="*.SL_from_mapped_reads.*" --exclude="*" \
 ${WD}/./ \
 . --dry-run
```

