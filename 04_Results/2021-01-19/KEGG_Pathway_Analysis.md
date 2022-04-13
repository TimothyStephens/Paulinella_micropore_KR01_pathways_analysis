# Analysis of *Paulinella micropore* KR01 KEGG pathways

`initial_sequences_blastp_UniProt_KEGG_seqs/`

- For each gene, compare using BLASTP against other sequences from UniProt annotated with the same KO numbers. 
- Collect metrics about query and subject hit coverage in order to assess how our proteins compare against other sequences of the same type.
  - If protein has low subject coverage then it might be false positive
  - If protein is much longer than subject (but had high coverage) then we could look more into it to see what the function of the extra sequence is. 

`initial_sequences_info/`

- For each gene in target KEGG pathways collect extra info and add it to a table.

`initial_sequences_plot_read_mapping/`

- Compare initial proteins against UniProt KEGG Orthologs. 
- Limit analysis to protein<->KO_orthologs that we expect based on the annotations we already have.
- Extract coverage information about the top hits to each KR01 protein. 
    - Query/Subject coverage can be used to infer false positives or if a protein has extended 5- or 3-prime ends. might indicate that it has a known or unknown crTP.

---

`manually_correct_gene_models/`

- Annotate each gene with BLASTP hits to UniProt orthologs and visualize using IGV. 
    - Will help identify any problems with the gene models (which were automatically predicted)
    - Fix identified problems manually and recreate gene models/proteins/CDS for each corrected gene (used for subsequent analysis).

`manually_correct_gene_models_blastp_UniProt_KEGG_seqs/`

- For each manually corrected gene, compare using BLASTP against other sequences from UniProt annotated with the same KO numbers. 
- Collect metrics about query and subject hit coverage in order to assess how our proteins compare against other sequences of the same type.
  - If protein has low subject coverage then it might be false positive
  - If protein is much longer than subject (but had high coverage) then we could look more into it to see what the function of the extra sequence is. 

`manually_correct_gene_models_crTP/`

- Search the "new" manually corrected proteins for crTP using the crTP alignment/HMM that we have available for KR01.

`manually_correct_gene_models_phylogenetic_analysis_FullProtein/`

- Perform phylogenetic analysis with the new manually corrected genes.
  - BLAST new proteins against a taxonomically diverse database
  - Down sample hits using taxonomic information
  - Align down sampled proteins
  - Build tree using down sampled proteins + new KR01 protein

`manually_correct_gene_models_phylogenetic_analysis_FullProtein_Grouped/`

- Perform phylogenetic analysis with the new manually corrected genes grouped by KO annotation (i.e., build tree with all proteins annotated with the same KO number)
  - Use BLAST results from `manually_correct_gene_models_phylogenetic_analysis_FullProtein/`
  - Down sample hits using taxonomic information
  - Combined BLAST results into groups based off of KO numbers. 
  - Remove redundant/duplicate sequences in combined sets. 
  - Align combined down sampled proteins
  - Build tree using down sampled proteins + new KR01 protein

`manually_correct_gene_models_phylogenetic_analysis_LongTermini/`

- Perform phylogenetic analysis with long termini of the new manually corrected genes.
  - Extract termini regions of new proteins that are longer then we would expect (based on BLAST comparison with UniProt KEGG proteins).
  - BLAST new proteins against a taxonomically diverse database
  - Down sample hits using taxonomic information
  - Align down sampled proteins
  - Build tree using down sampled proteins + new KR01 protein

`manually_correct_gene_models_read_mapping/`

- Align RNA-seq reads to KR01 CDS + manually corrected mRNA (CDS+UTRs); initial KEGG CDS were removed before the manually corrected sequences were added. 

`manually_correct_gene_models_transcript_visualization/`

- Aggregate files for manually corrected proteins for visualization. 
- Generate annotation files with coords of Pfam domains and BLAST hits to UniProt KEGG orthologs. 

