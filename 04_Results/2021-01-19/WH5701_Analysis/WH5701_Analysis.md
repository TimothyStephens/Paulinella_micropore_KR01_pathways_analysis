# Analysis of KEGG Pathways in ***Synechococcus*** sp. WH5701

Check which proteins from each of the pathways are present in *Synechococcus* sp. WH5701. This will give us an idea of which proteins are important to the chromatophore (and conversely which aren't). Use KAAS to annotate proteins; use GHOSTX and GHOSTZ run using the "for Eukartotes" and "for Prokaryotes" datasets.

## Setup analysis directory

Link WH5701 proteins which I had already downloaded previously.

```bash
ln -s ../../../../../genome_data/Alpha-Cyanobacteria/databases/Synechococcus_sp_wh_5701.ASM15304v1.pep.all.fa
```

Link to custom scripts for text parsing.

```bash
ln -s ../../../02_Scripts/groupby.py
ln -s ../../../02_Scripts/add_value_to_table.py
```

Setup bash environment.

```bash
conda activate py27
```

## Analysis 1

```bash
cat /scratch/timothy/databases/KEGG/ko00001.keg \
  | awk 'BEGIN{P=0} { if($1=="C"){ if($2=="03030"){P=1} else{P=0}} else{if(P==1){LINE=$3; for (i=4; i<NF; i++){LINE=LINE" "$i}; print $2"\t"LINE}} }' \
  | ./add_value_to_table.py -a <(./groupby.py --delim_groups ";" -i <(awk '$2!=""{print $2"\t"$1}' Synechococcus_sp_wh_5701.ASM15304v1.pep.all.KAAS_Combined.ko | sort)) \
  | awk -F'\t' '$3!=""' \
  | sort -k1,1
```

```bash
cat /scratch/timothy/databases/KEGG/ko00001.keg \
  | awk 'BEGIN{P=0} { if($1=="C"){ if($2=="00230"){P=1} else{P=0}} else{if(P==1){LINE=$3; for (i=4; i<NF; i++){LINE=LINE" "$i}; print $2"\t"LINE}} }' \
  | ./add_value_to_table.py -a <(./groupby.py --delim_groups ";" -i <(awk '$2!=""{print $2"\t"$1}' Synechococcus_sp_wh_5701.ASM15304v1.pep.all.KAAS_Combined.ko | sort)) \
  | awk -F'\t' '$3!=""' \
  | sort -k1,1
```

```bash
cat /scratch/timothy/databases/KEGG/ko00001.keg \
  | awk 'BEGIN{P=0} { if($1=="C"){ if($2=="00240"){P=1} else{P=0}} else{if(P==1){LINE=$3; for (i=4; i<NF; i++){LINE=LINE" "$i}; print $2"\t"LINE}} }' \
  | ./add_value_to_table.py -a <(./groupby.py --delim_groups ";" -i <(awk '$2!=""{print $2"\t"$1}' Synechococcus_sp_wh_5701.ASM15304v1.pep.all.KAAS_Combined.ko | sort)) \
  | awk -F'\t' '$3!=""' \
  | sort -k1,1
```

