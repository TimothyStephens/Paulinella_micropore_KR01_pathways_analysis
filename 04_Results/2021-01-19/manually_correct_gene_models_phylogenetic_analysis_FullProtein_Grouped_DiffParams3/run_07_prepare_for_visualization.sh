#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate r4_env; set -eu


export PATH="$PATH:/home/timothy/programs/AliStat_1.12"
export PATH="$PATH:/home/timothy/programs/TreeViewer-Linux-x64"

IDS="SuppTable_S2.transcript_info.txt"


##
## Run AliStat and crate annotation file for *FULL* alignments.
##
while read KO; do
	ALN="${KO}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln"
	TREE="${ALN}.contree"
	
	echo "## $ALN"
       	NSEQS=$(grep -c '>' "$ALN" || echo 0)
       	if [ $NSEQS -gt 3 ]; then
		#AliStat
		#alistat "${ALN}" 6  -o "${ALN}.AliStat" -t -i 1>"${ALN}.AliStat.log" 2>&1
		#R CMD BATCH "${ALN}.AliStat.Cumulative_Cc.R"
		#R CMD BATCH "${ALN}.AliStat.Histogram_Cc.R"
		#R CMD BATCH "${ALN}.AliStat.Histogram_Cr.R"
		
		#Annotation
		cat "${ALN}" | grep '>' \
		 | sed -e 's/>//' -e 's/^MSTRG/NA-NA-NA-NA-NA-MSTRG/' -e 's/^Paulinella_/NA-NA-NA-NA-NA-Paulinella_/' \
		 | awk '
		        BEGIN{
		                print "SeqName\tColor\tLevel_1\tLevel_2\tLevel_3\tLevel_4\tLevel_5\tLevel_6"
		        } {
		                split($1,a,"-");
		                     if(a[6]~"MSTRG" || a[6]~"Paulinella") {COLOR="#e41a1c"}
		                else if(a[1]=="Archaea")   {COLOR="#4daf4a"}
		                else if(a[1]=="Bacteria")  {COLOR="#377eb8"}
		                else if(a[1]=="Eukaryota") {COLOR="#f781bf"}
		                else if(a[1]=="Viruses")   {COLOR="#a65628"}
		                else                       {COLOR="#484848"};
		                print $1"\t"COLOR"\t"a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]
		        }' \
		| sed -e 's/^NA-NA-NA-NA-NA-MSTRG/MSTRG/' -e 's/^NA-NA-NA-NA-NA-Paulinella_/Paulinella_/' \
		| ./add_value_to_table.py -a <(awk -F',' 'BEGIN{print "SeqName\tAlistat_Cr"} NR>1 {print $2"\t"$4}' "${ALN}.AliStat.Table_1.csv") \
		| ./add_value_to_table.py -a <(awk -F'\t' 'BEGIN{print "SeqName\tcrTP\tmtTP"} $1!="seqid" {print $0}' ../organelle_targeting_peptides/*.organelle_TP.txt) -d $'No\tNo' \
		> "${ALN}.annots.txt"
       	else
       	       	echo "   - Too few seqs! (only $NSEQS seqs)"
       	fi
done < <(awk -F'\t' '$1~"MSTRG" {print $2}' "$IDS" | sort | uniq)



##
## Run AliStat and crate annotation file for *TRIMMED* alignments.
##
while read KO; do
        ALN="${KO}.blastp_ALL.outfmt6.parsed.filtered.uniq.faa.mafft.aln.trimal1"
        TREE="${ALN}.contree"

        echo "## $ALN"
        NSEQS=$(grep -c '>' "$ALN" || echo 0)
        if [ $NSEQS -gt 3 ]; then
                #AliStat
                #alistat "${ALN}" 6  -o "${ALN}.AliStat" -t -i 1>"${ALN}.AliStat.log" 2>&1
                #R CMD BATCH "${ALN}.AliStat.Cumulative_Cc.R"
                #R CMD BATCH "${ALN}.AliStat.Histogram_Cc.R"
                #R CMD BATCH "${ALN}.AliStat.Histogram_Cr.R"

                #Annotation
                cat "${ALN}" | grep '>' | sed -e 's/>//' -e 's/^MSTRG/NA-NA-NA-NA-NA-MSTRG/' -e 's/^Paulinella_/NA-NA-NA-NA-NA-Paulinella_/' \
                 | awk '
                        BEGIN{
                                print "SeqName\tColor\tLevel_1\tLevel_2\tLevel_3\tLevel_4\tLevel_5\tLevel_6"
                        } {
                                split($1,a,"-");
                                     if(a[6]~"MSTRG" || a[6]~"Paulinella") {COLOR="#e41a1c"}
                                else if(a[1]=="Archaea")   {COLOR="#4daf4a"}
                                else if(a[1]=="Bacteria")  {COLOR="#377eb8"}
                                else if(a[1]=="Eukaryota") {COLOR="#f781bf"}
                                else if(a[1]=="Viruses")   {COLOR="#a65628"}
                                else                       {COLOR="#484848"};
                                print $1"\t"COLOR"\t"a[1]"\t"a[2]"\t"a[3]"\t"a[4]"\t"a[5]"\t"a[6]
                        }' \
		| sed -e 's/^NA-NA-NA-NA-NA-MSTRG/MSTRG/' -e 's/^NA-NA-NA-NA-NA-Paulinella_/Paulinella_/' \
		| ./add_value_to_table.py -a <(awk -F',' 'BEGIN{print "SeqName\tAlistat_Cr"} NR>1 {print $2"\t"$4}' "${ALN}.AliStat.Table_1.csv") \
		| ./add_value_to_table.py -a <(awk -F'\t' 'BEGIN{print "SeqName\tcrTP\tmtTP"} $1!="seqid" {print $0}' ../organelle_targeting_peptides/*.organelle_TP.txt) -d $'No\tNo' \
                > "${ALN}.annots.txt"
        else
                echo "   - Too few seqs! (only $NSEQS seqs)"
        fi
done < <(awk -F'\t' '$1~"MSTRG" {print $2}' "$IDS" | sort | uniq)



echo ""; echo "Done preparing for tree visualization!"

