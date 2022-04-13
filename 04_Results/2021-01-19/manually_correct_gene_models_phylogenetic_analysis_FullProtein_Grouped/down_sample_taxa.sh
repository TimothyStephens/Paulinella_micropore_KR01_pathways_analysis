#!/bin/bash
i=$1
cat $i | awk '{print $2}' | awk -F '-' '{print $1}' | awk '!seen[$0]++' > $i.lvl1
cat $i | awk '{print $2}' | awk -F '-' '{print $1"-"$2}' | awk '!seen[$0]++' > $i.lvl2
cat $i | awk '{print $2}' | awk -F '_' '{print $1}' | awk '!seen[$0]++' > $i.lvl3
numClass=`wc -l $i.lvl2 | awk '{print $1}'`
numGen=$((180/numClass))
for z in `cat $i.lvl2`; do grep -- $z $i.lvl3 | grep -v -- "----" | grep -v -- "---" |head -${numGen}; done > $i.genera
num=`wc -l $i.genera | awk '{print $1}'`
if [ $num -lt 25 ]
  then
  for z in `cat $i.genera`; do grep "$z" $i | head -2 | awk '{print $2}'; done > $i.parsed
  else
  for z in `cat $i.genera`; do grep "$z" $i | head -1 | awk '{print $2}'; done > $i.parsed
fi
cat $i.parsed | awk '!seen[$0]++' > $i.tmp
mv $i.tmp $i.parsed
rm $i.lvl* $i.genera
