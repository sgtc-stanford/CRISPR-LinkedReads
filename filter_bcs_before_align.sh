#!/bin/bash

PERC=$1

for file in *R1*.gz; do zcat $file | sed -n '2~4p' | cut -c1-16 >> bcs.txt; done

cat bcs.txt | sort | uniq -c | sort -nr | sed "s/^[ \t]*//" | sed "s/\ /\t/" > bcs.counts.txt

cat bcs.counts.txt | grep -v 'N' > bcs.counts.noN.txt

filter_bc_by_perc.py bcs.counts.noN.txt $PERC


