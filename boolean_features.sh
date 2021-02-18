#!/bin/bash

# sorting the annotation files by gene
for chr in {1..22} X Y; do
	grep -v \# ../finalAnnot.chr${chr}.txt  | sort -k 4 > ../finalAnnot.SortByGene.chr${chr}.txt
done

exit
