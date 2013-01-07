#!/bin/bash
shopt -s nullglob
for f in *.fa
do
	gene=${f%.*}
	../make_hmmdb_from_gene.sh $gene
done
