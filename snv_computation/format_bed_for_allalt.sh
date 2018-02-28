#!/usr/bin/env bash

################################################################################
#
################################################################################

bed=$1;
output=$2;
genome=$3;
sizes=$4;
src=$5;

tmp=`mktemp`;
bedtools slop -b 60 -i $bed -g $sizes > $tmp;
fasta=`mktemp`;
bedtools getfasta -fi $genome -bed $tmp -fo $fasta -name;
rm $tmp;
python $src/format_for_allalt.py -b $bed -f $fasta > $output;
rm $fasta;
