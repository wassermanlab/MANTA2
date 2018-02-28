#!/usr/bin/env bash

################################################################################
# Script computing the potential impact of all SNVs that could occur in the
# TFBSs given in the input BED file. It uses a JASPAR PFM to compute impact
# scores. These scores are used to be stored in our MANTA database.
################################################################################

#src=PATH_TO_THE_REPOSITORY_CONTAINING_THIS_SCRIPT
src=/storage/scratch/anthoma/MANTA2/MANTA2/bin/MANTA_TFBS_snv_score_computation/
raven=$src/raven_scoring_snv_noalt_htt_tfbsrange_search_serial_allalt.pl

usage="""
\n
$0 -b <TFBS BED file> -p <JASPAR PFM> -g <genome fasta file> -s <genome size file> -o <output dir>\n
\n
    TFBS BED file     : BED file containing the positions of the TFBSs to
                        consider\n
    JASPAR PFM        : file containing the PFM in JASPAR PFM format\n
    genome fasta file : FASTA file containing the genome to consider\n
    genome size file  : file containing the size of each chromosome (as used by
                        BEDTOOLS)\n
    output dir        : name of the output directory that will contain the
                        results\n

\n
""";

TEMP=`getopt -o hb:p:o:g:s: -- "$@"`;
if [ $? != 0 ]
then
  echo "Terminating..." >&2;
  exit 1;
fi
eval set -- "$TEMP";

bed="";
pfm="";
outdir="";
fasta="";
sizes="";
while true
do
  case "$1" in
    -b) bed=$2; shift 2;;
    -p) pfm=$2; shift 2;;
    -o) outdir=$2; shift 2;;
    -g) fasta=$2; shift 2;;
    -s) sizes=$2; shift 2;;
    -h) echo -e $usage; exit;;
    --) shift; break;;
    *) echo "Internal error!"; exit 1;;
  esac
done

# Checking options and parameters
if [ ! -e $bed ]
then
    echo "File $bed does not exist!";
    exit 1;
fi;
if [ ! -e $pfm ]
then
    echo "PFM file $pfm does not exist!";
    exit 1;
fi;
if [ ! -e $fasta ]
then
    echo "Genome fasta file $fasta does not exist!";
    exit 1;
fi;
if [ ! -e $sizes ]
then
    echo "Genome size file $sizes does not exist!";
    exit 1;
fi;
if [ ! -e $outdir ]
then
    mkdir -p $outdir;
fi;


name=`basename $bed .bed`;
tmp=`mktemp -p $outdir`;
tmpbed=`mktemp -p $outdir`;
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_"NR,$5,$6}' $bed > $tmpbed;
bash $src/format_bed_for_allalt.sh $tmpbed $tmp $fasta $sizes $src;
rm $tmpbed;
mkdir -p $outdir/$mat_name;
out=$outdir/$mat_name/$name"_SNV_impacts.txt";
perl $raven -i $tmp -m $pfm -o $outdir -of $out -c $src;
