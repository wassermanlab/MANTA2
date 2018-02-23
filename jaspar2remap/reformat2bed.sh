#usage: bash reformat2bed.sh /path/to/where/all/jaspar/pwms/needed/are/stored

## ARGUMENT PARSING
out_dir=;
if [ -z "$1" ]
    then
        echo "No output folder for the bed files specified. Exiting...";
        exit 1
    else
        out_dir=$1
fi

#create out_dir unless exists
mkdir -p $out_dir

#reformat JASPAR TFBS predictions from TSV to BED
cd $out_dir
for matrix_id in $(ls | cut -d "." -f 1,2 | sort | uniq); do
    tsv_file=$matrix_id".tsv.gz"
    bed_file=$matrix_id".bed"
    echo "reformating $matrix_id..."
    gzcat $tsv_file | tail -n +2 | perl -e 'while(<>){chomp;@a=split("\t",$_);print $a[0]."\t".$a[1]-1."\t".$a[2]."\t".$a[3]."\t".$a[4]."\t".$a[-1]."\n";}' > $bed_file
done