#usage: bash sort_bed_simple.sh /path/to/the/aggregated/data /path/to/the/final/result/folder

##ARGUMENT PARSING
in_dir=;
if [ -z "$1" ]
    then
        echo "No input folder specified. Exiting...";
        exit 1
    else
        in_dir=$1
fi

out_dir=;
if [ -z "$2" ]
    then
        echo "No output folder specified. Exiting...";
        exit 1
    else
        out_dir=$2
fi

mkdir -p $out_dir

cd $in_dir
for bed in $(find ./ -name "*.bed")
   do 
     echo "Processing $bed"
     sort -k1,1 -k2,2n "$bed" > "$out_dir/${bed##*/}"
done

