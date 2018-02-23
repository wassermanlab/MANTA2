#usage: bash get_tfbs_predictions.sh /path/to/the/original/BED/files

## ARGUMENT PARSING
out_dir=;
if [ -z "$1" ]
    then
        echo "No output folder for the downloaded files specified. Exiting...";
        exit 1
    else
        out_dir=$1
fi

#create out_dir unless exists
mkdir -p $out_dir

#download TFBS predictions
url="http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg38/tsv/"
cd $out_dir

for file in $(curl -s $url       | 
              grep href          |
              sed 's/.*href="//' |
              sed 's/".*//'      |
              grep '^[a-zA-Z].*'); do
    echo "downloading $file..."
    curl -s -O $url$file
    bed_file=$(ls $file | cut -c 1-8)".bed"
    gzcat $file | tail -n +2 | awk '{print $1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $7}' > $bed_file
    gzip $bed_file
    rm $file
done