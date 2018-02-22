#usage: bash get_JASPAR_TFBS_predictions.sh /path/to/where/all/jaspar/pwms/needed/are/stored

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

#download JASPAR TFBS predictions
cd $out_dir
for file in $(curl -s http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg38/tsv/ |
              grep href |
              sed 's/.*href="//' |
              sed 's/".*//' |
              grep '^[a-zA-Z].*'); do
    echo "downloading $file..."
    curl -s -O http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2018/hg38/tsv/$file
done