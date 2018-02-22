#usage: bash get_JASPAR_CORE_vertebrate_profiles.sh /path/to/where/all/jaspar/core/vertebrates/profiles/are/stored

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
curl -O http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.zip
unzip JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.zip