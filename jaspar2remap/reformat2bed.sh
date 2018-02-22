#usage: bash reformat2bed.sh /path/to/where/all/jaspar/core/vertebrates/profiles/are/stored /path/to/where/all/jaspar/pwms/needed/are/stored

## ARGUMENT PARSING
profiles_dir=;
if [ -z "$1" ]
    then
        echo "No input folder containing the JASPAR CORE vertebrate profiles specified. Exiting...";
        exit 2
    else
        profiles_dir=$1
fi

pwm_dir=;
if [ -z "$2" ]
    then
        echo "No input folder containing the JASPAR TFBS predictions. Exiting...";
        exit 2
    else
        pwm_dir=$2
fi

#get current path
script_path=$PWD

#create out_dir unless exists
mkdir -p $pwm_dir

#download JASPAR TFBS predictions
cd $pwm_dir
for matrix_id in $(ls | cut -d "." -f 1,2 | sort | uniq); do
    bed_file=$matrix_id".bed"
    #echo "reformating $matrix_id..."
    echo "$script_path/JASPAR-UCSC-tracks/fetch_binding_sites.py -i $pwm_dir -p $profiles_dir -m $matrix_id -o $bed_file"
done