# usage:
#    for the ENCODE datasets: bash copy_folders_based_on_file.sh /path/to/the/ENCODE/peaks/folder /path/to/the/encode_datasets_for_manta.txt ENCODE /path/to/the/target/folder
#    for the GEO&AE datasets: bash copy_folders_based_on_file.sh /path/to/the/GEO&AE/peaks/folder /path/to/the/geo_datasets_for_manta.txt GEO /path/to/the/target/folder


## ARGUMENT PARSING
in_dir=;
if [ -z "$1" ]
    then
        echo "No input folder specified. Exiting...";
        exit 1
    else
        in_dir=$1
fi


datasets=;
if [ -z "$2" ]
    then
        echo "No file containing a list of datasets to copy was specified. Exiting...";
        exit 1
    else
        datasets=$2
fi

data_type=;
if [ -z "$3" ]
    then
        echo "No data type specified. Either ENCODE or GEO allowed. Exiting...";
        exit 1
    else
        data_type=$3
fi

target=;
if [ -z "$4" ]
    then
        echo "No target folder specified. Exiting...";
        exit 1
    else
        target=$4
fi

## DEFAULT DATA TYPE
dt="ENCODE"

mkdir -p $target

cd $in_dir
while IFS=$'\t' read -r -a line
    do
      if [ "$data_type" == "$dt" ]
         then 
           #encode 
           dataset="${line[10]}"
           folder=$(find . -mindepth 1 -maxdepth 1 -type d -iname "$dataset*")
      else  
           #geo
           dataset="${line[11]}"
           folder=$(find . -mindepth 1 -maxdepth 1 -type d -iname "$dataset")
     fi
       #  echo "$folder"  
     cp -R $folder $target
done < "$datasets"

