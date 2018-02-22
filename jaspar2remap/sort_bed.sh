# This script assumes that the BED files are zipped and will create a folder where the files are unzipped (as being the "raw" data to be used from now on) and another folder (out_dir) where the sorted files are output.

# This script is called by the sort_launcher.sh script

## ARGUMENT PARSING
in_dir=$1
data_dir=$2/${in_dir##*/}
out_dir=$3/${in_dir##*/}

echo "Processing folder: ${in_dir##*/}"

if [ ! -d "$data_dir" ]
    then
 
      #create folder structure 
       mkdir -p $data_dir
       mkdir -p $out_dir
          
       #unzip and sort
       for bed in $(find $in_dir -name "*.bed.gz")
           do
             gunzip -c $bed > $data_dir/${in_dir##*/}.bed                  
             sortBed -i $data_dir/${in_dir##*/}.bed > $out_dir/${in_dir##*/}.bed.sorted
       done

else
    echo "Folder $in_dir already exists in results. Skipping...";
fi

