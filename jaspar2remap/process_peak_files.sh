#This script is called by the process_launcher.sh script

## ARGUMENT PARSING
in_dir=$1
out_dir=$2

#create folder structure
if [ ! -d "$out_dir" ]
   then
       mkdir -p $out_dir
fi

# get the peak file from the current folder
find $in_dir -type f -name "*.narrowPeak" -exec cat {} > "$out_dir/${in_dir##*/}.bed" \;  
# sort it
sort -k1,1 -k2,2n "$out_dir/${in_dir##*/}.bed" > "$out_dir/${in_dir##*/}.bed.sorted"
#and merge the peaks
mergeBed  -i "$out_dir/${in_dir##*/}.bed.sorted" > "$out_dir/${in_dir##*/}.bed.sorted.merged"

