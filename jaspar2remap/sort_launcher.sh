#usage: bash sort_launcher.sh /path/to/original/bed/files/ /path/to/the/sort_bed.sh/script /path/to/the/folder/to/store/unzipped/files/ /path/to/the/folder/to/store/sorted/files

## ARGUMENT PARSING
in_dir=;
if [ -z "$1" ]
    then 
	echo "No input folder specified. Exiting..."; 
	exit 1
    else	
	in_dir=$1
fi

script=;
if [ -z "$2" ]
    then 
	echo "No script file specified. Exiting..."; 
	exit 1
    else	
	script=$2
fi

data_dir=;
if [ -z "$3" ]
    then 
	echo "No data folder where to store the unzipped files specified. Exiting.."; 
	exit 1
    else	
	data_dir=$3
fi

out_dir=;
if [ -z "$4" ]
    then
        echo "No output folder for the sorted files specified. Exiting...";
        exit 1
    else
        out_dir=$4
fi

#launch 50 processes in parallel
cd $in_dir
find $in_dir -maxdepth 1 -mindepth 1 -type d | xargs -I {} --max-proc=50 bash $script {} 
