#usage: bash process_launcher.sh /path/to/the/split/by/tf/peaks/folder process_peak_files.sh /path/to/the/result/folder 

## ARGUMENT PARSING
in_dir=;
if [ -z "$1" ]
    then 
	echo "No input file folder specified. Exiting..."; 
	exit 1
    else	
	in_dir=$1
fi

script_file=;
if [ -z "$2" ]
    then 
	echo "No script file specified. Exiting..."; 
	exit 1
    else	
	script_file=$2
fi

out_dir=;
if [ -z "$3" ]
    then 
	echo "No output folder specified. Exiting..."; 
	exit 1
    else	
	out_dir=$3
fi

## call the processing in parallel
cd $in_dir
find $in_dir -maxdepth 1 -mindepth 1 -type d | xargs -I {} --max-proc=50 bash $script_file {} $out_dir 
