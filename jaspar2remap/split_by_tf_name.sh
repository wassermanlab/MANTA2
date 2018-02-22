#usage: bash split_by_tf_name.sh /path/to/folder/where/peaks/where/copied ENCODE /path/to/the/output/folder

## ARGUMENT PARSING
in_dir=;
if [ -z "$1" ]
    then
        echo "No input folder specified. Exiting...";
        exit 1
    else
        in_dir=$1
fi

data_type=;
if [ -z "$2" ]
    then
        echo "No data type specified. Either ENCODE or GEO. Exiting...";
        exit 1
    else
        data_type=$2
fi


out_dir=;
if [ -z "$3" ]
    then
        echo "No output folder specified. Exiting...";
        exit 1
    else
        out_dir=$3
fi


## DEFAULT DATA TYPE
dt="ENCODE"

cd $in_dir   
for d in */
    do 
      if [ "$data_type" == "$dt" ]
         then
           # encode
           folderIN=(${d//_/ })
           data_set=${folderIN[0]%/*}
           tf_name=${folderIN[2]%/*}
      else
         #geo
         folderIN=(${d//./ })
         data_set=${folderIN[0]%/*}
         tf_name=${folderIN[1]%/*}
      fi
       
      #in case the TF name is composed with some specificity of the data set
      tf_nameIN=(${tf_name//-/ })
      if [ "${#tf_nameIN[@]}" -gt "1" ]
         then
           tf_name=${tf_nameIN[-1]}
      fi

      #make it upper case to avoid duplicated folders
      tf_name="$(echo "$tf_name" | tr '[:lower:]' '[:upper:]')"
       echo "TF name: $tf_name"

       #create folder structure 
       if [ ! -d "$out_dir/$tf_name" ]
         then
            mkdir -p $out_dir/$tf_name
       fi

       #move the peaks accordingly
       mv $d $out_dir/$tf_name    
done  
