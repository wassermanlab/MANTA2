#usage bash get_intersection_tfbs.sh /path/to/the/manta/rescored/tfbs /path/to/the/processed/ENCODE/peaks /path/to/the/GEO/peaks /path/to/output/folder

## ARGUMENT PARSING
in_dir=;
if [ -z "$1" ]
    then
        echo "No input folder specified. Exiting...";
        exit 1
    else
        in_dir=$1
fi

encode_peaks=;
if [ -z "$2" ]
    then
        echo "No input folder containing encode peaks specified. Exiting...";
        exit 1
    else
        encode_peaks=$2
fi

geo_peaks=;
if [ -z "$3" ]
    then
        echo "No input folder containing geo peaks specified. Exiting...";
        exit 1
    else
        geo_peaks=$3
fi

out_dir=;
if [ -z "$4" ]
    then
        echo "No result folder specified. Exiting...";
        exit 1
    else
        out_dir=$4
fi

mkdir -p $out_dir
cd $in_dir

for bed in $(find ./ -name "*.bed")
   do
     
     echo "Processing file: $bed"

     #get TF name to match it on
     IN=${bed##*/}
     folderIN=(${IN//_/ })
     manta_tf_name=${folderIN[1]%/*}
     folderIN=(${manta_tf_name//./ })
     manta_tf_name=${folderIN[0]%/*}
     manta_tf_name="$(echo "$manta_tf_name" | tr '[:lower:]' '[:upper:]')"

    #special case for NKX3
    if [[ "$manta_tf_name" == "NKX3-1" ]]
        then
           manta_tf_name="NKX3_1"
    fi

    manta_tf_name_alt=""
    #special case for STAT5
    if [[ "$manta_tf_name" == "STAT5A" ]]
        then
           manta_tf_name_alt="STAT5B"
    fi

           matches=()
          #go through all encode peaks
          for enc_peak in $(find $encode_peaks -name "*.narrowPeak")
             do
               #get TF name to match it on
               IN=${enc_peak##*/}
               folderIN=(${IN//./ })
               tf_name=${folderIN[1]%/*}
               tf_name="$(echo "$tf_name" | tr '[:lower:]' '[:upper:]')"

              if [[ "$tf_name" == "$manta_tf_name" ]] || [[ "$tf_name" == "$manta_tf_name_alt" ]]
                 then
                     matches+=($enc_peak)
              fi
          done

          #go through all encode peaks
          for geo_peak in $(find $geo_peaks -name "*.narrowPeak")
             do
               #get TF name to match it on
               IN=${geo_peak##*/}
               folderIN=(${IN//./ })
               tf_name=${folderIN[1]%/*}
               tf_name="$(echo "$tf_name" | tr '[:lower:]' '[:upper:]')"
 
              if [[ "$tf_name" == "$manta_tf_name" ]] || [[ "$tf_name" == "$manta_tf_name_alt" ]]
                 then
                     matches+=($geo_peak)
                   
              fi
          done
          echo "Matches found ${matches[@]}"
          intersectBed -wa -wb -filenames -f 1 -a $bed -b ${matches[@]} > $out_dir/${bed##*/}
          echo "---------------------------------------------------------------------------------"
done

