#!/bin/bash
# NOTE: this code was not parallelized due to the fact that the sequence scoring code is very fast
#usage: bash peak_processing_pwm.sh /path/to/the/input/folder /path/to/the/pwm/scoring/code/pwm_searchPFF /path/to/folder/containing/the/JASPAR/PWMs /path/to/the/whole/genome/hg38/fasta/file /path/to/the/output/folder


## ARGUMENT PARSING
in_dir=;
if [ -z "$1" ]
    then
        echo "No input folder specified. Exiting...";
        exit 1
    else
        in_dir=$1
fi

seq_scoring_code=;
if [ -z "$2" ]
    then
        echo "No sequence scoring code file pecified. Exiting...";
        exit 1
    else
        seq_scoring_code=$2
fi

pwm_folder=;
if [ -z "$3" ]
    then
        echo "No folder containing the PWMs was specified. Exiting...";
        exit 1
    else
        pwm_folder=$3
fi

hg_gen_file=;
if [ -z "$4" ]
    then
        echo "No reference genome file specified. Exiting...";
        exit 1
    else
        hg_gen_file=$4
fi

result_folder=;
if [ -z "$5" ]
    then
        echo "No result folder specified. Exiting...";
        exit 1
    else
        result_folder=$5
fi

mkdir -p $result_folder

#counters
no_pwm_counter=0

# get ready to process
cd $in_dir
#main loop over input folders
for bed in $(find ./$d -name "*.bed")
   do
   
   echo "Processing file $bed"
 
   #create folder structure and output file name
   out_dir=$result_folder
   out_file=$out_dir/$(basename $bed)

   fileIN=$(basename $bed)
   fileIN=(${fileIN//_/ })
   pwm_id=${fileIN[0]%/*}
   pwm=$pwm_folder/${pwm_id}.pwm

   
   fileIN=${fileIN[1]%/*}
   fileIN=(${fileIN//./ })
   tf_name=${fileIN[0]%/*}

   #check if the pwm exists
    if [ ! -e $pwm ]
       then
          echo "No pwm found for file $bed"
          ((no_pwm_counter++))
          exit 1
    fi
    
  # Map back to the genome and get a FASTA file ##################
  getFastaFromBed -s -fo ${out_file}.fa -fi $hg_gen_file -bed $bed

  # Scoring the sequences and output ###################
  $seq_scoring_code $pwm ${out_file}.fa 0.0 -b -n $tf_name > ${out_file}.score
      # 1. Chromosome:start-end
      # 2. Description
      # 3. TF name
      # 4. Strand
      # 5. Absolute score
      # 6. Relative score
      # 7. Start offset of the TFBS in the sequence from the window start (Note it is one based, so no need to +1)
      # 8. End of TFBS relative to start offset (8th column - 7th column gives the TFBS sequence length)
      # 9. Sequence (Note: that it is 1 based, so there is no need to +1 the start offset of the TFBS 

      ################# 6. Format file for centrimo and visualization tools ####################
      awk -v OFS='\t' '{print $1,$3,$5 ";" $6,$4}' ${out_file}.score | awk -v OFS='\t' -F":" '$1=$1' | awk -v OFS='\t' '{sub(/\-/,"\t",$2)};1' > $out_file
      # Output format:
      # 1. Chromosome
      # 2. Start of the TFBS on chromosome (chrStart from step 4 + column 7 from step 4)
      # 3. End of the TFBS on chromosome (resulting chrStart + length of TFBS (8th column - 7th column from step 4)
      # 4. TF name
      # 5. Absolutescore;Relative score
      # 6. Strand

      #Replace the strand with the one from the FASTA header
      awk -v OFS='\t' '{print $1,$2,$3,$4,$5,substr($3,length($3)-1,1)}' $out_file > ${out_file}.tmp
      
      #Remove the strand information from the chromosome end column (e.g., (+))
      sed "s/[(][^)]*[)]//g" ${out_file}.tmp > $out_file
      rm ${out_file}.tmp

     # Clean up 
#     rm ${out_file}.fa
#     rm ${out_file}.score

done #end for over folders

#ouput processing summary
echo "--------------------------------------------------------------------------------------"
echo "$no_pwm_counter data set(s) with no PWM(s)"
echo "--------------------------------------------------------------------------------------"
