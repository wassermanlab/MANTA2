This file contains all the processing steps for the MANTA2 data preparation. Note that this file contains indicative commands, some parts being just pure pseudo-code.

## Configuration
Set the paths below for data processing.

HG38_GENOME_FILE="/path/to/the/whole/reference/genome/fasta/file"

JASPAR_VERTEBRATE_PROFILES="/path/to/where/all/jaspar/core/vertebrates/profiles/are/stored"
JASPAR_PWM_FOLDER="/path/to/where/all/jaspar/pwms/needed/are/stored"

INITIAL_DATA_FOLDER="/path/to/the/original/BED/files"
RAW_DATA_FOLDER="/path/to/storage/of/unzipped/bed/files"
SORTED_BED_FOLDER="/path/to/where/BED/data/ready/for/processing/should/be"

ENCODE_PEAKS_FOLDER="/path/to/the/ENCODE/peaks"
GEO_PEAKS_FOLDER="/path/to/the/GEO/peaks"

PEAKS_FOLDER="/path/to/peak/storing" #temporary; folders will be moved at step 3
SPLIT_PEAKS_FOLDER="/path/to/where/peak/data/ready/for/processing/should/be"
PROCESSED_PEAKS_FOLDER="/path/to/where/the/processed/peaks/should/be/stored"
RESCORED_TFBS_FOLDER="/path/to/where/the/rescored/peaks/should/be/stored"

INTERSECTED_FOLDER="/path/to/where/the/result/of/intersected/peaks/and/tfbs/should/be/stored"
AGGREGATED_FOLDER="/path/to/where/the/aggregated/results/from/the/intersection/should/be/stored"
SORTED_RESULT_FOLDER="/path/to/wehre/the/final/results/are/to/be/stored"

## Intersect JASPAR w/ ReMap

### Step 1: Download JASPAR TFBS predictions
JASPAR genome-wide TFBS predictions as tab-separated values (TSV) files for the hg38 assembly can be downloaded as follows: `bash get_JASPAR_TFBS_predictions.sh $JASPAR_PWM_FOLDER`.

### Step 2: Reformat TSV to BED
JASPAR genome-wide TFBS predictions as tab-separated values files for the hg38 assembly can be downloaded as follows: `bash get_JASPAR_TFBS_predictions.sh`. Note that the bash script will create a folder (*i.e.*  `./JASPAR_TFBS_predictions`) in which to store the downloaded profiles.

# unzip and sort the initial BED files (from Oriol)
#
#dependencies: sort_launcher.sh
#              sort_bed.sh
# command: 
bash sort_launcher.sh $INITIAL_DATA_FOLDER sort_bed.sh $RAW_DATA_FOLDER $SORTED_BED_FOLDER
 
#### STEP 2 ######
# copy the associated datasets (only Remap 2018) 
#
#dependencies: copy_folders_based_on_file.sh
#              encode_datasets_for_manta.txt
#              geo_datasets_for_manta.txt
#commands:
bash copy_folders_based_on_file.sh $ENCODE_PEAKS_FOLDER encode_datasets_for_manta.txt ENCODE $PEAKS_FOLDER/ENCODE

bash copy_folders_based_on_file.sh $GEO_PEAKS_FOLDER geo_datasets_for_manta.txt GEO $PEAKS_FOLDER/GEO

#### STEP 3 #####
# split the data in folders based on the TF name
# 
# dependencies: split_by_tf_name.sh
#
#command: 
bash split_by_tf_name.sh $PEAKS_FOLDER/ENCODE ENCODE $SPLIT_PEAKS_FOLDER

bash split_by_tf_name.sh $PEAKS_FOLDER/GEO GEO $SPLIT_PEAKS_FOLDER

#### STEP 4 ######
# launch processing of peak files in parallel:
# i) concatenated all bed in a TF folder
# ii) sorted the resulting file
# iii) merge the peaks
#
#dependencies: process_peak_files.sh
#               process_launcher.sh 
#
#command:               
bash process_launcher $SPLIT_PEAKS_FOLDER process_peak_files.sh $PROCESSED_PEAKS_FOLDER

#### STEP 5 #####
# rescore the MANTA TFBSs using a modified version of the scoring algorithm: the backward score is ignored
#
# dependencies: sequence_scoring_pwm.sh
#               pwm_search.h
#               pwm_searchPFF
#
# command:
bash sequence_scoring_pwm.sh $SORTED_BED_FOLDER pwm_searchPFF $JASPAR_PWM_FOLDER $HG38_GENOME_FILE $RESCORED_TFBS_FOLDER

#### STEP 6 #####
#intersect all the coresponding peaks with each of the entries in each MANTA TFBS files. 
#
#dependencies: get_intersection.sh
#
#command:
bash get_intersection_tfbs.sh $RESCORED_TFBS_FOLDER $PEAKS_FOLDER/ENCODE $PEAKS_FOLDER/GEO $INTERSECTED_FOLDER

#### STEP 7 #####
# aggregate each of the resulting files to put all intersecting datasets on one line
#
# dependencies: aggregate_datasets.R
#
# command:
R
source("aggregate_datasets.R")
aggregate_datasets(input_folder = $INTERSECTED_FOLDER, output_folder = $AGGREGATED_FOLDER)
q()

#### STEP 8 #####
#sort the result files
#
#dependencies: sort_bed_simple.sh
#
#command:
bash sort_bed_simple.sh $AGGREGATED_FOLDER $SORTED_RESULT_FOLDER

