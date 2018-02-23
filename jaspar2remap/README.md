This folder encloses all the scripts to process the JASPAR TFBS predictions and ReMap ChIP-seq datasets in preparation for MANTA2. Note that the commands included within this `README` file are merely indicative, with some parts being pure pseudo-code.

## Configuration
Set the paths below for data processing.

`HG38_GENOME_FILE="/path/to/the/whole/reference/genome/fasta/file"`, `JASPAR_PWM_FOLDER="/path/to/where/all/jaspar/pwms/needed/are/stored"`,
`INITIAL_DATA_FOLDER="/path/to/the/original/BED/files"`,
`RAW_DATA_FOLDER="/path/to/storage/of/unzipped/bed/files"`,
`SORTED_BED_FOLDER="/path/to/where/BED/data/ready/for/processing/should/be"`,
`ENCODE_PEAKS_FOLDER="/path/to/the/ENCODE/peaks"`,
`GEO_PEAKS_FOLDER="/path/to/the/GEO/peaks"`,
`PEAKS_FOLDER="/path/to/peak/storing"` (temporary; folders will be moved at step 3),
`SPLIT_PEAKS_FOLDER="/path/to/where/peak/data/ready/for/processing/should/be"`,
`PROCESSED_PEAKS_FOLDER="/path/to/where/the/processed/peaks/should/be/stored"`,
`RESCORED_TFBS_FOLDER="/path/to/where/the/rescored/peaks/should/be/stored"`,
`INTERSECTED_FOLDER="/path/to/where/the/result/of/intersected/peaks/and/tfbs/should/be/stored"`,
`AGGREGATED_FOLDER="/path/to/where/the/aggregated/results/from/the/intersection/should/be/stored"`,
`SORTED_RESULT_FOLDER="/path/to/wehre/the/final/results/are/to/be/stored"`

## Intersect JASPAR TFBS predictions w/ ReMap ChIP-seq regions

### Step 1: Get TFBS predictions
For each profile in the JASPAR CORE vertebrates collection: 1) download TFBS predictions on the hg38 genome assembly as tab-separated values files; 2) reformat the TSV file to BED; and 3) compress the resulting BED file.

`bash get_JASPAR_TFBS_predictions.sh $INITIAL_DATA_FOLDER`

### Step 2: Get ChIP-seq datasets
Marius is to provide the scripts that download the data, and provide a description.

`STEP 2
copy the associated datasets (only Remap 2018) 

dependencies: copy_folders_based_on_file.sh
              encode_datasets_for_manta.txt
              geo_datasets_for_manta.txt
commands:
bash copy_folders_based_on_file.sh $ENCODE_PEAKS_FOLDER encode_datasets_for_manta.txt ENCODE $PEAKS_FOLDER/ENCODE

bash copy_folders_based_on_file.sh $GEO_PEAKS_FOLDER geo_datasets_for_manta.txt GEO $PEAKS_FOLDER/GEO
`
### Step 3: Unzip and sort the initial BED files

`bash sort_launcher.sh $INITIAL_DATA_FOLDER sort_bed.sh $RAW_DATA_FOLDER $SORTED_BED_FOLDER`

### Step 4: Split the ChIP-seq data per TF

`bash split_by_tf_name.sh $PEAKS_FOLDER/ENCODE ENCODE $SPLIT_PEAKS_FOLDER
bash split_by_tf_name.sh $PEAKS_FOLDER/GEO GEO $SPLIT_PEAKS_FOLDER`

### Step 5: Intersect TFBS predictions with ChIP-seq regions
For each ChIP-seq region: 1) concatenated all bed in a TF folder; 2) sorted the resulting file; and 3) merge the peaks. The intesected TFBS predictions conform the MANTA database. <--- I don't get this step, could you please check the wording?

`bash process_launcher $SPLIT_PEAKS_FOLDER process_peak_files.sh $PROCESSED_PEAKS_FOLDER`

### Step 6: Re-score the TFBS predictions
This step re-scores the MANTA TFBS predictions using a modified version of the scoring algorithm (the backward score is ignored).

`bash sequence_scoring_pwm.sh $SORTED_BED_FOLDER pwm_searchPFF $JASPAR_PWM_FOLDER $HG38_GENOME_FILE $RESCORED_TFBS_FOLDER`

### Step 7: Intersect
Please write description: Intersect ChIP-seq peaks with each of the entries in each MANTA TFBS files. 

`bash get_intersection_tfbs.sh $RESCORED_TFBS_FOLDER $PEAKS_FOLDER/ENCODE $PEAKS_FOLDER/GEO $INTERSECTED_FOLDER`

### Step 8: DESCRIPTION
Please write description: Aggregate each of the resulting files to put all intersecting datasets on one line

`R
source("aggregate_datasets.R")
aggregate_datasets(input_folder = $INTERSECTED_FOLDER, output_folder = $AGGREGATED_FOLDER)
q()`

### Step 9: Sort blablabla
Please write description.

`bash sort_bed_simple.sh $AGGREGATED_FOLDER $SORTED_RESULT_FOLDER`

