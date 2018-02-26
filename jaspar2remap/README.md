This folder contains the scripts to process the JASPAR 2018 TFBS predictions and ReMap 2018 ChIP-seq peaks in preparation for the MANTA2 database. **Note** that the steps/commands included within this `README.md` file are merely indicative, with some parts being pure pseudo-code.

## Dependencies
The scripts present in this folder have been developed to be used in a UNIX OS (for Microsoft Windows, the system should include support for Bash Unix Shell).
* The bash commands `gzip`, `tr`, `awk` and `sed`
* `R` statistical environment (>3.0)

## Requirements 
 As software dependencies, the system should have support for . It should have the . The base packages should suffice.

## Configuration
Set the paths below for data processing. The scripts will produce output at most of the steps listed below. **Important**: ensure that the **genome file** comes from the `hg38` human genome assembly.

```
HG38_GENOME_FILE="/path/to/the/whole/reference/genome/fasta/file"
JASPAR_PWM_FOLDER="/path/to/where/all/jaspar/pwms/needed/are/stored"
```

The JASPAR TFBSs
```
INITIAL_DATA_FOLDER="/path/to/the/original/BED/files"
RAW_DATA_FOLDER="/path/to/storage/of/unzipped/bed/files"
SORTED_BED_FOLDER="/path/to/where/BED/data/ready/for/processing/should/be"
```

The ChIP-seq peaks
```
ENCODE_PEAKS_FOLDER="/path/to/the/ENCODE/peaks"
GEO_PEAKS_FOLDER="/path/to/the/GEO/peaks"
```

Processed data
```
RESCORED_TFBS_FOLDER="/path/to/where/the/rescored/peaks/should/be/stored"
INTERSECTED_FOLDER="/path/to/where/the/result/of/intersected/peaks/and/tfbs/should/be/stored"
AGGREGATED_FOLDER="/path/to/where/the/aggregated/results/from/the/intersection/should/be/stored"
SORTED_RESULT_FOLDER="/path/to/wehre/the/final/results/are/to/be/stored"
```

## Intersect JASPAR TFBS predictions w/ ReMap ChIP-seq regions

### Step 1: Get the TFBS predictions
For each profile in the JASPAR 2018 CORE vertebrates collection: 1) download the TFBS predictions on the hg38 genome assembly as a tab-separated values (TSV) file; 2) reformat the TSV file to BED; and 3) compress the resulting BED file.

`bash get_JASPAR_TFBS_predictions.sh $INITIAL_DATA_FOLDER`

### Step 2: Get the ChIP-seq peak datasets
Retrieve the collection of ChIP-seq peaks from the ReMap 2018 resource from [here](http://tagc.univ-mrs.fr/remap/download/MACS/ReMap2_allPeaks.bed.gz).

### Step 3: Copy the associated datasets
After decompressing the ReMap 2018 ChIP-seq peaks, copy in a separate folder only the ChIP-seq peak datasets associated with the JASPAR TFBSs. **Note** that this step assumes that the ENCODE peaks and the Public (GEO & ArrayExpress) peaks are located in different folders. 

```
bash copy_folders_based_on_file.sh $ENCODE_PEAKS_FOLDER encode_datasets_for_manta.txt ENCODE $PEAKS_FOLDER/ENCODE
bash copy_folders_based_on_file.sh $GEO_PEAKS_FOLDER geo_datasets_for_manta.txt GEO $PEAKS_FOLDER/GEO
```

### Step 4: Unzip and sort the initial BED files (JASPAR TFBSs)

`bash sort_launcher.sh $INITIAL_DATA_FOLDER sort_bed.sh $RAW_DATA_FOLDER $SORTED_BED_FOLDER`

### Step 5: Re-score the TFBS predictions
TFBS predictions are re-scored using a modified version of the algorithm (original scores provided in the BED files are ignored).

`bash sequence_scoring_pwm.sh $SORTED_BED_FOLDER pwm_searchPFF $JASPAR_PWM_FOLDER $HG38_GENOME_FILE $RESCORED_TFBS_FOLDER`

### Step 6: Intersect the JASPAR TFBSs with the ReMap 2018 ChIP-seq peaks
Each TFBS file is intersected with **all** corresponding ReMap 2018 ChIP-seq peaks.

`bash get_intersection_tfbs.sh $RESCORED_TFBS_FOLDER $PEAKS_FOLDER/ENCODE $PEAKS_FOLDER/GEO $INTERSECTED_FOLDER`

### Step 7: Aggregating the intersected regions
The results of the intersection are post-processed in order to group the entries by TFBS. In other words, all entries corresponding to the same TFBS are aggregated such that all the datasets where an intersection occurred for that TFBS will be put in one single line, with the name of the datasets separated by a comma.

```
R
source("aggregate_datasets.R")
aggregate_datasets(input_folder = $INTERSECTED_FOLDER, output_folder = $AGGREGATED_FOLDER)
q()
```

### Step 8: Sort the results of the aggregation - optional
The resulting BED files can be sorted by chromosome name and chromosome start (this step is optional).

`bash sort_bed_simple.sh $AGGREGATED_FOLDER $SORTED_RESULT_FOLDER`
