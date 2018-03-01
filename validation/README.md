This folder contains the code for the `external validation of MANTA2 on allelic imbalance data`.

### Step 1: Download allelic imbalance data
Allelic imbalance data from TF binding events at heterozygous locations are distributed as part of the **Supplementary Data** [here](https://academic.oup.com/nar/article/44/21/10106/2628006).

### Step 2: Get ASB coordiantes in BED format
`grep -v cell ./gkw691_Supp/ASB_GM12878_HeLa_1based.txt | cut -f 1,2 | awk '{print $1 "\t" ($2 - 1) "\t" $2}' > coordinates.hg19.bed`

### Step 3: LiftOver from hg19 to hg38
This was performed using the [liftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver) from the UCSC genome browser. The file `coordinates.hg19.bed` was provided as input and the output file from the liftOver was saved as `coordinates.hg38.liftOver.bed`. Note that **conversion failed on 12 records**. After liftOver, failed conversions were commented (*i.e.* lines starting with `#` in the `coordinates.hg19.bed` file).

### Step 4: Create BED file for MANTA2
`./makeBED.py > ASB.hg38.liftOver.bed`

### Step 5: Execute MANTA2
`../search_manta2.py -d manta2 -H manta.cmmt.ubc.ca -u manta_r -p mantapw -i ASB.liftOver.hg38.bed -t bed > MANTA2.out`

### Step 6: Evaluate MANTA2 impact scores on allelic imbalance data
This step creates a CSV file and two figures in png and svg format, respectively.

`./evaluate.py`
