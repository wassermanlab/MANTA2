# MANTA2
The **M**ongoDB for the **AN**alysis of **T**ranscription factor (TF)-binding site (TFBS) **A**lterations (MANTA) was originally created in 2015 to study the impact of regulatory mutations in B-cell lymphomas [(Mathelier et al. 2015)](https://doi.org/10.1186/s13059-015-0648-7). The database stores TFBSs predicted in the human genome by combining ChIP-seq regions from [ReMap](http://remap.cisreg.eu) with [JASPAR](http://jaspar.genereg.net) TF binding profiles, as well as the potential impact scores for all possible single nucleotide variants (SNVs) at these TFBSs. This second release of the database, MANTA2, houses >48 million TFBS predictions for 225 TFs and covers ~8% of the human genome ([hg38](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26)).

## Content
The repository is organized as follows:
* The `examples` folder contains a shell script (*i.e.* `get_VCF_example.sh`) to generate a variant file in VCF format
* The `jaspar2remap` folder contains the scripts that were used to merge JASPAR TFBS predictions with ReMAP ChIP-seq peaks along with some instructions on how to use them
* The `manta2` folder contains the scripts related to loading the data to the MongoDB system and the web interface of MANTA2
* The `snv_computation` folder contains the scripts used to compute the impact scores of the SNVs within the MANTA2 TFBSs and instructions on how to use them
* The symbolic link to the `search_manta2.py` script provides programmatic access to MANTA2 (both locally and hosted at the Wasserman Lab)

## Dependencies
MANTA2 requires the following dependencies:
* [`MongoDB`](https://www.mongodb.com) (≥2.4; version 3.4 or higher is strongly recommended)
* `Python` (>2.7) with the [`Flask`](http://flask.pocoo.org) and [`PyMongo`](https://api.mongodb.com/python/current/) libraries

## Installation
To create a local build of MANTA2, follow the next steps:
* Install (if necessary) and initialize the MongoDB system: instructions for the main OS can be found [here](https://docs.mongodb.com/manual/administration/install-community/)
* Download, uncompress and unpack the MANTA2 [MongoDB dump](https://doi.org/10.5281/zenodo.1044747): this creates a folder (*i.e.* `manta2_mongodb_dump`) containing 5 different files: `experiments.bson`, `experiments.metadata.json`, `system.indexes.bson`, `tfbs_snvs.bson`, and `tfbs_snvs.metadata.json`
* Restore MANTA2 from the MongoDB dump: `mongorestore -d $DB_NAME $PATH_TO_MONGODB_DUMP`, where `$DB_NAME` specifies the database name for MANTA2 in the MongoDB system (e.g. `manta2`) and `$PATH_TO_MONGODB_DUMP` the path to the folder `manta2_mongodb_dump` from the previous step
* Create a user with `read` privileges to the MANTA2 database: `mongo --eval "db.getSiblingDB('$DB_NAME').createUser({user: '$USER', pwd: '$PASSWORD', roles: [{role: 'read', db: '$DB_NAME'}]})"`

Note that, in the commands from the previous steps, the strings `$DB_NAME`, `$PATH_TO_MONGODB_DUMP`, `$USER`, `$PASSWORD` should be replaced with your own.

To check whether MANTA2 was restored successfully, type `mongo` on a terminal to open a MongoDB shell and then type `show databases`. The system should list all available Mongo databases, including MANTA2 (*e.g.* `manta2  14.616GB`).

## Usage
The script `search_manta2.py` provides programmatic access to MANTA2. It requires the following inputs:
* The name of the MANTA2 database in the MongoDB system (option `-d`)
* The name of the server where the MongoDB system is hosted (option `-H`)
* A user with `read` privileges to the MANTA2 database (option `-u`)
* The password for the previous user (option `-p`)
* A file containing a list of variants in [VCF](https://genome.ucsc.edu/FAQ/FAQformat.html#format10.1), [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) or [GFF](https://genome.ucsc.edu/FAQ/FAQformat.html#format3) format (option `-i`)

Non-mandatory options include:
* The format of the input variant file (option `-t`; by default the script tries to identify the input format automatically)
* The name of a file to output the results (option `-o`; by default is set to the standard output stream (stdout))

As a usage example, the MANTA2 database hosted at the Wasserman Lab can be accessed as follows: `./search_manta2.py -d manta2 -H manta.cmmt.ubc.ca -u manta_r -p mantapw -i <variant file>`.

A variant file can be obtained by executing the shell script `get_VCF_example.sh` located in the `./examples/` folder. The resulting VCF file (*i.e.* `chr20.vcf`) contains high-confidence SNP, small indel, and homozygous reference calls on chromosome 20 from the Genome in a Bottle (version 3.3.2) sample HG001 [(Zook et al. 2014)](https://doi.org/10.1038/nbt.2835).

The `search_manta2.py` script returns all TFBS predictions potentially impacted by these variants as tab-separated values. For each TFBS alteration, the script provides the variant information along with the associated wild-type (reference) and mutated (alternative) TFBS information, including:
* the chromosome and position of the variant;
* the reference and alternative alleles at that genomic location;
* the mutation ID (if the input file format allowed for it, otherwise the field is displayed as `.`);
* the TF name and associated JASPAR profile ID;
* the start, end and strand, as well as the absolute (raw) and relative scores for both the reference and alternative TFBSs; and
* the impact score.

Users planning on performing large numbers of searches should create their local builds of the MANTA2 database (see [Installation](https://github.com/wassermanlab/MANTA2/blob/master/README.md#installation)).

The MANTA2 database hosted at the Wasserman Lab can also be accessed *via* a dedicated web server at [URL](URL). Similar to the `search_manta2.py` script, the server requires as input a list of variants in [VCF](https://genome.ucsc.edu/FAQ/FAQformat.html#format10.1), [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) or [GFF](https://genome.ucsc.edu/FAQ/FAQformat.html#format3) format, and it returns all TFBS predictions potentially impacted by these variants as a tab-separated values table. The table can be sorted on any column by clicking on the column header.
