This folder contains the code to compute the `impact` of all possible SNVs that could occur in a set of predicted TFBSs using a JASPAR position frequency matrix (PFM).

## Dependencies
The scripts present in this folder have been developed to be used in a UNIX OS.
* `awk`
* The [`BEDtools`](http://bedtools.readthedocs.io) accessible from the PATH
* `Python` version 2.7 and [`Biopython`](http://biopython.org)
* `Perl` and the following modules:
..* [`TFBS Perl module`](http://tfbs.genereg.net)
..* [`BioPerl`](http://bioperl.org)
..* [`Statistics-R`](http://search.cpan.org/~fangly/Statistics-R-0.34/lib/Statistics/R.pm)

## Configuration
Set the path to the directory containing the scripts in the file:

`compute_snv_impacts_launcher_allalt.sh`.

## Calculate SNVs
Please refer to the test case script `./test/test.sh` for further details on how to launch the computations.
