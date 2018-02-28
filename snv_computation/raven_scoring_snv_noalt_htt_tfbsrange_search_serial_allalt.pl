#!/bin/env perl

=head1 DESCRIPTION

For each variation in the input VCF file, compute which TFBSs are affected
by the variation and output the results. TFBS profile is specified
by a matrix file (PFM)

=cut

use warnings;
use strict;
use threads;

use Getopt::Long;
use Pod::Usage;
use TFBS::Matrix::PFM;
use TFBS::DB::JASPAR5;
use Bio::SeqIO;
use Statistics::R;
use POSIX;
use Time::HiRes qw( gettimeofday );

use constant MINIMAL_TFBS_THRESHOLD     => "10%";
use constant MINIMAL_TFBS_THRESHOLD_NUM => "0.10";
use constant TFBS_THRESHOLD             => "85%";
use constant TFBS_SCORE_DIFF            => "0%";

use Parallel::Iterator qw/iterate_as_array/;

my $start_time = time();

my $src;
my $in_file;
my $matrix_file;
my $out_dir;
my $out_file;
GetOptions(
    'i=s'  => \$in_file,
    'm=s'  => \$matrix_file,
    'o=s'  => \$out_dir,
    'of=s' => \$out_file,
    'c=s'  => \$src,
);

unless ($out_dir) {
    $out_dir = ".";
}
unless ( -e $out_dir ) {
    mkdir "$out_dir";
}

#
# Three directories used to store/pass information between threads. Do not use these directories
# for other purposes. These directories may all be made children of another temporary directory,
# which can be specified by the user (has yet to be implemented). Consider using a database instead.
# SQLite is a favourable option, for its speed compared to MySQL.
#
my $out_dir_scoring = "$out_dir/tempScores";
my $temp_matrix_dir = "$out_dir/tempMatrix";
my $temp_scores_dir = "$out_dir/allScoresDir";
unless ( -e $out_dir_scoring ) {
    mkdir "$out_dir_scoring";
}
unless ( -e $temp_matrix_dir ) {
    mkdir "$temp_matrix_dir";
}
unless ( -e $temp_scores_dir ) {
    mkdir "$temp_scores_dir";
}
unless ($in_file) {
    pod2usage(
        -verbose => 1,
        -msg     => "No input variation file specified\n"
    );
}

unless ($matrix_file) {
    pod2usage(
        -verbose => 1,
        -msg     => "Please specify a matrix file name (-m)\n"
    );
}

my $threshold_str      = TFBS_THRESHOLD;
my $min_score_diff_str = TFBS_SCORE_DIFF;

my $threshold;
my $threshold_type = 'absolute';
if ( $threshold_str =~ /(\S+)%/ ) {
    $threshold      = $1 / 100;
    $threshold_type = 'relative';
}
else {
    $threshold = $threshold_str;
}

my $min_score_diff;
my $score_diff_type = 'absolute';
if ( $min_score_diff_str =~ /(\S+)%/ ) {
    $min_score_diff  = $1 / 100;
    $score_diff_type = 'relative';
}
else {
    $min_score_diff = $min_score_diff_str;
}

print "\nFetching matrix set...\n";

my $matrix_set = get_matrix_set();
unless ( $matrix_set && $matrix_set->size > 0 ) {
    die "Error fetching matrix set\n";
}

my $max_profile_width = max_profile_width($matrix_set);

print "\nChecking VCF file...\n";

check_VCF($in_file);

print "Searching for TFBSs affected by variations...\n";

vcf_search_tfbss( $in_file, $matrix_set );

print "Done.\n\n";
my $end_time = time();
print "Time taken: ", $end_time - $start_time, "\n";

exit;

#
# Process VCF file one line at a time (to save memory).
#
# For each line in the file, extract the corresponding reference sequence with
# enough flanking nucleotides to search all TFBS matrices for each position
# where the matrix overlaps the variant nucleotides. Scan both the reference
# sequence and the variant sequence. In cases where the matrix scores above
# threshold for one sequence and below threshold for the other and the score
# difference is >= than the minimum score difference, write this information
# to the output file.
#
# Keep in mind that the min_score_diff, as well as chromosome M positions, are not
# supported yet. Chromosome M positions seem to induce errors upon sequence
# extraction.
#
sub vcf_search_tfbss {
    my ( $file, $matrix_set ) = @_;

    my $timeadp0 = gettimeofday();

    open( FH, $file ) || die "Error opening input VCF file $file\n";

    my %tf_names;
    my %tab_scores;
    my %tab_info;

    my $timeread0     = gettimeofday();
    my $matrixfilestr = "";
    my $matrixIDstr   = "";
    my $startsitestr  = "";
    my $endsitestr    = "";
    my $miter         = $matrix_set->Iterator();

    #
    # Precomputation before iterating through the file. Necessary to
    # prevent unnecessary conditional checks within the critical portion,
    # as well as to supply parameters for all_scores_two.
    #
    while ( my $matrix = $miter->next ) {
        my $tf_name = $matrix->ID;
        $tf_names{$tf_name} = 1;
        unless ( -e "$out_dir/$tf_name" ) {
            mkdir "$out_dir/$tf_name" or die $!;
        }
        unless ( -e "$out_dir/$tf_name/Rel_score" ) {
            mkdir "$out_dir/$tf_name/Rel_score" or die $!;
        }
        my $raw        = $matrix->rawprint;
        my $matrixpath = "$temp_matrix_dir/" . "$tf_name" . ".txt";
        unless ( -e $matrixpath ) {
            open( OUTMATRIX, ">" . $matrixpath );
            print OUTMATRIX ">" . "$tf_name\n" . "$raw";
            close(OUTMATRIX);
        }

        #
        # Computes the reference start and end sites based SOLELY on matrix
        # length --> this means that these are not the actual start or end
        # sites. The values of these will be updated when processing the file.
        #
        my $ref_start_site = 2 - $matrix->length;
        my $ref_end_site   = $matrix->length - 2;
        #
        # Computes information to pass on to the threads, for execution of the
        # all_scores_two procedure. Possibly a better way to implement this.
        # Do NOT use arrays, as they are passed by reference to the individual
        # threads, and the threads each retrieve copies of the original array
        # upon each iteration. Scalars are the way to go. OPTIMIZE ME.
        #
        if ( $matrixfilestr eq "" ) {
            $matrixfilestr = $matrixpath;
            $matrixIDstr   = $tf_name;
            $startsitestr  = "$ref_start_site";
            $endsitestr    = "$ref_end_site";
        }
        else {
            $matrixfilestr .= "," . $matrixpath;
            $matrixIDstr   .= "," . $tf_name;
            $startsitestr  .= ",$ref_start_site";
            $endsitestr    .= ",$ref_end_site";
        }

        system("rm -r $out_dir/$tf_name");
    }
   #
   # Initialize thread info. Allows for batching of the individual scoring jobs.
   # Very crude implementation. OPTIMIZE ME.
   #
    my $outputs          = "";
    my $erroutputs       = "";
    my $chroms           = "";
    my $poss             = "";
    my $rel_poss         = "";
    my $refs             = "";
    my $alts             = "";
    my $alt_top_scores   = "";
    my $ref_top_scores   = "";
    my $ref_seqs         = "";
    my $tf_names         = "";
    my $matrixpaths      = "";
    my $is_indels        = "";
    my $strands          = "";
    my $line_nums        = "";
    my $min_score_diffs  = "";
    my $score_diff_types = "";
    my $out_dirs         = "";
    my $integer          = 0;
#
# Splits into array now, to avoid doing this computation later. Note that this should not,
# and cannot, be done without making a copy of the array later, as scalar values are preserved
# for each worker using the Storable module.
#
    my @startsitestrarr = split ",", $startsitestr;
    my @endsitestrarr   = split ",", $endsitestr;

#
# Iterate over all the lines in the file, producing the reference and alternative
# scores. Meanwhile, batch + submit jobs to the cluster compute nodes. Note that
# batching is necessary due to the fixed time interval after which the cluster can
# dequeue jobs off of the wait queue.
#
# This portion is parallelized on the head node. You can specify more worker threads;
# read the implementation of Parallel::Iterator. The default is 10.
#
# After a particular integer limit, which approximates min(# of (pos,matrix) 2-tuples
# completed by worker i), batching will stop, and sequential submission of jobs will
# resume. This should be changed, so as to allow some thread-to-thread communication.
#
    my $line_num = 0;
    open( OUT, ">" . $out_file );
    while ( my $line = <FH> ) {
        ++$line_num;
        if ( $line_num % 1000 == 0 ) {
            print STDERR $line_num, "\n";
        }
        chomp $line;
        my @elem          = split /\t/, $line;
        my $chrom         = $elem[0];
        my $ref_strand    = $elem[1];
        my $position      = $elem[2];
        my $end_pos       = $elem[3];
        my $tf_file       = $elem[4];
        my $ref_rel_score = $elem[5];
        my $ref_abs_score = $elem[6];
        my $ref_species   = $elem[8];
        my $ref_seq       = $elem[9];
        my $ref_len       = 1;

        if ( $chrom =~ /^chr(\S+)/ ) {
            $chrom = $1;
        }
        my $rel_position     = 2 * $max_profile_width + $ref_len;
        my $seq_start        = $position - $rel_position + 1;
        my $seq_end          = $end_pos + $rel_position - 1;
        my $mismatched_ref   = 0;
        my @startsitestrcopy = @startsitestrarr;
        my @endsitestrcopy   = @endsitestrarr;
        foreach my $startsite (@startsitestrcopy) {
            $startsite += $rel_position - $ref_len;
        }
        foreach my $endsite (@endsitestrcopy) {
            $endsite += $rel_position + $end_pos - $position + $ref_len;
        }
        my $startsitestr2 = join ",", @startsitestrcopy;
        my $endsitestr2   = join ",", @endsitestrcopy;

      #
      # Retrieve a hash of (matrix,array) pairs corresponding to the calculation
      # of all_scores for particular matrices.
      #
        my $scorefile = all_scores_two(
            $matrixfilestr, $ref_seq, $startsitestr2, $endsitestr2,
            $matrixIDstr,   $chrom,   $position,      $end_pos,
            $tf_file,       $ref_strand
        );

        #		print STDOUT $scorefile,"\n";
        my %lookuphash   = all_scores_hash($scorefile);
        my $strand       = 1;
        my $ref_from_seq = substr( $ref_seq, $position - $seq_start, $ref_len );
        my $ref          = $ref_from_seq;

        my $timetot = 0.000000;
        ######
  # Iterate through all matrices in the matrix_set, for the given position. File
  # format may be updated in the future, to avoid a global set of matrices for
  # all positions.
  #
        my $miter = $matrix_set->Iterator();
        while ( my $matrix = $miter->next ) {
            my $tf_name     = $matrix->ID;
            my $contains_tf = 1;
            ( my $mat_filename = $tf_name ) =~ s{^>}{};
            if ( $mat_filename =~ /$tf_file/ ) {
                $contains_tf = 1;
            }
            if ( !$contains_tf ) {
                next;
            }
            ## Diagnostic ##
            my $matrixpath = "$temp_matrix_dir/" . "$tf_name" . ".txt";
            my $ref_start_site;
            my $ref_end_site;
            if ( $strand == -1 ) {
                $ref_start_site =
                  $rel_position - $matrix->length + 2 - $ref_len;
                $ref_end_site = $rel_position + $matrix->length - 1;
            }
            else {
                $ref_start_site = $rel_position - $matrix->length + 1;
                $ref_end_site =
                  $rel_position +
                  $end_pos -
                  $position +
                  $ref_len +
                  $matrix->length - 2;
            }
            my $ref_top_score_raw_point;
            my @ref_top_score_raw_arr;
            my @ref_seq_array;
            my $ref_top_score;
            my @ref_top_score_arr;
            my @ref_top_strand_arr;
            my @ref_top_pos_arr;
            my $ref_top_strand_point;
            my $ref_top_pos_point;
            my @ref_top_rel_score_arr;
            my @ref_abs_pos_arr;
            my @ref_abs_end_arr;
            my $is_indel = 0;

# Indels must be treated differently, due to their incompatibility with the all_scores
# procedure. May be implemented in C++ in the future, to allow for faster indel computation.
#
            unless ($is_indel) {
                $| = 1;
                if ( !$lookuphash{ $matrix->ID } ) {
                    print STDERR "ERROR OCCURRING:\t", $chrom, "\t", $position,
                      "\t", $matrix->ID, "\n";
                }
                @ref_seq_array = @{ $lookuphash{ $matrix->ID } };
            }
            else {
                $ref_top_score = top_score( $matrix, $ref_seq, $ref_start_site,
                    $ref_end_site );
            }
#
# Iterate through all of the alternate alleles for the given position and matrix.
# Indels that are larger than the matrix length are not guaranteed to compute correctly.
#
            my @alt_top_scores;
            my @alt_top_poss;
            my @alt_top_strands;
            my @alt_top_rel_scores;
            my @poss;
            my @ref_nucs;
            my @alt_array;

            foreach my $pos (
                $rel_position .. ( $rel_position + $matrix->length - 1 ) )
            {
                my $ref_nuc = substr( $ref_seq, $pos - 1, 1 );
                if ( $ref_nuc eq "N" ) {
                    next;
                }
                my @alts = ();
                foreach my $nuc ( "A", "C", "G", "T" ) {
                    if ( !( $nuc eq $ref_nuc ) ) {
                        push( @alts, $nuc );
                    }
                }
                foreach my $alt (@alts) {
                    ++$integer;
                    my $alt_len = length $alt;
                    my $alt_seq = $ref_seq;
                    if ( $strand == -1 ) {
                        substr(
                            $alt_seq, $pos - 1 - $ref_len + 1,
                            $ref_len, reverse_complement($alt)
                        );
                    }
                    else {
                        substr( $alt_seq, $pos - 1, $ref_len, $alt );
                    }

# if $alt_len is too big this will fail, obviously, so we restrict $alt_len to be less than $matrix->length for now
#my $timealt0 = gettimeofday();
                    my $alt_top_score;
                    my $alt_top_pos;
                    my $alt_top_strand;
                    unless ($is_indel) {
                        if ( $strand == 1 ) {
                            ( $alt_top_score, $alt_top_strand, $alt_top_pos ) =
                              top_score_search(
                                \@ref_seq_array,
                                $pos - $matrix->length + 1,
                                $pos + $alt_len + $matrix->length - 2,
                                $alt,
                                $ref_nuc,
                                $pos,
                                $matrix
                              );
                            $alt_top_score = sprintf "%.3f", $alt_top_score;
                        }
                        else {
                            ( $alt_top_score, $alt_top_strand, $alt_top_pos ) =
                              top_score_search(
                                \@ref_seq_array,
                                $pos - $matrix->length + 1,
                                $pos + $alt_len + $matrix->length - 2,
                                reverse_complement($alt),
                                reverse_complement($ref_nuc),
                                $pos,
                                $matrix
                              );
                            $alt_top_score = sprintf "%.3f", $alt_top_score;
                        }
                    }
                    push( @alt_top_scores,  $alt_top_score );
                    push( @alt_top_strands, $alt_top_strand );
                    push( @alt_top_poss,    $alt_top_pos );
                    my $alt_top_rel_score =
                      ( $alt_top_score - $matrix->min_score ) /
                      ( $matrix->max_score - $matrix->min_score );
                    push( @alt_top_rel_scores, $alt_top_rel_score );
                    push( @alt_array,          $alt );
                    push( @poss,               $pos );
                    push( @ref_nucs,           $ref_nuc );
                }
            }
            my $alt_score_total = 0;
            my $alt_score_mean;
            my $alt_score_stdev;
            my $alt_counts = 0;
            foreach my $alt_score (@alt_top_scores) {
                $alt_score_total += $alt_score;
                $alt_counts++;
            }
            $alt_score_mean = $alt_score_total / $alt_counts;
            my $alt_score_prevar = 0;
            foreach my $alt_score (@alt_top_scores) {
                $alt_score_prevar += ( $alt_score - $alt_score_mean )**2;
            }
            $alt_score_stdev = ( $alt_score_prevar / ( $alt_counts - 1 ) )**0.5;
            foreach my $i ( 0 .. scalar(@poss) - 1 ) {
                my $alt_pos_start =
                  $alt_top_poss[$i] +
                  $poss[$i] +
                  $position -
                  $rel_position -
                  $matrix->length + 1;
                my $alt_pos_end = $alt_pos_start + $matrix->length - 1;
                my $z_score;
                if ( $alt_score_stdev != 0 ) {
                    $z_score = ( $alt_top_scores[$i] - $alt_score_mean ) /
                      $alt_score_stdev;
                }
                else {
                    $z_score = "N/A";
                }
                my $abs_pos = $poss[$i] - $rel_position + $position;
                print OUT $chrom, "\t", $abs_pos, "\t", $poss[$i], "\t",
                  $ref_nucs[$i], "\t", $alt_array[$i], "\t", $position, "\t",
                  $end_pos, "\t", $ref_abs_score, "\t", $ref_rel_score, "\t",
                  $ref_strand, "\t", $tf_name, "\t", $alt_top_scores[$i], "\t",
                  $alt_top_rel_scores[$i], "\t", $alt_top_strands[$i], "\t",
                  $alt_pos_start, "\t", $alt_pos_end, "\t", $z_score, "\t",
                  $alt_counts, "\t", $ref_species, "\n";
            }
        }

        $erroutputs       = "";
        $chroms           = "";
        $poss             = "";
        $rel_poss         = "";
        $refs             = "";
        $alts             = "";
        $alt_top_scores   = "";
        $ref_top_scores   = "";
        $ref_seqs         = "";
        $tf_names         = "";
        $matrixpaths      = "";
        $is_indels        = "";
        $strands          = "";
        $line_nums        = "";
        $min_score_diffs  = "";
        $score_diff_types = "";
        $out_dirs         = "";
        ## So allScoresDir doesn't grow too big
        system("rm $scorefile");
    }
    close(FH);
    close(OUT);

#
# Collect the outputs from the computation. Interpret these outputs and store them.
#
    system("rm $file");
    system("rm -r $temp_scores_dir");
    system("rm -r $temp_matrix_dir");
    system("rm -r $out_dir_scoring");
    print STDOUT "All done - Lets have a beer.\n";

    #exit;
}

# Read through VCF file and check whether it passes overall criteria for
# analysus. Currently the only check is whether any two variations are within
# the maximum matrix width nucleotides from one another. If so output a warning
# message, as the script does not currently take into account cases where more
# than one variant positions may affect a single TFBS.
#
# NOTE: This check ASSUMES the VCF file is sorted by chromosomal position!!!
#
sub check_VCF {
    my ($file) = @_;

    my $ok = 1;

    open( FH, $file ) || die "Error opening input VCF file $file\n";

    my $last_chrom = ' ';
    my $last_pos   = -1;
    my $line_num   = 0;
    while ( my $line = <FH> ) {
        $line_num++;

        chomp $line;

        my @elem = split /\t/, $line;

        my $chrom    = $elem[0];
        my $position = $elem[2];

        if ( $chrom eq $last_chrom ) {
            if ( $position - $last_pos + 1 <= $max_profile_width ) {
                print
                  "\n\nWARNING: VCF contains at least one case where",
                  " multiple variations may affect a single TFBS.",
                  " This script does not currently take into",
                  " consideration such cases, i.e. only the effect of",
                  " individual variations on a single TFBS are reported.",
                  "\n\n";
                print "$position\t";
                print "$last_pos\n";

                $ok = 0;

                last;
            }
        }

        $last_chrom = $chrom;
        $last_pos   = $position;
    }

    close(FH);

    return $ok;
}

#
# Retrieve the matrix
#
sub get_matrix_set {
    my $matrix_set;

    $matrix_set = read_PWMs($matrix_file);
    die "Error reading PWMs from $matrix_file\n" unless $matrix_set;

    return $matrix_set;
}

#
# Reads in matrices from the specified matrix file. Matrices are in position frequency
# matrix form.
#
sub read_PFMs {
    my ($file) = @_;

    open( FH, $file ) || die "Error opening PFM file $file\n";

    my $matrix_set = TFBS::MatrixSet->new();

    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    while ( my $line = <FH> ) {
        chomp $line;
        next unless $line;
        if ( $line =~ /^>\s*(\S+)/ ) {
            $name = $1;
        }
        else {
            if ( $line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/ ) {

                # line of the form: A [ # # # ... # ]
                $matrix_string .= "$1\n";
            }
            elsif ( $line =~ /^\s*\d+/ ) {

                # line of the form: # # # ... #
                $matrix_string .= "$line\n";
            }
            else {
                next;
            }
            $line_count++;

            if ( $line_count == 4 ) {
                my $pfm = TFBS::Matrix::PFM->new(
                    -matrixstring => $matrix_string,
                    -name         => $name,
                    -ID           => $name
                );

                $matrix_set->add_Matrix( $pfm->to_PWM );

                $line_count    = 0;
                $name          = '';
                $matrix_string = '';
            }
        }
    }
    close(FH);

    return $matrix_set;
}

#
# Computes the maximum matrix width (# of columns) among all of the matrices
# in the matrix set.
#
# NOTE: right now we always return 30 as the max size for a JASPAR PFM is 15
#
sub max_profile_width {
    my ($matrix_set) = @_;

    my $max_width = 0;

    my $iter = $matrix_set->Iterator;
    while ( my $matrix = $iter->next ) {
        if ( $matrix->length > $max_width ) {
            $max_width = $matrix->length;
        }
    }

    return 30;
}
#
# (VERY SLOW) Score calculation function. Give it sequence informations and a pwm matrix,
# it returns an item top_site with all needed score information.
#
sub top_site {
    my ( $matrix, $seq, $search_start, $search_end ) = @_;
    my $site_set = $matrix->search_seq(    # this right here is the trouble
        -seqstring => $seq,
        -threshold => MINIMAL_TFBS_THRESHOLD,
        -subpart   => {
            -start => $search_start,
            -end   => $search_end
        }
    );
    unless ( $site_set && $site_set->size > 0 ) {
        return undef;
    }
    my $top_site;
    my $top_score = -999;
    my $iter      = $site_set->Iterator;
    while ( my $site = $iter->next ) {
        if ( $site->score > $top_score ) {
            $top_score = $site->score;
            $top_site  = $site;
        }
    }
    return $top_site;
}

#
# Computes the reverse, and the complement, of a particular DNA sequence. Should
# be self explanatory.
#
sub reverse_complement {
    my $seq = shift;

    $seq =~ tr/acgtACGT/tgcaTGCA/;

    return reverse $seq;
}

#
# Really simple functions to sort numeric items.
# (They are actually not used in the script but
# i keep them for checking data if necessary)
#
sub par_num { return $a <=> $b }

#
# Keep in mind that this sort will not appreciate NA's in the
# output. Does not affect the output values themselves, only the
# ordering of them. NA's arise from p-value calculations in which
# the reference/alternative is the lowest/highest score.
#
sub numeric_sort {
    if    ( $a < $b )  { return -1; }
    elsif ( $a == $b ) { return 0; }
    else               { return 1; }
}

#
# Prototype procedure, to do indel computations quickly, in C++.
# Not implemented yet, because of certain issues...
#
#sub top_scores_two {
#	my ($matrixfiles, $seq, $startlist, $endlist, $namestr, $chrom, $pos) = @_;
#	`./littleprogramindels $matrixfiles $seq $startlist $endlist $namestr $chrom $pos $temp_scores_dir`;
#	return $temp_scores_dir."/"."$chrom"."$pos".".txt";
#}
#
# Returns a path to a list of arrays for each matrix. Each array stores
# the computation for the score at every position along the sequence,
# with the given matrix.
#
sub all_scores_two {
    my (
        $matrixfiles, $seq, $startlist, $endlist,  $namestr,
        $chrom,       $pos, $end,       $chipfile, $strand
    ) = @_;

#print "littleprogramrange $matrixfiles $seq $startlist $endlist $namestr $chrom $pos $temp_scores_dir\n\n";
#print STDOUT $endlist,"\n\n";
    my $tmp =
`$src/littleprogramrange $matrixfiles $seq $startlist $endlist $namestr $chrom $pos $temp_scores_dir $end $chipfile $strand`;

    #print $tmp;
    return
        $temp_scores_dir . "/" . "chr"
      . "$chrom" . "_" . "$pos" . "_" . "$end" . "_"
      . "$chipfile" . "_"
      . "$strand" . ".txt";
}

#
# Retrieves the appropriate entry from the file created by the procedure above.
# Note that this is desperately slow, due to grep having to iterate through. Do not use.
#
sub all_scores_retrieve {
    my ( $scorefile, $matrix_name ) = @_;
    my $line = `grep $matrix_name $scorefile`;

    #print STDOUT "LINE: ", $line,"\n";
    my @retarr = split ",", $line;
    shift(@retarr);
    return @retarr;
}
#
# Improvement to the above procedure -- reads the file returned by all_scores_two, storing
# the entries into a (matrix, *array) hash.
#
sub all_scores_hash {
    my ($scorefile) = @_;
    my %ret_hash;
    open( scoreFH, $scorefile );
    while ( my $line = <scoreFH> ) {
        next unless $line;
        chomp $line;
        my @array = split ",", $line;
        my $tf = shift(@array);
        $ret_hash{$tf} = \@array;
    }
    close(scoreFH);
    return %ret_hash;
}
#
# Perl computation of the all_scores procedure. Once again, do not use unless absolutely necessary.
# Does not work for indels.
#
sub all_scores {
    my ( $matrix, $seq, $search_start, $search_end ) = @_;
    my $base;
    ## note that limit is EXCLUSIVE, NOT INCLUSIVE, as $position was originally 1-indexed from database
    my $mat = $matrix->matrix;
    my $fscore;
    my $bscore;
    my @seq_array;
    my @sequence = split "", $seq;
    for (
        $base = $search_start ;
        $base <= ( $search_end - $matrix->length + 1 ) ;
        ++$base
      )
    {
        my $pos;
        my $nuc;
        $fscore = 0;
        $bscore = 0;
        for ( $pos = 0 ; $pos < $matrix->length ; ++$pos ) {
            $nuc = $sequence[ $base + $pos - 1 ];
            $nuc =~ tr/ACGT/0123/;
            $fscore += $mat->[$nuc][$pos];
            $bscore += $mat->[ 3 - $nuc ][ $matrix->length - 1 - $pos ];
        }
        $seq_array[ 2 * $base ] = $fscore;
        $seq_array[ 2 * $base + 1 ] = $bscore;
    }

    return @seq_array;
}
#
# Searches for the top score among a specified array, reference/alternative nucleotides, indices, and
# matrices. This does a simple score-swap computation, in O(m) time.
#
sub top_score_search {
    my ( $arrayref, $start, $end, $base, $ref_base, $index, $matrix ) = @_;
    my $matrix_entry;
    my $max_entry = -1024;
    my $j;
    my $matrixref = $matrix->matrix;
    my $comp_matrix_entry;
    my $max_pos;
    my $entryFor;
    my $entryRev;
    my $entry0;
    my $entry1;
    my $best_pos;
    my $best_strand;
    my @good_poss;
    my @good_strand;
    my @good_score;
    $base =~ tr/ACGTacgt/01230123/;
    $ref_base =~ tr/ACGTacgt/01230123/;

    for ( $j = $start ; $j <= ( $end - $matrix->length + 1 ) ; ++$j ) {
        $matrix_entry      = $index - $j;
        $comp_matrix_entry = $matrix->length - 1 - $matrix_entry;
        $entryFor          = @$arrayref[ 2 * $j ];
        $entryRev          = @$arrayref[ 2 * $j + 1 ];
        $entry0 =
          $entryFor -
          $matrixref->[$ref_base][$matrix_entry] +
          $matrixref->[$base][$matrix_entry];
        $entry1 =
          $entryRev -
          $matrixref->[ 3 - $ref_base ][$comp_matrix_entry] +
          $matrixref->[ 3 - $base ][$comp_matrix_entry];
        if ( $entry0 > $entry1 && $entry0 > $max_entry ) {
            $max_entry   = $entry0;
            $best_strand = 1;
            $best_pos    = $j - $start;
        }
        elsif ( $entry1 > $max_entry ) {
            $max_entry   = $entry1;
            $best_strand = -1;
            $best_pos    = $j - $start;
        }
    }

    return ( $max_entry, $best_strand, $best_pos );

}

#
# Computation for the top score for a given position, matrix, start and end. Is faster than
# the slow top_site computation. Note however that top_site returns the coordinates of the
# top site (in a top_site object wrapper) -- this merely returns the score. Simple O(m^2)
# implementation.
# Works for indels. This should be the only time when this is used.
#
sub top_score {
    my ( $matrix, $seq, $search_start, $search_end ) = @_;
    my $base;
    ## note that limit is EXCLUSIVE, NOT INCLUSIVE, as $position was originally 1-indexed from database
    my $limit         = $search_end - $matrix->length;
    my $mat           = $matrix->matrix;
    my $top_score_ret = -1024;
    my $fscore;
    my $bscore;
    ## search_start is counted from a 1-indexed perspective (as $position was), and hence,
    ## the -1 is necessary.
    for ( $base = $search_start ; $base <= $limit + 1 ; ++$base ) {
        my $pos;
        my $nuc;
        $fscore = 0;
        $bscore = 0;
        for ( $pos = 0 ; $pos < $matrix->length ; ++$pos ) {
            $nuc = substr( $seq, $base + $pos - 1, 1 );
            if ( $nuc eq "N" ) {

                # fasta format N equals any base
                $fscore +=
                  ( $mat->[0][$pos] +
                      $mat->[1][$pos] +
                      $mat->[2][$pos] +
                      $mat->[3][$pos] ) / 4;
                $bscore +=
                  ( $mat->[0][ $matrix->length - 1 - $pos ] +
                      $mat->[1][ $matrix->length - 1 - $pos ] +
                      $mat->[2][ $matrix->length - 1 - $pos ] +
                      $mat->[3][ $matrix->length - 1 - $pos ] ) / 4;
            }
            else {
                $nuc =~ tr/ACGT/0123/;
                $fscore += $mat->[$nuc][$pos];
                $bscore += $mat->[ 3 - $nuc ][ $matrix->length - 1 - $pos ];
            }
        }
        if ( $fscore > $bscore && $fscore > $top_score_ret ) {
            $top_score_ret = $fscore;
        }
        elsif ( $bscore > $top_score_ret ) {
            $top_score_ret = $bscore;
        }
    }
    $top_score_ret = sprintf "%.3f", $top_score_ret;

    return $top_score_ret;
}

#
# Reads in position weight matrices from the specified matrix file.
#
sub read_PWMs {
    my ($file) = @_;

    open( FH, $file ) || die "Error opening PWM file $file\n";

    my $matrix_set = TFBS::MatrixSet->new();

    my $name          = '';
    my $matrix_string = '';
    my $line_count    = 0;
    while ( my $line = <FH> ) {
        chomp $line;
        next if !$line;
        if ( $line =~ /^>\s*(\S+)/ ) {
            $name = $1;
        }
        else {
            if ( $line =~ /^\s*[ACGT]\s*\[\s*(.*)\s*\]/ ) {

                # line of the form: A [ # # # ... # ]
                $matrix_string .= "$1\n";
            }
            elsif ( $line =~ /^\s*\d+/ ) {

                # line of the form: # # # ... #
                $matrix_string .= "$line\n";
            }
            else {
                next;
            }
            $line_count++;

            if ( $line_count == 4 ) {
                my $pwm = TFBS::Matrix::PWM->new(
                    -matrixstring => $matrix_string,
                    -name         => $name,
                    -ID           => $name
                );

                $matrix_set->add_Matrix($pwm);

                $line_count    = 0;
                $name          = '';
                $matrix_string = '';
            }
        }
    }
    close(FH);

    return $matrix_set;
}

