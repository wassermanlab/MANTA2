import os
import sys
import re

REF_ALLELE_REGEX    = re.compile('[ACGTNactgn]')
ALT_ALLELE_REGEX    = re.compile('[ACGTactg]')
VCF_COMMENT_REGEX   = re.compile('#.*fileformat=VCF.*')

def determine_file_format(filename):
    """ Try to determine the input variant file format based on the columns
    within and possibly the extension.
    """

    filetype = None
    errors = []

    # First try to use the extension to determine file type
    fname, ext = os.path.splitext(filename)

    if ext:
        ext = ext.lstrip('.')
        ext = ext.lower()

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith("#"):
                if VCF_COMMENT_REGEX.match(line):
                    filetype = 'vcf'
                    break
                else:
                    # skip other comment lines
                    continue

            line = line.rstrip()

            #
            # NOTE: VCF and GFF formats explicitly state that the columns
            # should be tab delimited, however BED format allows columns to
            # be delimited on any whitespace. First try to split on tabs.
            #
            if ext == 'vcf' or ext == 'gff':
                cols = line.split('\t')

                if len(cols) < 4:
                    # Less than 4 columns when splitting on tabs.
                    errors.append("Provided variant file has less than 4 columns")
                    break
            else:
                # Either a BED file or unknown by extension. First try to
                # split by tabs.
                cols = line.split('\t')
                if len(cols) < 4:
                    cols = line.split('\s+')

                    if len(cols) < 4:
                        errors.append("Provided variant file has less than 4 columns")
                        break

            # File has at least 4 columns
            if is_int(cols[1]):
                # Both VCF and BED files have integer start position in
                # the 2nd column 
                if is_int(cols[2]):
                    # BED files have integer end coordinate in the 2nd
                    # column but it's possible it could still be a VCF
                    # where the ID field contains an integer ID.
                    # If it's a BED the end coordinate should be 1 greater
                    # than the start coordinate
                    if (int(cols[2]) - int(cols[1]) == 1 and len(cols) == 4
                        and len(cols[3]) == 1):
                        # Should be a BED file
                        filetype = "bed"
                    else:
                        # Should be a VCF file. Check it has at least 5
                        # columns and the 4th and 5th columns look like
                        # single nucleotides
                        if (len(cols) >= 5 and len(cols[3]) == 1
                            and len(cols[4]) == 1):
                            filetype = "vcf"
                else:
                    # Should be a VCF file.
                    # Check we have at least 5 columns
                    if (len(cols) >= 5 and len(cols[3]) == 1
                        and len(cols[4]) == 1):
                        filetype = "vcf"
            else:
                # Should be a GFF file. Check we have at least 9 columns
                if len(cols) >= 9:
                    filetype = "gff"
                else:
                    errors.append("Provided variant file appears to be a GFF file but has less than 9 columns")

            break   # stop reading lines from the file

        return {'filetype' : filetype, 'errors' : errors}



def read_variants_file(filename, filetype):
    """Read the input variants file and return the variants as a list.
    """

    #sys.stderr.write("Reading variants file {0} of type {1}\n".format(filename, filetype))

    variants = []
    errors = []

    with open(filename, 'r') as f:
        if filetype == 'vcf':
            for line in f:
                if line.startswith("#"):
                    continue

                line = line.rstrip()

                cols = line.split('\t')

                if len(cols) < 5:
                    errors.append("For VCF formatted input, <em>please make sure</em> that at least the first <b>5</b> columns are provided. <em>PLEASE NOTE</em> that the VCF format specification requires that the columns be <b>tab-separated</b> (this error often occurs if the VCF input columns were space-separated). <br><br>Also, <em>PLEASE NOTE</em> that if you do not have a value for the <b>ID</b> field, please use a dot (\".\") instead of leaving it blank.")
                    return {'variants' : variants, 'errors' : errors}

                chrom       = cols[0]
                position    = cols[1]
                id          = cols[2]
                ref_allele  = cols[3]
                alt_allele  = cols[4]

                if not check_alleles(ref_allele, alt_allele):
                    continue

                chrom = chrom.lstrip('chr')

                var = {
                    'id'    : id,
                    'chrom' : chrom,
                    'position' : int(position),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        elif filetype == 'gff':
            line_num = 0
            for line in f:
                line_num += 1

                if line.startswith("#"):
                    continue

                line = line.rstrip()

                cols = line.split('\t')

                if len(cols) < 9:
                    errors.append("For GFF formatted input, <em>please make sure</em> that <b>9</b> columns are provided. <em>PLEASE NOTE</em> that the GFF format specification requires that the columns be <b>tab-separated</b> (this error often occurs if the GFF input columns were space-separated). <br><br>Also, <em>PLEASE NOTE</em> that if you do not have a value for one of the \"optional\" fields, e.g. the <samp>source</samp> field, please use a dot (\".\") instead of leaving it blank.")
                    return {'variants' : variants, 'errors' : errors}

                chrom, data_source, feature_type, start, end, score, strand, frame, attribute_list = cols[:9]

                if start != end:
                    # In this case, assume we have an indel rather than
                    # an SNV and ignore it
                    continue
                
                ref_allele = None
                alt_allele = None

                if attribute_list:
                    attributes = attribute_list.split('; ')
                    for attr in attributes:
                        attr_name, attr_val = attr.split('=')
                        if attr_name == 'ref_allele' or attr_name == 'reference_allele':
                            ref_allele = attr_val
                        if attr_name == 'alt_allele' or attr_name == 'alt_allele':
                            alt_allele = attr_val
                else:
                    errors.append("GFF file has blank attributes field at line {}".format(line_num))
                    return {'variants' : variants, 'errors' : errors}

                if not ref_allele or not alt_allele:
                    errors.append("GFF attributes field contains no ref_allele or alt_allele information at line {}".format(line_num))

                    return {'variants' : variants, 'errors' : errors}

                if not check_alleles(ref_allele, alt_allele):
                    continue

                chrom = chrom.lstrip('chr')

                #
                # XXX
                # If strand is given and negative, should we reverse
                # complement the alleles?
                # XXX
                #

                #
                # Assume the name field holds the alt. allele. We don't
                # know the ref. allele, so set it to 'N'
                #
                var = {
                    'chrom' : chrom,
                    'id'    : '.',
                    'position' : int(start),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        elif filetype == 'bed':
            for line in f:
                if line.startswith("#"):
                    continue

                line = line.rstrip()

                #sys.stderr.write("BED line: {0}\n".format(line))

                #
                # NOTE: The BED format does not explicitly state that
                # the columns should be tab delimited and in fact when
                # loading BED files into the UCSC genome browser, the
                # browser splits on any whitespace. But since it may be
                # commonly assumed that BED files should be tab-delimited
                # try to first split on tabs and if that 'fails', attempt
                # to split on whitespace.
                #
                cols = line.split('\t')

                if len(cols) < 4:
                    cols = line.split('\s+')

                    if len(cols) < 4:
                        errors.append("For BED formatted input, <em>please make sure</em> that at least the first <b>4</b> columns are provided. The fields in the BED lines may be tab or space separated.")

                        return {'variants' : variants, 'errors' : errors}

                chrom = cols[0]
                start = int(cols[1])
                end   = int(cols[2])
                # NOTE: using the name field as alt. allele
                alt_allele = cols[3]

                # Don't have ref. allele info so set to '.'
                # This is used later to determine if the input file was BED
                # and therefore the ref. allele is actually unknown.
                ref_allele = '.'

                if start < end - 1:
                    # In this case, assume we have an indel rather than
                    # an SNV and ignore it
                    continue

                if len(cols) >= 6:
                    score  = cols[4]
                    strand = cols[5]

                #
                # XXX
                # If strand is given and negative, should we reverse
                # complement the alleles?
                # XXX
                #

                if not check_alleles(ref_allele, alt_allele):
                    continue

                chrom = chrom.lstrip('chr')

                #
                # Assume the name field holds the alt. allele. We don't
                # know the ref. allele, so set it to 'N'
                #
                var = {
                    'chrom' : chrom,
                    'id'    : '.',
                    'position' : int(end),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        elif filetype == 'simple':
            #
            # My own test file type. Not used...
            #
            for line in f:
                if line.startswith("#"):
                    continue

                line = line.rstrip()

                # NOTE: assuming 1-based coordinates for position
                cols = line.split('\t')

                if len(cols) < 4:
                    errors.append("The provided SNV file appears to contain less than 4 columns. Please make sure file is tab-delimited.")

                    return {'variants' : variants, 'errors' : errors}

                chrom      = cols[0]
                position   = cols[1]
                ref_allele = cols[2]
                alt_allele = cols[3]

                chrom = chrom.lstrip('chr')

                var = {
                    'chrom' : chrom,
                    'id'    : '.',
                    'position' : int(position),
                    'ref_allele' : ref_allele,
                    'alt_allele' : alt_allele
                }

                variants.append(var)

        else:
            errors.append("Unknown SNV file type {0}".format(filetype))
            return {'variants' : variants, 'errors' : errors}

    return {'variants' : variants, 'errors' : errors}


def write_snv_impacts(filename, snv_impacts):
    """Write the TFBS impacts for each variant to the given output file.
    """

    if filename:
        fh = open(filename, 'w')
    else:
        fh = sys.stdout

    #
    # Note, the values in the snv_impacts structure are now already
    # pre-formatted for consistency in both the web and flat file display.
    #
    for si in snv_impacts:
        fh.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\n".format(si['chrom'], si['position'], si['ref_allele'], si['alt_allele'], si['snv_id'], si['jaspar_tf_name'], si['matrix_id'], si['start1'], si['end1'], si['strand1'], si['abs_score1'], si['rel_score1'], si['start2'], si['end2'], si['strand2'], si['abs_score2'], si['rel_score2'], si['impact']))

    fh.close()

    return


def check_alleles(ref_allele, alt_allele):
    """Check alleles to make sure this looks like an SNV
    """

    if len(ref_allele) != 1:
        return False
    if len(alt_allele) != 1:
        return False
    if ref_allele.upper() == alt_allele.upper():
        return False
    if not (REF_ALLELE_REGEX.match(ref_allele) or ref_allele == '.'):
        return False
    if not ALT_ALLELE_REGEX.match(alt_allele):
        return False

    return True


def is_int(val):
    """Check if value is an interger.
    """

    try:
        int(val)
    except:
        return False
    else:
        return True
