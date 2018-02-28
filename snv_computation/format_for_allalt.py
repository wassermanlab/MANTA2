#!/usr/bin/python
#*-* coding: utf-8 *-*

###############################################################################
#
###############################################################################

import sys
import getopt
from Bio import SeqIO


def get_sequences(fasta_file):
    with open(fasta_file) as stream:
        seqs = {}
        for record in SeqIO.parse(stream, "fasta"):
            seqs[record.id] = str(record.seq).upper()
        return seqs


def construct_output(sequences, bedfile):
    with open(bed_file) as stream:
        for line in stream:
            spl = line.rstrip().split()
            abs_score, rel_score = spl[4].split(';')
            seq = sequences[spl[3]]
            strand = "1"
            if spl[5] == "-":
                strand = "-1"
            print("{0}\t{1}\t{2:d}\t{3}\t{4}\t{5}\t{6}\t".format(spl[0],
                    strand, eval(spl[1]) + 1, spl[2], spl[3],
                    eval(rel_score)/100., abs_score)),
            print(".\tspecies\t{0}".format(seq))


###############################################################################
#                               MAIN
###############################################################################
if __name__ == "__main__":
    usage = '''
    %s -b <bed file> -f <fasta file>
    '''%(sys.argv[0])

    try:
        opts, args = getopt.getopt(sys.argv[1:], "b:f:h")
    except getopt.GetoptError:
        sys.exit(str(getopt.GetoptError) + usage)

    bed_file = None
    fasta_file = None
    for o, a in opts:
        if o == "-b":
            bed_file = a
        elif o == "-f":
            fasta_file = a
        else:
            sys.exit(usage)
    if not(bed_file and fasta_file):
        sys.exit(usage)

    sequences = get_sequences(fasta_file)
    construct_output(sequences, bed_file)
