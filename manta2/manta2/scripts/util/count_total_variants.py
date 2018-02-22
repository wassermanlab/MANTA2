#!/usr/bin/env python2.7

import os
import Bio

import Bio.motifs.jaspar as jaspar
from Bio.motifs.jaspar.db import JASPAR5

JASPAR_DB_HOST = 'vm5.cmmt.ubc.ca'
JASPAR_DB_NAME = 'JASPAR_2018'
JASPAR_DB_USER = 'jaspar_r'
JASPAR_DB_PASS = ''

DFLT_COLLECTION = 'CORE'
jdb = JASPAR5(
        host=JASPAR_DB_HOST,
        name=JASPAR_DB_NAME,
        user=JASPAR_DB_USER,
        password=JASPAR_DB_PASS
    )

BED_FILE_DIR = '/raid9/dave/MANTA2/data/20171101_variant_score_uppercase'

os.chdir(BED_FILE_DIR)

total_variants = 0
total_bed_files = 0
total_binding_sites = 0

for file in os.listdir(BED_FILE_DIR):
    if file.endswith(".bed"):
        total_bed_files += 1

        line_count = 0
        with open(file) as f:
            for i, l in enumerate(f):
                pass

        line_count = i + 1

        total_binding_sites += line_count
            
        bed_root_file = os.path.splitext(os.path.basename(file))[0]
        (matrix_id, tf_name) = bed_root_file.split('_')

        motif = jdb.fetch_motif_by_id(matrix_id)

        total_variants += len(motif) * 3 * line_count

print "Totals: motifs = {}\tbinding_sites = {}\tvariants = {}".format(
        total_bed_files, total_binding_sites, total_variants)




