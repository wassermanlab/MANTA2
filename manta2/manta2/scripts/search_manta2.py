#!/usr/bin/env python2.7

import os
import sys
import argparse
import re
import pymongo
import bson

from pymongo import MongoClient

script_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.join(os.path.abspath(script_path), os.pardir)
sys.path.append(lib_path)

import fileIO
from dbIO import search_variants


###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    """Search the MANTA2 database for all TFBSs inpacted by SNVs provided
    in the input file. The file may be in VCF, BED, or GFF format.

    Usage: search_impacts.py -i variants_file -o out_file
                -H db_host -d db_name -u db_user -p db_password

    Where:
        -i FILE     - Input file containing a list of variants. The file
                      should contain the following tab delimited columns:
                          chromosome  position  ref_allele  alt_allele

        -o FILE     - Ouput file listing the TFBSs impacted by the given
                      variants.

        -H host     - The MANTA2 database host name.

        -d name     - The name of the MANTA2 database.

        -u host     - The MANTA2 database user name.

        -p password - The MANTA2 database password.
    """

    parser = argparse.ArgumentParser(
        description='Parse the input variant position file to search the database for damaging impacts to TFBSs'
    )

    parser.add_argument(
        '-d', '--db', nargs='?', required=True, help='MANTA2 database name'
    )

    parser.add_argument(
        '-H', '--host', nargs='?', required=True,
            help='MANTA2 database host name'
    )

    parser.add_argument(
        '-u', '--user', nargs='?', required=True,
            help='MANTA2 database user name'
    )

    parser.add_argument(
        '-p', '--password', nargs='?', required=True,
            help='MANTA2 database password'
    )

    parser.add_argument(
        '-i', '--variant_file', nargs='?', required=True, help='Input tab delimited file containing variants to search the database'
    )

    parser.add_argument(
        '-t', '--filetype', nargs='?', help='Input file type - one of: vcf, bed, gff.'
    )

    parser.add_argument(
        '-o', '--out_file', nargs='?', help='If specified, output tab delimited results to this file, otherwise write results to standard output.'
    )

    args = parser.parse_args()

    variant_file = args.variant_file
    out_file = args.out_file
    filetype = args.filetype
    manta2_db = args.db
    manta2_host = args.host
    manta2_user = args.user
    manta2_password = args.password

    if filetype == None:
        ft_result = fileIO.determine_file_format(variant_file)

        filetype = ft_result['filetype']
        ft_errors = ft_result['errors']

        if not filetype or ft_errors:
            sys.stderr.write("Could not determine input variant file type.\n")

            if ft_errors:
                sys.stderr.write("The following errors were encountered trying to determine the input variant file type:\n")
                for e in ft_errors:
                    sys.stderr.write("{}\n".format(e))

            sys.exit("Please fix the input file and/or explicity specify the file type using the -t option\n")

    filetype = filetype.lower()

    var_result = fileIO.read_variants_file(variant_file, filetype)

    uri = "mongodb://{0}:{1}@{2}/{3}".format(
        manta2_user, manta2_password, manta2_host, manta2_db
    )

    try:
        client = MongoClient(uri)
    except:
        sys.exit("Could not connect to the MANTA2 database. Please make sure you have specified the database host name, db name, user and password correctly")

    db = client[manta2_db]

    variants = var_result['variants']
    var_errors = var_result['errors']

    search_result = search_variants(db, variants)

    snv_impacts = search_result['snv_impacts']
    snv_errors = search_result['errors']

    client.close()

    fileIO.write_snv_impacts(out_file, snv_impacts)

    if var_errors or snv_errors:
        print("\nProblems were detected while searching the MANTA2 database:")

    if var_errors:
        for err in var_errors:
            print("\n{}".format(err))

    if snv_errors:
        for err in snv_errors:
            print("\n{}".format(err))

    if out_file and out_file != '-':
        print("\nSNV impact results written to file {}.\n".format(out_file))
