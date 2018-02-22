#!/raid2/local/python2.7.3/bin/python2.7
#*-* coding: utf-8 *-*

import argparse
import os
import glob

LOAD_SCRIPT = '/home/dave/projects/manta2/scripts/load_manta2_db.py'

def create_qsub_files(data_dir, qsub_dir, db_user, db_pass):
    for file in os.listdir(data_dir):
        if file.endswith(".bed"):
            create_bed_qsub_file(data_dir, qsub_dir, file, db_user, db_pass)


def create_bed_qsub_file(data_dir, qsub_dir, bed_file, db_user, db_pass):

    bed_base_name = os.path.basename(bed_file)  # without path
    bed_root_file = os.path.splitext(bed_base_name)[0]  # without extension

    jaspar_matrix_id, tf_name = bed_root_file.split('_')

    log_file = "manta2_load_{0}.log".format(jaspar_matrix_id)
    log_path = os.path.join(qsub_dir, 'log', log_file)  

    qsub_sh_file = "manta2_load_{0}.sh".format(jaspar_matrix_id)
    qsub_out_file = "manta2_load_{0}.out".format(jaspar_matrix_id)
    qsub_err_file = "manta2_load_{0}.err".format(jaspar_matrix_id)

    qsub_sh_path = os.path.join(qsub_dir, 'sh', qsub_sh_file)
    qsub_out_path = os.path.join(qsub_dir, 'out',  qsub_out_file)
    qsub_err_path = os.path.join(qsub_dir, 'err', qsub_err_file)
    jobname = "m2{0}".format(jaspar_matrix_id)

    db_load_cmd = "{} -d {} -f {} -u {} -p {}".format(
                    LOAD_SCRIPT, data_dir, bed_base_name, db_user, db_pass)

    fh = open(qsub_sh_path, 'w')

    fh.write("#!/bin/bash\n")
    fh.write("#$ -N {0}\n".format(jobname))
    fh.write("#$ -r n\n")
    fh.write("#$ -V\n")
    fh.write("#$ -o {0}\n".format(qsub_out_path))
    fh.write("#$ -e {0}\n".format(qsub_err_path))
    fh.write("#$ -m a\n")
    fh.write("#$ -M dave@cmmt.ubc.ca\n")
    fh.write("#$ -l mem_free=4G,h_vmem=4G,h_rt=4:0:0\n")
    fh.write("\necho $HOSTNAME\n".format(db_load_cmd))
    fh.write("\npython2.7 {}\n".format(db_load_cmd))

    fh.close()


###############################################################################
#                               MAIN
###############################################################################
if __name__ == '__main__':
    '''
    Create qsub files for each JASPAR TF ID for loading the CRV database.

    Usage: create_load_db_qsub_files.py -dd data_dir -qd qsub_dir

    Expects the data_dir to contain a set of BED files containing
    the list of experiments that the individual qsub job will process. All
    experiments in a CSV file should be related to one specific JASPAR TF ID.
    '''

    parser = argparse.ArgumentParser(
        description='Create qsub files for each of the individual JASPAR TF data sets files.'
    )

    parser.add_argument(
        '-dd', '--data-dir', nargs='?', required=True, help='Top level directory containing CSV, TFBS BED files and subdirectories with the SNV files.'
    )

    parser.add_argument(
        '-qd', '--qsub-dir', nargs='?', required=True, help='Output directory to which the qsub files are written' 
    )

    parser.add_argument(
        '-u', '--db-user', nargs='?', required=True, help='MANTA2 DB user' 
    )

    parser.add_argument(
        '-p', '--db-pass', nargs='?', required=True, help='MANTA2 DB password' 
    )

    args = parser.parse_args()

    data_dir = args.data_dir
    qsub_dir = args.qsub_dir
    db_user = args.db_user
    db_pass = args.db_pass

    create_qsub_files(data_dir, qsub_dir, db_user, db_pass)
