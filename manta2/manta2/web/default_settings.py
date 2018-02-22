# Define various constants needed by the MANTA2 web application (including
# standalone scripts). This file is under version control and should not
# contain anything that is not supposed to be exposed to the public.

import os
import re

APPLICATION_ROOT = "/manta2"

# Allowed extensions for input variant files.
ALLOWED_EXTENSIONS  = ['vcf', 'bed', 'gff', 'txt']

# If true, prompt user to choose variant file format explicitly 
ASK_FOR_FILETYPE = False

#
# The following setttings depend on the MANTA_ROOT_DIR being set. This is set
# in another config file which is private (not stored in the application tree
# and not under version control) and which is loaded via the enviroment
# variable MANTA2_PRIVATE_SETTINGS which is set via:
#   export MANTA2_PRIVATE_SETTINGS=/path/to/private_settings.cfg
# Then in the web application:
#   app.config.from_envvar('MANTA2_PRIVATE_SETTINGS')
#
LOG_FILE        = 'manta2.log'
ERROR_LOG       = os.path.join('logs', LOG_FILE)

UPLOAD_REL_DIR   = 'upload'

DOWNLOAD_REL_DIR = 'download'
RESULTS_REL_DIR  = os.path.join(DOWNLOAD_REL_DIR, 'results')

# The results filename is saved in a directory which consists of the
RESULTS_FILENAME = 'manta2_results.txt'

# Clean up results files older than this number of days
REMOVE_RESULTFILES_OLDER_THAN   = 7 # days

# Flask-PyMongo connects using PREFIX_HOST, PREFIX_DBNAME etc.
MANTA2_HOST     = 'manta.cmmt.ubc.ca'
MANTA2_DBNAME   = 'manta2'
MANTA2_USERNAME = 'manta_r'
MANTA2_PASSWORD = 'mantapw'

# Used to check that provided variants are in the accepted nucleotide alphabet
REF_ALLELE_REGEX    = re.compile('[ACGTNactgn]')
ALT_ALLELE_REGEX    = re.compile('[ACGTactg]')
