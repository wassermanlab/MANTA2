#!/usr/bin/env python2.7
import os, sys, re

def parse_file(file_name):
    """
    This function parses any file and yields lines one by one.
    
    @input:
    file_name {string}
    @return:
    line {string}

    """
 
    if os.path.exists(file_name):
        # Initialize #
        f = None
        # Open file handle #
        try: f = open(file_name, "rt")
        except: raise ValueError("Could not open file %s" % file_name)
        # For each line... #
        for line in f:
            yield line.strip()
        f.close()
    else:
        raise ValueError("File %s does not exist!" % file_name)

# Initialize #
asb = []
coordinates = []
liftOver = []
script_path = os.path.dirname(os.path.realpath(__file__))
asb_file = os.path.join(script_path, "gkw691_Supp", "ASB_GM12878_HeLa_1based.txt")
coordinates_file = os.path.join(script_path, "coordinates.hg19.bed")
coordinates_liftOver_file = os.path.join(script_path, "coordinates.hg38.liftOver.bed")

# Read ASB data #
for line in parse_file(asb_file):
    if line.startswith("chr\tstart\tref"): continue
    asb.append(line.split("\t"))
# Read coordinates #
for line in parse_file(coordinates_file):
    coordinates.append(line)
# Read liftOver coordinates #
for line in parse_file(coordinates_liftOver_file):
    liftOver.append(line.split("\t"))
# For each ASB event... #
for i in range(len(asb)):
    if coordinates[i].startswith("#"): continue
    a = liftOver.pop(0)
    print("%s\t%s\t%s\t%s" % (a[0], a[1], a[2], asb[i][3]))
