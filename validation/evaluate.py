#!/usr/bin/env python2.7
import os, sys, re
from matplotlib.offsetbox import AnchoredText
#from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import pandas
from scipy.stats import pearsonr

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
manta2 = {}
data = []
script_path = os.path.dirname(os.path.realpath(__file__))
asb_file = os.path.join(script_path, "gkw691_Supp", "ASB_GM12878_HeLa_1based.txt")
coordinates_file = os.path.join(script_path, "coordinates.hg19.bed")
coordinates_liftOver_file = os.path.join(script_path, "coordinates.hg38.liftOver.bed")
manta2_file = os.path.join(script_path, "MANTA2.out")
data_frame_file = os.path.join(script_path, "data_frame.csv")
figure_file = os.path.join(script_path, "evaluation")

# If data frame does not exist... #
if not os.path.exists(data_frame_file):
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
    # Read MANTA2 output #
    for line in parse_file(manta2_file):
        line = line.split("\t")
        manta2.setdefault(line[0], [])
        manta2[line[0]].append(line)
    # For each ASB event... #
    for i in range(len(asb)):
        if coordinates[i].startswith("#"): continue
        # Initialize #
        a = liftOver.pop(0)
        ai = str(float(asb[i][7])/(int(asb[i][6]) + int(asb[i][7]))) # allelic imbalance on the alternate allele
        # Skip if wrong key #
        if a[0][3:] not in manta2: continue
        # For each MANTA2 result... #
        for j in manta2[a[0][3:]]:
            if a[2] == j[1] and asb[i][5].upper() in j[5] and asb[i][2] == j[2] and asb[i][3] == j[3]:
                data.append([j[0], j[1], j[2], j[3], asb[i][4], asb[i][5], asb[i][-1], j[6], ai, j[-1]])
    # Make data frame #
    data_frame = pandas.DataFrame(data)
    data_frame.columns = ["chr", "pos", "ref", "alt", "cell", "tf", "asb", "jaspar_id", "allelic_imbalance", "impact_score"]
    # Save data frame to CSV #
    data_frame.to_csv(data_frame_file)

# Make figure #
data_frame = pandas.read_csv(data_frame_file)
fig = plt.figure(tight_layout=True)
ax = fig.add_subplot(111)
# Plot scatter #
#df = data_frame[(data_frame["asb"]=="nonASB")] # nonASB events
#ax.scatter(df["allelic_imbalance"].tolist(), df["impact_score"].tolist(), color="#4477AA", alpha=0.5)
#df = data_frame[(data_frame["asb"]=="ASB")] # ASB events
#ax.scatter(df["allelic_imbalance"].tolist(), df["impact_score"].tolist(), color="#CC6677", alpha=0.5)
ax.scatter(data_frame["allelic_imbalance"].tolist(), data_frame["impact_score"].tolist(), color="#4477AA", alpha=0.5)
# Get correlation #
corr = pearsonr(data_frame["allelic_imbalance"].tolist(), data_frame["impact_score"].tolist())
# Add legend #
#handles = [Patch(color="#CC6677", label="ASB"),
#           Patch(color="#4477AA", label="nonASB"),
#           Patch(color="none", label="R = %.3f" % corr[0]),
#           Patch(color="none", label='$\it{p}$' + " = %.1e" % corr[1])]
#ax.legend(handles=handles, frameon=False, loc="best")
anchored_text = "R = %.3f\n" % corr[0] + '$\it{p}$' + " = %.1e" % corr[1]
ax.add_artist(AnchoredText(anchored_text, loc=2, frameon=False))
# Set x/y labels #
ax.set_xlabel("allelic imbalance (ChIP-seq)")
ax.set_ylabel("impact score (MANTA2)")
# Set aspect ratio #
x0, x1 = ax.get_xlim()
y0, y1 = ax.get_ylim()
ax.set_aspect((x1 - x0) / (y1 - y0))
# Save figure #
fig.savefig("%s.svg" % figure_file, dpi=300, format='svg', transparent=True)
fig.savefig("%s.png" % figure_file, dpi=300, format='png', transparent=True)
plt.close("all") # closes all plots
