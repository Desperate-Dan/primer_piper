#!/user/bin/env/ python3

import re
import seaborn as sns
import matplotlib as mpl
import sys
from matplotlib import pyplot as plt
from collections import defaultdict
import gzip

input_file = sys.argv[1]

if len(sys.argv) > 2:
    out_name = sys.argv[2]
else:
    out_name = "primercounter"

if re.search(".gz", input_file):
    with gzip.open (input_file, "rt") as file:
        primer_dict = defaultdict(set)
        individual_primers_dict = {}
        individual_primers_amplicons_dict = defaultdict(set)
        ref_amplicons_dict = defaultdict(set)
        pairs_dict = defaultdict(set)
        for line in file:
            if re.search(">Amp_", line):
    #This bit of code counts the number of times each primer is part of an amplicon
    #Probably not a meaningful metric as will inflate count of static primers that match with a redundant primer binding site. 
                if line.split()[1] in individual_primers_dict.keys():
                    individual_primers_dict[line.split()[1]] += 1
                else:
                    individual_primers_dict[line.split()[1]] = 1
                if line.split()[3] in individual_primers_dict.keys():
                    individual_primers_dict[line.split()[3]] += 1
                else:
                    individual_primers_dict[line.split()[3]] = 1
    #This bit of code counts specific primer involvement in amplicons 
                individual_primers_amplicons_dict[line.split()[1]].add(line.split()[5])
                individual_primers_amplicons_dict[line.split()[3]].add(line.split()[5])
    #This bit of code is looking at what ref the predicted amplicon is hitting and which region in that ref
                ref_amplicons_dict[line.split()[5].split(":")[0]].add(line.split()[5].split(":")[1])
                
    #This bit of code below gets the primer, removes the suffix of which number it is, and adds the amplicon position to a set
    #This in theory gives a number to how many amplicons a primer contributes to, while removing multiple hits to the same site
                primer_dict[line.split()[1].split(".")[0]].add(line.split()[5])
                primer_dict[line.split()[3].split(".")[0]].add(line.split()[5])
    #This bit of code compares the f and r primers of an amplicon
                pair = line.split()[1].split("_")[0] + "_F_and_" + line.split()[3].split("_")[0] + "_R"
                pairs_dict[pair].add(line.split()[5])
        
else:
    with open (input_file, "r") as file:
        primer_dict = defaultdict(set)
        individual_primers_dict = {}
        individual_primers_amplicons_dict = defaultdict(set)
        ref_amplicons_dict = defaultdict(set)
        pairs_dict = defaultdict(set)
        for line in file:
            if re.search(">Amp_", line):
    #This bit of code counts the number of times each primer is part of an amplicon
    #Probably not a meaningful metric as will inflate count of static primers that match with a redundant primer binding site. 
                if line.split()[1] in individual_primers_dict.keys():
                    individual_primers_dict[line.split()[1]] += 1
                else:
                    individual_primers_dict[line.split()[1]] = 1
                if line.split()[3] in individual_primers_dict.keys():
                    individual_primers_dict[line.split()[3]] += 1
                else:
                    individual_primers_dict[line.split()[3]] = 1
    #This bit of code counts specific primer involvement in amplicons 
                individual_primers_amplicons_dict[line.split()[1]].add(line.split()[5])
                individual_primers_amplicons_dict[line.split()[3]].add(line.split()[5])
    #This bit of code is looking at what ref the predicted amplicon is hitting and which region in that ref
                ref_amplicons_dict[line.split()[5].split(":")[0]].add(line.split()[5].split(":")[1])

    #This bit of code below gets the primer, removes the suffix of which number it is, and adds the amplicon position to a set
    #This in theory gives a number to how many amplicons a primer contributes to, while removing multiple hits to the same site
                primer_dict[line.split()[1].split(".")[0]].add(line.split()[5])
                primer_dict[line.split()[3].split(".")[0]].add(line.split()[5])
    #This bit of code compares the f and r primers of an amplicon
                pair = line.split()[1].split("_")[0] + "_F_and_" + line.split()[3].split("_")[0] + "_R"
                pairs_dict[pair].add(line.split()[5])

primer_count_dict = {}
for k, v in primer_dict.items():
    primer_count_dict[k] = len(v)
with open (out_name + "_primer_counts.csv", "w") as file:
    file.write("primer,count\n")
    for k, v in primer_count_dict.items():
        file.write(k + "," + str(v) + "\n")


individual_primer_count_dict = {}
for k, v in individual_primers_amplicons_dict.items():
    individual_primer_count_dict[k] = len(v)
with open (out_name + "_individual_primer_amplicon_counts.csv", "w") as file:
    file.write("primer,count\n")
    for k, v in individual_primer_count_dict.items():
        file.write(k + "," + str(v) + "\n")


pairs_count_dict = {}
for k, v in pairs_dict.items():
    pairs_count_dict[k] = len(v)
with open (out_name + "_pairs_count.csv", "w") as file:
    file.write("primer_pair,count\n")   
    for k, v in pairs_count_dict.items():
        file.write(k + "," + str(v) + "\n")

with open (out_name + "_ref_amplicons.csv", "w") as file:
    file.write("reference,amplicon_pos,length\n")
    for k, v in ref_amplicons_dict.items():
        for pos in v:
            length = str(int(pos.split("-")[1]) - int(pos.split("-")[0])) 
            file.write(k + "," + str(pos) + "," + length + "\n")


sns.barplot(primer_count_dict)
plt.xticks(rotation=45, ha='right')
plt.savefig(out_name + "_figure_primer.svg", bbox_inches="tight")
plt.close()

sns.barplot(pairs_count_dict)
plt.xticks(rotation=45, ha='right')
plt.savefig(out_name + "_pairs_primer.svg", bbox_inches="tight")
plt.close()
