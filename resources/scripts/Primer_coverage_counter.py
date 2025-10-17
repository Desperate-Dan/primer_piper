#!/usr/bin/env python3

import sys
from Bio import SeqIO
from collections import defaultdict

def input_combiner(ref_counts):
    ref_hits_dict = defaultdict(set)
    for in_file in ref_counts:
        with open(in_file,"r") as file:
            for line in file:
                if line.split(",")[0] == "reference":
                    continue
                ref_hits_dict[in_file.split(".csv")[0]].add(line.split(",")[0])
    #print(ref_hits_dict)
    return ref_hits_dict

def fun_with_sets(ref_hits_dict):
    print(max(ref_hits_dict, key=ref_hits_dict.values))
    
if __name__ == "__main__":
    input_alignment = sys.argv[1]
    input_ref_counts = sys.argv[2:]
    ref_hits_dict = input_combiner(input_ref_counts)
    fun_with_sets(ref_hits_dict)