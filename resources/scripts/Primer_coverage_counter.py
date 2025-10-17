#!/usr/bin/env python3

# Then plan here is to take in the refs that each newly generated primer covers, count the number of refs the best primer covers, see if that matches
# a target value (ie 75% of family members; may have to define that separately to feed the nextflow "recurse" method), write a file out with those
# primers in, "pop" out the references that have already been covered from the alignment file, then allow the primer generation to restart.

import sys
from Bio import SeqIO
from collections import defaultdict

def alignment_parser(input_alignment):
    records = list(SeqIO.parse(input_alignment, "fasta"))
    alignment_size = len(records)
    #print(alignment_size)
    return alignment_size


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

def fun_with_sets(ref_hits_dict,target):
    size_dict = defaultdict(list)
    for k,v in ref_hits_dict.items():
        size_dict[len(v)].append(k)
    max_primer_list = size_dict[max(size_dict)]
    for i in max_primer_list:
        # This will currently just be the last iters primers - is that a concern?
        refs_covered_count = len(ref_hits_dict[i])
        ref_set = ref_hits_dict[i]

    return max_primer_list, refs_covered_count, ref_set

def primer_output(max_primer_list):
    with open("output.csv","a") as output:
        for primer_file in max_primer_list:
            with open(f"{primer_file}.csv", "r") as file:
                for line in file:
                    line = line.rstrip("\n")
                    output.write(f"{line},{primer_file}\n")

def alignment_popper(input_alignment, ref_set):
    records_list = list(SeqIO.parse(input_alignment, "fasta"))
    #print(records_list)
    with open("example.fasta", "w") as output_handle:
        for record in records_list:
            if record.id in ref_set:
                continue
            else:
                SeqIO.write(record,output_handle,"fasta")

if __name__ == "__main__":
    # Parse the inputs, expect `input_alignment target_proportion previous_coverage primer1 primer2...`
    input_alignment = sys.argv[1]
    target = sys.argv[2]
    previous_coverage = sys.argv[3]
    input_ref_counts = sys.argv[4:]
    
    alignment_size = alignment_parser(input_alignment)
    ref_hits_dict = input_combiner(input_ref_counts)
    max_primer_list, refs_covered_count, ref_set = fun_with_sets(ref_hits_dict,target)
    primer_output(max_primer_list)
    alignment_popper(input_alignment, ref_set)