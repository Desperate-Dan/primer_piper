#!/user/bin/env/ python3
import sys
from Bio import SeqIO

def primer_splitter(input_primers):
    for seq_rec in SeqIO.parse(input_primers, "fasta"):
        with open("potential_primers_" + seq_rec.id.split("_")[1] + ".fasta", "a") as file:
            SeqIO.write(seq_rec, file, "fasta")

if __name__ == "__main__":
    input_primers = sys.argv[1]
    primer_splitter(input_primers)