#!/usr/bin/env python3
import sys
import os
import signal





def nextflow_runner(alignment,tree,ambig,iteration):
    # Need to repeatedly run the netflow pipe until a set proportion of species are covered.
    os.system(f"nextflow run -profile conda ~/repositories/primer_piper/main.nf --alignment {alignment} --tree {tree} --ambig {int(ambig)} --iter {iteration}")




if __name__ == "__main__":
    # Parse the inputs, expect `input_alignment primer1_ref_amplicons primer2_ref_amplicons...`

    try:
        input_alignment = sys.argv[1]
        input_tree = sys.argv[2]
        input_target = sys.argv[3]
        input_ambig = sys.argv[4]
        
        iteration = 1
        # Need to run it once to see what coverage we are dealing with
        nextflow_runner(input_alignment,input_tree,input_ambig,iteration)

        refs_covered_list = []
        with open(f"{os.getcwd()}/output/Coverage_poportion_{iteration}.csv", "r") as iteration_file:
            for line in iteration_file:
                refs_covered_list.append(float(line.split(",")[2]))
                if (line.split(",")[0] == "Iteration_1"):
                    alignment_size = float(line.split(",")[3])
        cumulative_covered = sum(refs_covered_list)
        proportion = (cumulative_covered/alignment_size) * 100

        while True:
            with open(f"{os.getcwd()}/output/Coverage_poportion_{iteration}.csv", "r") as iteration_file:
                    for line in iteration_file:
                        refs_covered_list.append(float(line.split(",")[2]))
                        if (line.split(",")[0] == "Iteration_1"):
                            alignment_size = float(line.split(",")[3])
            cumulative_covered = sum(refs_covered_list)
            proportion = (cumulative_covered/alignment_size) * 100
            print(f"\n\n\n\nAt Iteration {iteration} the proportion of references covered is: {proportion}\n\n\n\n")
            if proportion < float(input_target):
                
                updated_alignment = ("%s/output/updated_alignment_%s.aln" % (os.getcwd(),iteration))
                iteration += 1
                nextflow_runner(updated_alignment,input_tree,input_ambig,iteration)
            else:
                print(f"\n\n\n\nAll Done!\n\n")
                sys.exit(0)
    except KeyboardInterrupt:
        print('Muerte!')