process pickPrimers {
    conda "${HOME}/miniconda3/envs/varvamp"
    publishDir "output/varvamp/"

    debug true

    input:
    val ambig_val
    path in_alignment
    
    output:
    path("${in_alignment.simpleName}")
    path("${in_alignment.simpleName}/primers.fasta"), emit: new_primers

    script:
    """
    varvamp single -a ${ambig_val} ${in_alignment} ${in_alignment.simpleName}
    """
}

process splitPrimers {
    //If there are multiple prospective primer pairs, we need to split them out before checking family coverage
    conda "${HOME}/miniconda3/envs/code_hole"
    publishDir "output/varvamp"

    debug true

    input:
    path new_primers

    output:
    path("potential_primers_*.fasta"), emit: new_split

    script:
    """
    python3 ${projectDir}/resources/scripts/primer_splitter.py ${new_primers}
    """
}

process checkFamilycov {
    //Running MFEprimer to check if the new primers produce an amplicon from the input alignment
    publishDir "output/mfeprimer"

    debug true
    
    input:
    path in_alignment
    path new_split

    output:
    path("${in_alignment.simpleName}_*.txt"), emit: hits_file

    script:
    """
    ~/repositories/MFEprimer/mfeprimer-3.3.1-linux-amd64 index -i ${in_alignment}
    ~/repositories/MFEprimer/mfeprimer-3.3.1-linux-amd64 spec -i ${new_split} -d ${in_alignment} -c 6 -o ${in_alignment.simpleName}_${new_split.simpleName}.txt
    """
}

process countCovprop {
    conda "${HOME}/miniconda3/envs/code_hole"
    publishDir "output/"

    debug true

    input:
    path hits_file

    output:
    path("${hits_file.simpleName}_ref_amplicons.csv"), emit: ref_counts

    script:
    """
    python3 ${projectDir}/resources/scripts/primer_counter_waspp_ref_counter.py ${hits_file} ${hits_file.simpleName}
    """
}

process checkCovprop {
    //Probably a python script to do this?
    //Plan to check how many unique refs are in the XXX_ref_amplicons.csv file, then pop those out of the original alignment file before running the primer picking pipeline again.
    //Going to check the idea of "clustering" seqs before primer picking too.
    //Potentially complicated by the fact that there could be multiple potential primer pairs from the primer picking stage - what's a sensible way to handle?
    //Could check if the multiple primers "add" anything, ie cover more species.
    conda "${HOME}/miniconda3/envs/code_hole"
    publishDir "output/", mode: 'copy'

    debug true

    input:
    path in_alignment 
    path in_refCounts
    val iteration

    output:
    path("updated_alignment_*.aln")
    path("Coverage_poportion_*.csv")


    script:
    """
    python3 ${projectDir}/resources/scripts/Primer_coverage_counter.py ${in_alignment} ${iteration} ${in_refCounts}
    """
}

process treeBuilder {
    //This process assumes you already have a family tree built
    conda "${HOME}/miniconda3/envs/code_hole"
    publishDir "output/trees"

    debug true

    input:
    path in_tree
    path hits_file
    path in_meta
    path in_host

    output:
    path("${hits_file.simpleName}*")

    script:
    //This won't work anymore, need to adust inputs.
    """
    python3 ${projectDir}/resources/scripts/primer_counter_waspp_ref_counter.py ${hits_file} ${hits_file.simpleName}
    python3 ${projectDir}/resources/scripts/tree_builder_script.py ${in_tree} ./${hits_file.simpleName}_ref_amplicons.csv ${in_meta} ${in_host} ${hits_file.simpleName}
    """
}