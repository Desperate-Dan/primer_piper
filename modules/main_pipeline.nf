process pickPrimers {
    conda "${HOME}/miniconda3/envs/varvamp"
    publishDir "output/varvamp/"

    debug true

    input:
    path in_alignment
    
    output:
    path("${in_alignment.simpleName}")
    path("${in_alignment.simpleName}/primers.fasta"), emit: new_primers

    script:
    """
    varvamp single -a 10 ${in_alignment} ${in_alignment.simpleName}
    """
}

process checkFamilycov {
    publishDir "output/mfeprimer"

    debug true
    
    input:
    path in_alignment
    path new_primers

    output:
    path("${in_alignment.simpleName}_family_coverage.txt"), emit: hits_file

    script:
    """
    ~/repositories/MFEprimer/mfeprimer-3.3.1-linux-amd64 index -i ${in_alignment}
    ~/repositories/MFEprimer/mfeprimer-3.3.1-linux-amd64 spec -i ${new_primers} -d ${in_alignment} -c 6 -o ${in_alignment.simpleName}_family_coverage.txt
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
    path("primercounter*")

    script:
    """
    python3 ${projectDir}/resources/scripts/primer_counter_waspp_ref_counter.py ${hits_file}
    python3 ${projectDir}/resources/scripts/tree_builder_script.py ${in_tree} ./primercounter_ref_amplicons.csv ${in_meta} ${in_host}
    """
}