#!/usr/bin/env nextflow

//Get the modules for the pipeline
include { pickPrimers; splitPrimers; checkFamilycov; countCovprop; checkCovprop; treeBuilder } from './modules/main_pipeline.nf'

workflow primer_picker_wf {
    //Define the input channels
    inAlignment_ch = Channel.value("${params.alignment}")
    inAmbig_ch = Channel.value("${params.ambig}")
    inIteration_ch = Channel.value("${params.iter}")
    //Run varvamp on input alignments
    pickPrimers(inAmbig_ch, inAlignment_ch)
    //Split primer output to pairs
    splitPrimers(pickPrimers.out.new_primers)
    //mfeprimer index and spec on generated primers
    checkFamilycov(inAlignment_ch, splitPrimers.out.new_split.flatten())
    countCovprop(checkFamilycov.out.hits_file)
    checkCovprop(inAlignment_ch, countCovprop.out.ref_counts.collect(), inIteration_ch)
    
}

workflow tree_builder_wf {
    inHosts_ch = Channel.value("${projectDir}/resources/metadata/host_accession_curated.tsv")
    inTree_ch = Channel.value("${params.tree}")
    inTreemeta_ch = Channel.value("${projectDir}/resources/metadata/all_names.txt")
    //Need to thing of a way to five the final hits_file to this process (could just run the python script separately?)
    treeBuilder(inTree_ch, checkFamilycov.out.hits_file, inTreemeta_ch, inHosts_ch)
}

workflow {
    primer_picker_wf()
}