#!/usr/bin/env nextflow

//Get the modules for the pipeline
include { pickPrimers; splitPrimers; checkFamilycov; treeBuilder } from './modules/main_pipeline.nf'

workflow primer_picker_wf {
    //Define the input channels
    inAlignment_ch = Channel.value("${params.alignment}")
    //Run varvamp on input alignments
    pickPrimers(inAlignment_ch)
    //Split primer output to pairs
    splitPrimers(pickPrimers.out.new_primers)
    //mfeprimer index and spec on generated primers
    checkFamilycov(inAlignment_ch, splitPrimers.out.new_split.flatten())
    //slurp in relevent files
    inTree_ch = Channel.value("${params.tree}")
    inTreemeta_ch = Channel.value("${projectDir}/resources/metadata/all_names.txt")
    inHosts_ch = Channel.value("${projectDir}/resources/metadata/host_accession.tsv")
    treeBuilder(inTree_ch, checkFamilycov.out.hits_file, inTreemeta_ch, inHosts_ch)
}

workflow {
    primer_picker_wf()
}