#!/usr/bin/env nextflow

//Get the modules for the pipeline
include { pickPrimers; checkFamilycov; treeBuilder } from './modules/main_pipeline.nf'

workflow primer_picker_wf {
    //Define the input channels
    inAlignment_ch = Channel.value("${params.alignment}")
    //Eun varvamp on input alignments
    pickPrimers(inAlignment_ch)
    //mfeprimer index and spec on generated primers
    checkFamilycov(inAlignment_ch, pickPrimers.out.new_primers)
    //slurp in relevent files
    inTree_ch = Channel.value("${params.tree}")
    inTreemeta_ch = Channel.value("${projectDir}/resources/metadata/all_names.txt")
    inHosts_ch = Channel.value("${projectDir}/resources/metadata/host_accession.tsv")
    treeBuilder(inTree_ch, checkFamilycov.out.hits_file, inTreemeta_ch, inHosts_ch)
}

workflow {
    primer_picker_wf()
}