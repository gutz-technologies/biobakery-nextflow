#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process sample2markers {
    tag "$sample"
    publishDir "${params.outdir}/strainphlan/markers", mode: 'copy'

    input:
    tuple val(sample), path('alignment.sam.bz2')
    val database

    output:
    tuple val(sample), path("${sample}.json.bz2"), emit: markers

    script:
    """
    sample2markers.py -i alignment.sam.bz2 -o . -d "${database}" -n ${task.cpus}
    mv -- alignment.json.bz2 "${sample}.json.bz2"
    """
}

process strainphlan {
    tag "$clade"
    publishDir "${params.outdir}/strainphlan/${clade}", mode: 'copy'

    input:
    val clade
    path markers, stageAs: 'markers/*'
    val database
    path references, stageAs: 'references/*'

    output:
    path 'output/', emit: results
    path "output/RAxML_bestTree.${clade}.StrainPhlAn4.tre", emit: tree
    path 'output/*.tsv', emit: table, optional: true
    tuple val("${task.process}"), val('strainphlan'), eval('strainphlan --version 2>&1 | cut -d" " -f3'), topic: versions

    script:
    def ref_flag = references ? '-r references/*' : ''
    def extra_args = params.strainphlan_options ?: ''
    """
    mkdir -p output
    strainphlan \\
        -s markers/* \\
        -d "${database}" \\
        -o output \\
        -n ${task.cpus} \\
        -c ${clade} \\
        --phylophlan_mode ${params.strainphlan_phylophlan_mode} \\
        ${ref_flag} ${extra_args}
    test -s "output/RAxML_bestTree.${clade}.StrainPhlAn4.tre"
    """
}
