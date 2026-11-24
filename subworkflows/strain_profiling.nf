#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { sample2markers } from '../modules/strainphlan/main.nf'
include { strainphlan }    from '../modules/strainphlan/main.nf'

workflow STRAIN_PROFILING {
    take:
    sam_bzip  // [sample, .sam.bz2]

    main:
    def clades = (params.strainphlan_clades ?: '').split(',', -1).collect { it.trim() }
    if (clades.any { !(it ==~ /t__SGB[0-9]+/) } || clades.unique(false) != clades)
        error 'StrainPhlAn requires distinct t__SGB IDs in --strainphlan_clades (e.g. t__SGB15318).'
    if (!(params.strainphlan_phylophlan_mode in ['fast', 'accurate']))
        error 'StrainPhlAn requires fast or accurate in --strainphlan_phylophlan_mode.'

    def db_name = "${params.metaphlan_index}.pkl"
    def db = file(params.strainphlan_db ?: "${params.metaphlan_db}/${db_name}", glob: false, checkIfExists: true)
    if (!db.isFile() || db.name != db_name)
        error "StrainPhlAn requires the matching ${db_name} beside its Bowtie2 shards."
    // The tree CLI resolves the Bowtie2 shards beside this path. Do not stage only the .pkl.
    def database = Channel.value(db.toAbsolutePath().normalize().toString())

    def references = params.strainphlan_references ? params.strainphlan_references.split(',', -1).collect {
        def path_text = it.trim()
        if (!path_text)
            error 'StrainPhlAn reference FASTA paths must not be empty.'
        def fasta = file(path_text, glob: false, checkIfExists: true)
        if (!fasta.isFile())
            error "StrainPhlAn reference must be a FASTA file: ${fasta}"
        fasta
    } : []
    if (references.collect { it.name }.unique().size() != references.size())
        error 'StrainPhlAn reference FASTA filenames must be distinct.'

    markers_out = sample2markers(sam_bzip, database)
    all_markers = markers_out.markers.map { sample, marker -> marker }.toSortedList().map { markers ->
        if (markers.size() + references.size() < 4)
            error 'StrainPhlAn requires at least four total samples plus references. Use --strainphlan_references to add reference FASTAs.'
        markers
    }
    strainphlan_out = strainphlan(Channel.fromList(clades), all_markers, database, Channel.value(references))

    emit:
    results = strainphlan_out.results
    trees   = strainphlan_out.tree
}
