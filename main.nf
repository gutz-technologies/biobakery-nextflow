#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Workflow imports ───────────────────────────────────────────────────────
include { MGX }          from './workflows/mgx.nf'
include { MTX }          from './workflows/mtx.nf'
include { MGX_MTX }      from './workflows/mgx_mtx.nf'
include { SIXTEENS }     from './workflows/sixteens.nf'
include { VIS }          from './workflows/vis.nf'
include { STATS }        from './workflows/stats.nf'
include { ASSEMBLY } from './workflows/assembly.nf'

// ── Router ─────────────────────────────────────────────────────────────────
workflow {

    if (params.run_medi) {
        if (params.workflow != 'mgx')
            error "ERROR: run_medi is supported only with workflow mgx; MEDI quantifies food DNA (not ${params.workflow})."
        if (!(params.medi_input in ['humann', 'kneaddata']))
            error "ERROR: medi_input must be humann or kneaddata (got ${params.medi_input})."
        if (params.medi_input == 'kneaddata' && !params.run_qc)
            error "ERROR: medi_input kneaddata requires run_qc true."
        if (params.medi_input == 'humann' && !params.run_functional_profiling)
            error "ERROR: medi_input humann requires run_functional_profiling true."
        if (params.medi_input == 'humann' && params.humann_version != 'humann_v4a')
            error "ERROR: medi_input humann requires humann_version humann_v4a."
        def missing = ['medi_db', 'medi_food_matches', 'medi_food_contents'].findAll { !params[it] }
        if (missing)
            error "ERROR: run_medi requires: " + missing.collect { "--${it}" }.join(', ')
    }

    // With medi_input humann, MEDI consumes HUMAnN's retained nucleotide-unaligned reads.
    if (params.run_medi && params.medi_input == 'humann' && (params.humann_options ?: '').contains('--remove-temp-output')) {
        error "ERROR: run_medi requires HUMAnN temp reads; remove --remove-temp-output from humann_options."
    }

    // Validate required params for read-based workflows. mgx_mtx is the
    // exception: it takes two input folders instead of one, and validates them
    // itself, mirroring wmgx_wmtx.py replacing --input with --input-metagenome
    // and --input-metatranscriptome.
    if (!params.readsdir && params.workflow in ['mgx', 'mtx', '16s', 'assembly']) {
        error "ERROR: --readsdir is required. Example: --readsdir /path/to/fastqs"
    }

    // --check_inputs_only runs the read-based workflow as far as its input
    // check and stops there (see subworkflows/check_inputs.nf). It runs inside
    // the workflow rather than beside it so the check tasks carry the same
    // names, and the run that follows with -resume reuses them.
    if (params.check_inputs_only && !(params.workflow in ['mgx', 'mtx', 'mgx_mtx', 'assembly'])) {
        error "ERROR: --check_inputs_only applies to the read-based workflows " +
              "(mgx | mtx | mgx_mtx | assembly), not '${params.workflow}'."
    }
    if (params.check_inputs_only && !params.check_inputs) {
        error "ERROR: --check_inputs_only true needs --check_inputs true."
    }

    switch (params.workflow) {

        case 'mgx':
            MGX()
            break

        case 'mtx':
            MTX()
            break

        case 'mgx_mtx':
            MGX_MTX()
            break

        case '16s':
            SIXTEENS()
            break

        case 'vis':
            VIS(Channel.value(report_input_dir(params.vis_input, 'vis')))
            break

        case 'stats':
            STATS(Channel.value(report_input_dir(params.stats_input, 'stats')))
            break

        case 'assembly':
            ASSEMBLY()
            break

        default:
            error "Unknown workflow '${params.workflow}'. Choose: mgx | mtx | mgx_mtx | 16s | vis | stats | assembly"
    }
}

// The folder a standalone vis or stats run reports on. Chained runs get the
// folder stage_report_input builds instead, so this is only for --workflow
// vis|stats.
def report_input_dir(setting, label) {
    def dir = file(setting ?: params.outdir)

    if (!dir.exists())
        error "ERROR: ${label} input folder not found: ${dir}\n" +
              "Set --${label}_input to a bioBakery output folder."

    return dir
}
