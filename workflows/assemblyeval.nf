/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowAssemblyeval.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPARE_INPUT } from '../subworkflows/local/prepare_input'
include { READ_MAPPING } from '../subworkflows/local/read_mapping'
include { CONTIGUITY_ASM } from '../subworkflows/local/contiguity'
include { COMPLETENESS_ASM } from '../subworkflows/local/completeness'
include { CORRECTNESS_ASM } from '../subworkflows/local/correctness'
include { KMER_PROFILE } from '../subworkflows/local/kmerprofile'
include { CONTAMINATION_ASM } from '../subworkflows/local/contamination'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PARSE_RESULTS } from '../modules/local/parse_results'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ASSEMBLYEVAL {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    PREPARE_INPUT (
        Channel.fromPath(params.input)
    )
    
    ch_versions = ch_versions.mix(PREPARE_INPUT.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    illumina_ch = PREPARE_INPUT.out.illumina
    // illumina_ch.view()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        illumina_ch
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
     

    //
    // MODULE: Contamination
    //
    // CONTAMINATION_ASM( illumina_ch, PREPARE_INPUT.out.assemblies )

    //
    // MODULE: Completeness Metrics
    //
    COMPLETENESS_ASM (
        illumina_ch, PREPARE_INPUT.out.assemblies
    )

    // MODULE: Contiguity Metrics
    //
    CONTIGUITY_ASM (
        PREPARE_INPUT.out.assemblies
    )

    //
    // MODULE: Correctness Metrics
    //
    // SUB-MODULE: Mapping reads
    //
    READ_MAPPING (
        illumina_ch, PREPARE_INPUT.out.assemblies
    )

    // READ_MAPPING.out.bam.view{ "*** ETAPA 2 - BAM: $it"}
    // READ_MAPPING.out.bai.view{ "*** ETAPA 2 - BAI: $it"}

    
    CORRECTNESS_ASM (
        READ_MAPPING.out.asm, READ_MAPPING.out.bam, READ_MAPPING.out.bai
    )

    KMER_PROFILE (
        illumina_ch, PREPARE_INPUT.out.assemblies
    )

    // PREPARE_INPUT.out.assemblies.view{ "ASSEMBLIES: $it" }
    // parse_results_to_table(args.genomes_ids, args.ale_res, args.reapr_res, args.busco_re_summary, args.quast_res, args.file_out)
    
    out_asm_ch = CORRECTNESS_ASM.out.reapr.map{ meta, reapr -> meta['id'] }.collect( sort: true ).view{ "ASSEMBLIES: $it" }
    out_ale_ch = CORRECTNESS_ASM.out.ale.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_reapr_ch = CORRECTNESS_ASM.out.reapr.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_busco_re_summary_ch = COMPLETENESS_ASM.out.busco_short_summaries_txt.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_quast_ch = CONTIGUITY_ASM.out.quast_tsv.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_merfin_qv_ch = KMER_PROFILE.out.merfin_logs.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }
    out_merfin_completeness_ch = KMER_PROFILE.out.merfin_completeness.toSortedList( { a, b -> a[0]['id'] <=> b[0]['id'] } ).flatMap().collect{ it[1] }

    // out_table_ch = Channel.fromPath("out_table.csv")
    
    COMPLETENESS_ASM.out.busco_short_summaries_txt.view{ "REAPR: $it" }
    
    PARSE_RESULTS ( out_asm_ch, out_ale_ch, out_reapr_ch, out_busco_re_summary_ch, out_quast_ch, out_merfin_qv_ch, out_merfin_completeness_ch  )

    // // CORRECTNESS_ASM.out.reapr.collect( {it[1]}, sort: {it.getName()} ).set{ new_collect }
    

    // ch_versions = ch_versions.mix(FASTQC.out.versions.first()).mix(CONTIGUITY_ASM.out.versions.first())


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowAssemblyeval.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowAssemblyeval.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CONTIGUITY_ASM.out.quast.collect{it[1]}.ifEmpty([]))
    if ( params.busco ) {
        ch_multiqc_files = ch_multiqc_files.mix(COMPLETENESS_ASM.out.busco_short_summaries_txt.collect{it[1]}.ifEmpty([]))
    }
    ch_multiqc_files = ch_multiqc_files.mix( PARSE_RESULTS.out.res.collect{it[1]}.ifEmpty([]))


    // ch_multiqc_files.collect().view{ "QUAST_TSV: $it" }

    ch_multiqc_custom_config.mix( Channel.fromPath("./assets/section_name_with_slash.yml")).set{ ch_multiqc_custom_config }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
        
    sendMail(to: 'theodoro.biotec@gmail.com', from: 'theodoro.biotec@gmail.com', subject: 'My pipeline execution', body: msg)
}

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
