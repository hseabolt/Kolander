/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check already if long reads or fasta 'reads' are provided
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

if(hasExtension(params.input, "csv")){
    Channel
        .from(file(params.input))
        .splitCsv(header: true)
        .map { row ->
                if (row.size() != 6) {
                    log.error "Input samplesheet contains row with ${row.size()} column(s). Expects 6."
                    System.exit(1)        
                }
            }
}


// Validate input parameters
WorkflowKolander.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { POOL_SINGLE_READS as POOL_SHORT_SINGLE_READS        } from '../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                                   } from '../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                } from '../modules/local/pool_single_reads'
include { POOL_FASTA                                          } from '../modules/local/pool_fasta'
include { KRAKEN2_DB_PREPARATION                              } from '../modules/local/kraken2_db_preparation'
include { KRAKEN2                                             } from '../modules/local/kraken2'
include { KRONA_DB                                            } from '../modules/local/krona_db'
include { KRONA                                               } from '../modules/local/krona'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                         } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow KOLANDER {

    ch_versions           = Channel.empty()
    ch_reads_for_taxonomy = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ()
    ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
    ch_raw_long_reads  = INPUT_CHECK.out.raw_long_reads
    ch_fasta           = INPUT_CHECK.out.fasta
    ch_versions        = ch_versions.mix(INPUT_CHECK.out.versions)

    // Concatenate reads if required
    if ( params.concatenate_reads ) {
        // short reads
        // group and set group as new id
        ch_short_reads_grouped = ch_short_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                    def meta = [:]
                    meta.id          = "group-$group"
                    meta.group       = group
                    meta.single_end  = params.single_end
                    if (!params.single_end) [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
                    else [ meta, reads.collect { it }, [] ]
                }
        // long reads
        // group and set group as new id
        ch_long_reads_grouped = ch_long_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                def meta = [:]
                meta.id          = "group-$group"
                meta.group       = group
                [ meta, reads.collect { it } ]
            }
        // fasta files
        // group and set group as new id
        ch_fasta_grouped = ch_fasta
            .map { meta, fasta -> [meta.group, meta, fasta ] }
            .groupTuple(by: 0)
            .map { group, metas, fastas -> 
                def meta = [:]
                meta.id = "group-$group"
                meta.group = group
                [ meta, fastas.collect { it } ]
            }

        // pool short reads
        if ( params.single_end ){
            POOL_SHORT_SINGLE_READS ( ch_short_reads_grouped )
            ch_short_reads = POOL_SHORT_SINGLE_READS.out.reads
        } else {
            POOL_PAIRED_READS ( ch_short_reads_grouped )
            ch_short_reads = POOL_PAIRED_READS.out.reads
        }

        // long reads
        POOL_LONG_READS ( ch_long_reads_grouped )
        ch_long_reads = POOL_LONG_READS.out.reads

        // fasta files
        POOL_FASTA ( ch_fasta_grouped )
        ch_fasta = POOL_FASTA.out.fasta

    } else {
         ch_short_reads_grouped = ch_short_reads
            .map { meta, reads ->
                    if (!params.single_end){ [ meta, [reads[0]], [reads[1]] ] }
                    else [ meta, [reads], [] ] }
    }

    ch_reads_for_taxonomy = ch_reads_for_taxonomy.mix(

    )

    

    // Pass reads through Kraken2 
    KRAKEN2_DB_PREPARATION (
        ch_kraken2_db
    )
    KRAKEN2 ( 
        ch_reads_for_taxonomy, KRAKEN2_DB_PREPARATION.out.db
    )
    ch_versions = ch_versions.mix(KRAKEN2.out.versions)
    ch_profiles = ch_profiles.mix(KRAKEN2.out.report)
    ch_results_for_krona = ch_results_for_krona.mix(KRAKEN2.out.results_for_krona)

    // Generate a Krona HTML plot
    if ( !params.skip_krona ){
        KRONA_DB ()
        KRONA (
            ch_results_for_krona,
            KRONA_DB.out.db.collect()
        )
        ch_versions = ch_versions.mix(KRONA.out.versions)
    }


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowKolander.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowKolander.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

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
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
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
