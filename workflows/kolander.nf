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
                if (row.size() != 5) {
                    log.error "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                    System.exit(1)        
                }
            }
}


// Validate input parameters
WorkflowKolander.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.kraken2_db, params.multiqc_config ]
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
include { KRAKEN2 as KRAKEN2_ORIG                             } from '../modules/local/kraken2'
include { KRAKEN2 as KRAKEN2_SIEVED                           } from '../modules/local/kraken2'
include { KOLANDER as KOLANDER_SCRIPT                         } from '../modules/local/kolander'
include { SEQKIT_GREP                                         } from '../modules/local/seqkit_grep'
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
    ch_kraken2_db         = Channel.value( file("${params.kraken2_db}") )

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ()
    ch_short_reads = INPUT_CHECK.out.raw_short_reads
    ch_long_reads  = INPUT_CHECK.out.raw_long_reads

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

    } else {
         ch_short_reads = ch_short_reads
            .map { meta, reads ->
                    if (!params.single_end){ [ meta, [ reads[0], reads[1] ] ] }
                    else [ meta, [reads] ] }
    }
    ch_reads_for_taxonomy = ch_reads_for_taxonomy.mix(ch_short_reads, ch_long_reads)

    // Pass reads through Kraken2 
    ch_profiles = Channel.empty()
    ch_results_for_kolander = Channel.empty()
    ch_results_for_krona = Channel.empty()
    KRAKEN2_DB_PREPARATION (
        ch_kraken2_db
    )
    KRAKEN2_ORIG ( 
        ch_short_reads, KRAKEN2_DB_PREPARATION.out.db
    )
    ch_versions = ch_versions.mix(KRAKEN2_ORIG.out.versions)
    ch_profiles = ch_profiles.mix(KRAKEN2_ORIG.out.report)
    ch_results_for_kolander = ch_results_for_kolander.mix(KRAKEN2_ORIG.out.kraken_output)

    // Sieve reads through the Kolander program
    ch_seived_readids = Channel.empty()
    ch_seived_reports = Channel.empty()
    KOLANDER_SCRIPT (
        ch_results_for_kolander, params.taxids
    )
    ch_versions = ch_versions.mix(KRAKEN2_ORIG.out.versions)
    ch_seived_reports = ch_seived_reports.mix(KOLANDER_SCRIPT.out.filtered_report)
    ch_seived_readids = ch_seived_readids.mix(KOLANDER_SCRIPT.out.results_for_seqkit)

    // Separate out FastQ files based on Kolander results using SeqKit
    ch_seived_readids_and_reads = ch_seived_readids.join(ch_reads_for_taxonomy)
    SEQKIT_GREP (
        ch_seived_readids_and_reads
    )
    ch_versions = ch_versions.mix(KRAKEN2_ORIG.out.versions)
    
    // Re-run Kraken2 for final outputs
    ch_profiles_final = Channel.empty()
    KRAKEN2_SIEVED ( 
        SEQKIT_GREP.out.filtered_fastx, KRAKEN2_DB_PREPARATION.out.db
    )
    ch_profiles_final = ch_profiles_final.mix(KRAKEN2_SIEVED.out.report)
    ch_results_for_krona = ch_results_for_krona.mix(KRAKEN2_SIEVED.out.results_for_krona)

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
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_SIEVED.out.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRONA.out.html.collect{it[1]}.ifEmpty([]))


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
