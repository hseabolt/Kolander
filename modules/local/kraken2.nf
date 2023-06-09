// Credit: nf-core/mag (modified here)
process KRAKEN2 {
    tag "${meta.id}-${db_name}"
    label 'process_high'
    label 'process_high_memory'

    conda "bioconda::kraken2=2.0.8_beta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.0.8_beta--pl526hc9558a2_2' :
        'quay.io/biocontainers/kraken2:2.0.8_beta--pl526hc9558a2_2' }"

    input:
    tuple val(meta), path(reads)
    tuple val(db_name), path("database/*")

    output:
    tuple val(meta), path("*results.krona")                , emit: results_for_krona
    tuple val(meta), path("*kraken2_report.txt")           , emit: report
    tuple val(meta), path("*.kraken2.kraken")              , emit: kraken_output
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gzip = reads[0].endsWith('gz') ? "--gzip-compressed" : ""
    def input = meta.single_end ? "\"${reads}\"" :  "--paired \"${reads[0]}\" \"${reads[1]}\""
    """
    kraken2 \
        --report-zero-counts \
        --threads ${task.cpus} \
        --db database \
        --report ${prefix}.kraken2_report.txt \
        --output ${prefix}.kraken2.kraken \
        $gzip \
        $input
        
    cat ${prefix}.kraken2.kraken | cut -f 2,3 > ${prefix}.results.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //' | sed 's/ Copyright.*//')
    END_VERSIONS
    """
}