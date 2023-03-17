process KOLANDER {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::perl-parallel-forkmanager=2.02"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-parallel-forkmanager:2.02--pl5321hdfd78af_1' :
        'quay.io/biocontainers/perl-parallel-forkmanager:2.02--pl5321hdfd78af_1' }"

    input:
    tuple val(meta), path(report)
    val taxids

    output:
    tuple val(meta), path("*.kolander.txt")                , emit: filtered_report
    tuple val(meta), path("*.readids.txt")                 , emit: results_for_seqkit
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def VERSION = "1.0"  // WARN: Version is not provided by script on the command-line.
    """
    kolander \\
        --input ${report} \\
        --taxids ${taxids} \\
        --threads ${task.cpus} \\
        $args \\
        --output ${prefix}.kolander.txt
    cut -f2 ${prefix}.kolander.txt > ${prefix}.readids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kolander: ${VERSION}
    END_VERSIONS
    """
}