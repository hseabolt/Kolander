process SEQKIT_GREP {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::seqkit=2.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_02' }"

    input:
    tuple val(meta), path(list_of_ids), path(reads)

    output:
    tuple val(meta), path("*.fast*.gz")                    , emit: filtered_fastx
    path "versions.yml"                                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: ${reads[0]}.simpleName
    def args = task.ext.args ?: ""
    def extension = ${reads[0]}.endsWith('a') ? "fasta" : "fastq"
    if ( !${reads[1]} ) {
        """
        seqkit grep \\
            --pattern-file ${list_of_ids} \\
            --threads ${task.cpus} \\
            $args \\
            --out-file ${prefix}.${extension}.gz \\
            ${input}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
    }   else    {
        """
        seqkit grep \\
            --pattern-file ${list_of_ids} \\
            --threads ${task.cpus} \\
            $args \\
            --out-file ${prefix}_1.${extension}.gz \\
            ${reads[0]}

        seqkit grep \\
            --pattern-file ${list_of_ids} \\
            --threads ${task.cpus} \\
            $args \\
            --out-file ${prefix}_2.${extension}.gz \\
            ${reads[1]}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
    }
}