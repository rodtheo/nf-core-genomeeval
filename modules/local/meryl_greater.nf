process MERYL_GREATER_THAN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::meryl=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.3--h87f3376_1':
        'biocontainers/meryl:1.3--h87f3376_1' }"

    input:
    tuple val(meta), path(in_meryl_db)

    output:
    tuple val(meta), path("*.meryl"), emit: meryl_db
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
        meryl greater-than 1\\
            threads=$task.cpus \\
            $args \\
            $in_meryl_db \\
            output ${prefix}.gt1.meryl


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}