process HEADER_FASTA_REAPR {
    tag "$meta"
    label 'process_single'

    conda "bioconda::biopython=1.81"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81' :
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta), path(asm)

    output:
    tuple val(meta), path('*_header_valid.fasta'), emit: asm
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    check_header_fasta.py \\
        $asm \\
        ${prefix}_header_valid.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
