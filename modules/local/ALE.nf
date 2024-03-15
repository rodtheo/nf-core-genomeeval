process ALE {
    tag "$meta.id"
    label 'process_medium'

    // conda "${moduleDir}/environment.yml"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/busco:5.5.0--pyhdfd78af_0':
    //     'rodtheo_genomics:eval_assem_ale_reapr' }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ale:20180904--py27ha92aebf_0':
        'biocontainers/ale:20180904--py27ha92aebf_0' }"

    input:
    tuple val(meta_asm), path(asm)        // Required:    One genome assembly to evaluate
    tuple val(meta), path(bam)            // Required:    Corresponding aligned file
    

    output:
    tuple val(meta_asm), path('*.ale')       , emit: ale
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def prefix = task.ext.prefix ?: "${meta_asm.id}"
    def args   = task.ext.args ?: ''
    def kmer_size   = task.ext.kmer ?: "${meta_asm.kmer_size}"
    """
    ALE \\
        ${args} \\
        --kmer ${kmer_size} \\
        $bam \\
        $asm \\
        ${prefix}.ale

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ALE: 20220503)
    END_VERSIONS
    """
}
