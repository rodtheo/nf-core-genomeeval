process REAPR {
    tag "$meta.id"
    label 'process_high'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.5.0--pyhdfd78af_0':
        'rodtheo_genomics:eval_assem_ale_reapr' }"

    input:
    tuple val(meta_asm), path(asm)        // Required:    One genome assembly to evaluate
    tuple val(meta), path(bam)            // Required:    Corresponding aligned file
    

    output:
    tuple val(meta), path("*-REAPR"), emit: reapr_dir
    tuple val(meta), path("*-REAPR/05.summary.report.txt"), emit:reapr_summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def prefix = task.ext.prefix ?: "${meta_asm.id}"
    def args   = task.ext.args ?: ''
    """
    reapr pipeline \\
        $asm \\
        $bam \\
        ${prefix}-REAPR

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        REAPR: XXX)
    END_VERSIONS
    """
}
