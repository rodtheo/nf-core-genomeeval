process COMPLEASM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/compleasm:0.2.5--pyh7cba7a3_0 ':
        'biocontainers/compleasm:0.2.5--pyh7cba7a3_0 ' }"

    input:
    tuple val(meta), path(asm)
    // val mode                              // Required:    One of genome, proteins, or transcriptome
    val lineage                           // Required:    lineage to check against, "auto" enables --auto-lineage instead
    // path busco_lineages_path              // Recommended: path to busco lineages - downloads if not set
    // path config_file                      // Optional:    busco configuration file
    

    output:
    tuple val(meta), path("*-compleasm"), emit: busco_dir
    tuple val(meta), path("*-compleasm/summary.txt"), emit: short_summaries_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def args   = task.ext.args ?: ''
    // def busco_config = config_file ? "--config $config_file" : ''
    def compleasm_lineage = lineage.equals('auto') ? '--auto-lineage' : "-l ${lineage}"
    // def busco_lineage_dir = busco_lineages_path ? "--download_path ${busco_lineages_path}" : ''
    """
    compleasm run \\
        -t$task.cpus \\
        $compleasm_lineage \\
        $args \\
        -a $asm \\
        -o ${prefix}-compleasm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compleasm: 0.2.5
    END_VERSIONS
    """
}
