process PARSE_RESULTS {
    tag "$meta"
    label 'process_single'
    cache false

    conda "anaconda::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    // [ ...
    //      [ sample_id, [ ale_res ], [reapr_res], [busco_re_summary], [quast_res]],
    //      [ sample_id_B, [ ale_res_B ], [reapr_res_B], [busco_re_summary_B], [quast_res_B]]
    // ... ]
    input:
    val meta
    val ale_res
    val reapr_res
    val busco_re_summary
    val quast_res
    val merfin_qv_res
    val merfin_completeness_res

    output:
    tuple val(meta), path('table_data_mqc.out'), emit: res
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/assemblyeval/bin/
    def args = task.ext.args ?: ''
    def is_busco = params.busco ? '--busco ' : ''
    """
    parse_results.py \\
        $is_busco\\
        --genomes_ids "$meta" \\
        --ale_res "$ale_res" \\
        --reapr_res "$reapr_res" \\
        --busco_re_summary "$busco_re_summary" \\
        --quast_res "$quast_res" \\
        --merfin_qv_res "$merfin_qv_res" \\
        --merfin_comp_res "$merfin_completeness_res" \\
        -f table_data_mqc.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
