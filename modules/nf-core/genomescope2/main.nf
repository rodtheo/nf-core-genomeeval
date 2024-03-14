process GENOMESCOPE2 {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomescope2:2.0--py310r41hdfd78af_5':
        'biocontainers/genomescope2:2.0--py310r41hdfd78af_5' }"

    input:
    tuple val(meta), path(histogram)

    output:
    tuple val(meta), path("*_linear_plot.png")            , emit: linear_plot_png
    tuple val(meta), path("*_transformed_linear_plot.png"), emit: transformed_linear_plot_png
    tuple val(meta), path("*_log_plot.png")               , emit: log_plot_png
    tuple val(meta), path("*_transformed_log_plot.png")   , emit: transformed_log_plot_png
    tuple val(meta), path("*_model.txt")                  , emit: model
    tuple val(meta), path("*_summary.txt")                , emit: summary
    tuple val(meta), path("lookup_table.txt")             , emit: lookup_table
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def kmer = task.ext.kmer_size ?: "${meta.kmer_size}"
    def ploidy = task.ext.ploidy ?: "${meta.ploidy}"
    """
    genomescope2 \\
        --input $histogram \\
        --kmer_length $kmer \\
        --ploidy $ploidy \\
        --fitted_hist \\
        $args \\
        --output . \\
        --name_prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    '${task.process}':
        genomescope2: \$( genomescope2 -v | sed 's/GenomeScope //' )
    END_VERSIONS
    """
}
