process PLOT {
    label 'plot'
    label 'ultra_low'
    tag "$name"

    input:
    tuple val(metadata), path(file), val(value)

    output:
    tuple val(metadata), file("*.png"), emit: png
    tuple val(metadata), file("*.pdf"), emit: pdf

    script:
    name = metadata["name"]

    def args    = task.ext.args ?: ''
    def script  = task.ext.script ?: ''
    def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
    """
    $script $args $file $value ${name}.${suffix}
    """
}

process PLOT_HAMMING_HISTO {
    label 'plot'
    label 'ultra_low'
    tag "$name"

    input:
    tuple val(metadata), path(ordered), path(shuffled)

    output:
    tuple val(metadata), file("*.png"), emit: png
    tuple val(metadata), file("*.pdf"), emit: pdf

    script:
    name = metadata["name"]

    def args    = task.ext.args ?: ''
    def script  = task.ext.script ?: ''
    def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
    """
    $script $args $ordered $shuffled ${name}.${suffix}
    """
}

process PLOT_SPATIAL_UMI {
    label 'plot'
    label 'ultra_low'
    tag "$name"

    input:
    tuple val(metadata), path(mtx), path(coords)

    output:
    tuple val(metadata), file("*.png"), emit: png
    tuple val(metadata), file("*.pdf"), emit: pdf

    script:
    name = metadata["name"]

    def args    = task.ext.args ?: ''
    def script  = task.ext.script ?: ''
    def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
    """
    $script $args $mtx $coords ${name}.${suffix}
    """
}
