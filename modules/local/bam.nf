process BAM_METRICS {
    label "python"
    label "bam_merics"
    label "process_low"
    tag { "${name}" }

    input:
    tuple val(metadata), path(bam), path(bai)

    output:
    tuple val(metadata), path("*.csv")

    script:		
    name = metadata["name"]
    def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
    def script  = task.ext.script ?: ''
    """
    $script $bam ${name}.${suffix}.csv
    """
}

process BAM_FILTER {
    label "samtools"
    label 'process_medium'

    tag { "${name}" }

    input:
    tuple val(metadata), path(bam), path(bai)

    output:
    tuple val(metadata), path("*.bam"), path("*.bam.bai")

    script:		
    name = metadata["name"]
    def suffix  = task.ext.suffix ?: 'NO_SUFFIX'
    def expr    = task.ext.expr ?: ''
    """
    echo "Filtering..."
    samtools view --expr '${expr}' --output "${name}.${suffix}.bam" $bam
    echo "Indexing..."
    samtools index "${name}.${suffix}.bam"
    """
}
