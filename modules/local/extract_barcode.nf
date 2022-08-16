process EXTRACT_BARCODE {
    label "sequencing"
    label "ultra_low"
    tag { "${name}" }

    input:
    tuple val(metadata), path(read1), path(read2)

    output:
    tuple val(metadata), path("${name}.fastq.gz"), emit: fastq
    tuple val(metadata), path("${name}.extract_barcode.csv"), emit: metrics
    tuple val(metadata), path("${name}.extract_barcode.up_distances.csv"), emit: distances

    script:		
        
        name = metadata["name"]
        read_structure = metadata["read_structure"]

        """
        extract_barcode \
            --sample "${name}" \
            --fastq "${name}.fastq.gz" \
            --read-structure "${read_structure}" \
            --max-distance ${params.up_errors_threshold} \
            $read1 $read2
        """
}

