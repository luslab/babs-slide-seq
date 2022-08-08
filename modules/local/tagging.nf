import java.nio.file.Paths

process MARK_DUPLICATES {
    label "sequencing"
    label 'process_high'

    tag { "${name}" }

    publishDir Paths.get( params.outdir ),
        mode: "copy",
        overwrite: "true",
        saveAs: { filename ->
            if ( filename.indexOf(".txt") != -1 )
            {
                "qc/${filename}"
            }
            else
            {
                "files/${filename}"
            }
        }

    input:
        tuple val(metadata), path(bam)
    
    output:
        tuple val(metadata), path("${name}.dup.bam"), emit: bam
        tuple val(metadata), path("${name}.dup.txt"), emit: metrics
    
    script:
        
        name = metadata["name"]

        """
        picard-tools \
            MarkDuplicates \
                I=$bam \
                O=${name}.dup.bam \
                M=${name}.dup.txt
        """
}

process TAG_BAM {
    label "sequencing"
    label "process_low"
    tag { "${name}" }

    input:
    tuple val(metadata), path(bam)

    output:
    tuple val(metadata), path("*.bam"), path("*.bam.bai")

    script:
    name = metadata["name"]
    suffix  = task.ext.suffix ?: 'tagged'
    """
    echo "Tagging..."
    tag_bam $bam "${name}.${suffix}.bam"
    echo "Indexing..."
    samtools index "${name}.${suffix}.bam"
    """
}

process UMIS_PER_BARCODE {
    label "sequencing"
    label "process_low"
    tag { "${name}" }

    input:
    tuple val(metadata), path(bam), path(bai)

    output:
    tuple val(metadata), path("*.bam"), path("*.bam.bai")

    script:		
    name = metadata["name"]
    suffix = "umis"
    threshold = params.umis_threshold

    """
    umis_per_barcode --threshold $threshold $bam "${name}.${suffix}.bam"
    echo "Indexing..."
    samtools index "${name}.${suffix}.bam"
    """
}

process HTSEQ {
    label "sequencing"
    label "process_low"
    tag { "${name}" }
    
    input:
    tuple val(metadata), path(bam), path(bai)
    path gtf

    output:
    tuple val(metadata), path("${out_bam}"), path("${out_bam}.bai"), emit: bam
    tuple val(metadata), path("${out_txt}"), emit: txt

    script:		
    name = metadata["name"]
    out_sam = "${name}.htseq.sam"
    out_bam = "${name}.htseq.bam"
    out_txt = "${name}.htseq.txt"

    """
    htseq-count --format bam --samout sample.sam $bam $gtf > $out_txt
    samtools view -H $bam > $out_sam
    cat sample.sam >> $out_sam
    samtools view -S -b $out_sam > $out_bam
    samtools index $out_bam
    """
}

process SELECT {
    label "sequencing"
    tag { "${name}" }

    input:
    tuple val(metadata), path(bam)

    output:
    tuple val(metadata), path("*.bam"), path("*.bam.bai"), emit: bam
    tuple val(metadata), val("unique"), path("*.unique.csv"), emit: unique_reads
    tuple val(metadata), val("resolved"), path("*.resolved.csv"), emit: resolved_reads
    tuple val(metadata), val("unresolved"), path("*.unresolved.csv"), emit: unresolved_reads

    script:	
    name    = metadata["name"]
    suffix  = task.ext.suffix ?: 'NO_SUFFIX'
    """
    select "${name}.${suffix}" ${bam} "${name}.${suffix}.bam"
    echo "Indexing..."
    samtools index "${name}.${suffix}.bam"
    """
}

