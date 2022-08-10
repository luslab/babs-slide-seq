import java.nio.file.Paths

process STAR {
    tag { "${name}" }
    label 'process_high'
    container = 'quay.io/biocontainers/star:2.6.1a--1'

    publishDir Paths.get( params.outdir ),
        mode: "copy",
        overwrite: "true",
        saveAs: { filename ->
            if ( filename.indexOf(".Log.") != -1 )
            {
                "qc/align/${filename}"
            }
            else
            {
                "files/align/${filename}"
            }
        }

    input:
    tuple val(metadata), path(fastq)
    path(index)
    
    output:
    tuple val(metadata), path("${name}.Log.*"), emit: log
    tuple val(metadata), path("${name}.SJ.out.tab"), emit: tab
    tuple val(metadata), path("${name}.Aligned.sortedByCoord.out.bam"), emit: bam

    script:		
    name = metadata["name"]

    """
    STAR \
        --runMode alignReads \
        --runThreadN $task.cpus \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes NH HI AS nM NM MD \
        --genomeDir $index \
        --outFileNamePrefix "${name}." \
        --readFilesCommand zcat \
        --readFilesIn $fastq
    """
}

