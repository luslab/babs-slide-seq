import java.nio.file.Paths

process fastqc {

	label "sequencing"
	//label "quality_control"
	
	tag { "${name}" }

	publishDir Paths.get( params.out_dir , "qc" , "fastqc" ),
		mode: "copy",
		overwrite: "true"

	input:
		tuple val(metadata), path(fastq1), path(fastq2)

	output:
		tuple val(metadata), file("*_fastqc.html"), emit: html
		tuple val(metadata), file("*_fastqc.zip"), emit: zip

	script:		

		name = metadata["name"]
		
		"""
		fastqc $fastq1 $fastq2
		"""
}

// number of reads for each the barcodes matched or not, mapped or not
process reads_per_barcode_fastq {

	label "python_2"
	label "quality_control"
	
	tag { "${name}" }

	input:
		tuple val(metadata), path(fastq), path(script)

	output:
		tuple val(metadata), file("${name}.reads_per_barcode_fastq.csv")

	script:		

		name = metadata["name"]
		
		"""
		python3 $script $fastq "${name}.reads_per_barcode_fastq.csv"
		"""
}

//// the number of reads UP-matched and mapped
//process reads_up_matching {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.reads_up_matching.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view --expr '[ls] == "LONG_ENOUGH"' $bam \
//			| sed 's/^\\([^\\t]\\+\\).*us:Z:\\([A-Z]\\+\\).*as:Z:\\([A-Z]\\+\\).*/\\1\\t\\2\\t\\3/g' \
//			| sort \
//			| uniq \
//			| awk '{ print \$2 "\\t" \$3 }' \
//			| sort \
//			| uniq -c \
//			| sort -rn \
//			| awk '{ printf "%s,%s,%s\\n", \$2, \$3, \$1 }' \
//			> "${name}.reads_up_matching.csv"
//		"""
//}
//
//// the number reads of UP-matchec reads that were barcode-matched and mapped 
//process reads_barcode_matching {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.reads_barcode_matching.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view --expr '[us] == "MATCHED"' $bam \
//			| sed 's/^\\([^\\t]\\+\\).*as:Z:\\([A-Z]\\+\\).*bs:Z:\\([A-Z]\\+\\).*/\\1\\t\\2\\t\\3/g' \
//			| sort \
//			| uniq \
//			| awk '{ print \$3 "\\t" \$2 }' \
//			| sort \
//			| uniq -c \
//			| sort -rn \
//			| awk '{ printf "%s,%s,%s\\n", \$2, \$3, \$1 }' \
//			> "${name}.reads_barcode_matching.csv"
//		"""
//}
//
//// number of UMIs for each barcode matched, mapped and selected, i.e. the
//// barcodes used for final counting
//process reads_per_barcode_umi {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.reads_per_barcode_umi.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view \
//			--expr '[bs] == "MATCHED" && ([cs] == "UNIQUE" || [cs] == "INCLUDED")' \
//			$bam \
//			| sed 's/^\\([^\\t]\\+\\).*bc:Z:\\([A-Z]\\+\\).*mi:Z:\\([A-Z]\\+\\).*/\\2 \\3 \\1/g' \
//			| sort \
//			| cut -d " " -f 1,2 \
//			| uniq -c \
//			| sort -nr \
//			| awk '{ printf "%s,%s,%s\\n", \$2, \$3, \$1 }' \
//			> "${name}.reads_per_barcode_umi.csv"
//		"""
//}
//
//// number of dulicates reads for barcode-matched reads
//process count_duplicates {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.count_duplicates.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view \
//			--expr '[bs] == "MATCHED" && [as] == "MAPPED"' \
//			$bam \
//			| sed 's/^\\([^\\t]\\+\\).*ds:Z:\\([A-Z]\\+\\).*/\\1 \\2/g' \
//			| sort \
//			| uniq \
//			| awk '{ print \$2 }' \
//			| sort \
//			| uniq -c \
//			| sort -rn \
//			| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
//			> "${name}.count_duplicates.csv"
//		"""
//}
//
//// number of genome mappings for barcode-matched and primary reads
//process count_mappings {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.count_mappings.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view \
//			--expr '[bs] == "MATCHED" && [as] == "MAPPED" && [ds] == "PRIMARY"' \
//			$bam \
//			| awk '{ print \$1 }' \
//			| sort \
//			| uniq -c \
//			| awk '{ print \$1 }' \
//			| sort \
//			| uniq -c \
//			| sort -rn \
//			| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
//			> "${name}.count_mappings.csv"
//		"""
//}
//
//// count the number of resolved reads among the multi mapped reads
//process count_resolved {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.count_resolved.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view \
//			--expr '[as] == "MAPPED" && [bs] == "MATCHED"' \
//			$bam \
//			| sed 's/^\\([^\\t]\\+\\).*mm:Z:\\([A-Z]\\+\\).*/\\1 \\2/g' \
//			| sort \
//			| uniq \
//			| awk '{ print \$2 }' \
//			| sort \
//			| uniq -c \
//			| sort -rn \
//			| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
//			> "${name}.count_resolved.csv"
//		"""
//}
//
//// count the number of resolved reads among the reads that a UMI mapping
//// different genes
//process count_select {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.count_select.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view \
//			--expr '[bs] == "MATCHED" && ( [mm] == "UNIQUE" || [mm] == "INCLUDED" )' \
//			$bam \
//			| sed 's/^.*cs:Z:\\([A-Z]\\+\\).*/\\1/g' \
//			| sort \
//			| uniq -c \
//			| sort -rn \
//			| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
//			> "${name}.count_select.csv"
//		"""
//}
//
//// number of regions for barcodes matched, mapped and selected, i.e. the
//// barcodes used for final counting
//process reads_per_function {
//
//	label "samtools"
//	label "quality_control"
//	
//	tag { "${name}" }
//
//	input:
//		tuple val(metadata), path(bam)
//
//	output:
//		tuple val(metadata), file("${name}.reads_per_function.csv")
//
//	script:		
//
//		name = metadata["name"]
//		
//		"""
//		samtools view \
//			--expr '[bs] == "MATCHED" && ([cs] == "UNIQUE" || [cs] == "INCLUDED")' \
//			$bam \
//			| tr "\\t" "\\n" \
//			| grep "gf:Z:" \
//			| sed 's/gf:Z://g' \
//			| sort \
//			| uniq -c \
//			| sort -rn \
//			| awk '{ printf "%s,%s\\n", \$2, \$1 }' \
//			> "${name}.reads_per_function.csv"
//		"""
//}

