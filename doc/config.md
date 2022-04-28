

# Configuration

```
name,fastq_1,fastq_2,puck,read_structure,genome,gtf
sample_name,/path/to/read1/fastq,/path/to/read2/fastq,/path/to/coordinates/csv,read_structure_string,/path/to/STAR/index,/path/to/gtf
```

 * `name`: the sample name
 * `fastq_1` the path (relative or absolute) to the `FASTQ` file containing reads 1
 * `fastq_2` the path (relative or absolute) to the `FASTQ` file containing reads 2
 * `puck` the path (relative or absolute) to the `CSV` file containing the barcodes coordinates
 * `read_structure` the read structure string
 *  `genome` the path (relative or absolute) to the `STAR` genome index
 *  `gtf`: the path (relative or absolute) to the `GTF` file

The `FASTQ` files with the same sample name will be merged.

```
name,fastq_1,fastq_2,puck,read_structure,genome,gtf
sample1,test/s ample1_L001_R1.fastq.gz,test/sample1_L001_R2.fastq.gz,test/puck1.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample1,test/sample1_L002_R1.fastq.gz,test/sample1_L002_R2.fastq.gz,test/puck1.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample1,test/sample1_L003_R1.fastq.gz,test/sample1_L003_R2.fastq.gz,test/puck1.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample1,test/sample1_L004_R1.fastq.gz,test/sample1_L004_R2.fastq.gz,test/puck1.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample1,test/sample1_L005_R1.fastq.gz,test/sample1_L005_R2.fastq.gz,test/puck1.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample2,test/sample2_L001_R1.fastq.gz,test/sample2_L001_R2.fastq.gz,test/puck2.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample2,test/sample2_L002_R1.fastq.gz,test/sample2_L002_R2.fastq.gz,test/puck2.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample2,test/sample2_L003_R1.fastq.gz,test/sample2_L003_R2.fastq.gz,test/puck2.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample2,test/sample2_L004_R1.fastq.gz,test/sample2_L004_R2.fastq.gz,test/puck2.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
sample2,test/sample2_L005_R1.fastq.gz,test/sample2_L005_R2.fastq.gz,test/puck2.csv,8C18U7C2X8M,/camp/svc/reference/Genomics/babs/mus_musculus/ensembl/GRCm38/release-95/genome_idx/star/100bp,/camp/svc/reference/Genomics/babs/mus_musculus/
ensembl/GRCm38/release-95/gtf/Mus_musculus.GRCm38.95.gtf
```

![Oligos](oligos.png)

