----
## Exogene

A workflow for detecting viral integrations from both short read and long read sequencing data.

----
## usage
`docker pull exogene_v1`
`docker run -it -v ${HOME}:${HOME} exogene_v1`

### Create human + viral reference sequence:

(from inside the container)
`/home/scripts/init_ref.sh -i /path/to/hg38.fa -o /path/to/hg38_plus_viral.fa`

### Running Exogene-SR (with BAM input)

`/home/Exogene-SR.sh -b input.bam -r hg38_plus_viral.fa -o outDir/`

### Running Exogene-SR (with FQ input)

`/home/Exogene-SR.sh -f1 read1.fq.gz -f2 read2.fq.gz -r hg38_plus_viral.fa -o outDir/`

input FASTQ files must be gzipped, and readnames must end with either "/1" or "/2". Currently Exogene-SR does not support single-end reads.

### Running Exogene-LR (with FASTQ input, e.g. PacBio HiFi reads)

TBD

### Running Exogene-LR (with FASTA input, e.g. PacBio CLR reads)

TBD

### Intersecting Exogene-SR and Exogene-LR results:

TBD

### Output Report Description:

**Exogene-SR:**

* TBD

**Exogene-LR:**

* TBD

**Combined Reports:**

* TBD
