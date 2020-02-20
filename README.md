----
## Exogene

A workflow for detecting viral integrations from both short read and long read sequencing data.

----
## usage

`docker pull exogene_v1`

`docker run -it -v ${HOME}:${HOME} exogene_v1`

### Create human + viral reference sequence:

(from inside the container)

`/home/scripts/init_ref.sh \ `  
`    -i /path/to/hg38.fa \ `  
`    -o /path/to/hg38_plus_viral.fa`  

### Running Exogene-SR (with BAM input)

`/home/Exogene-SR.sh \ `  
`    -b input.bam \ `  
`    -r hg38_plus_viral.fa \ `  
`    -o outDir/`  

### Running Exogene-SR (with FQ input)

`/home/Exogene-SR.sh \ `  
`    -f1 read1.fq.gz \ `  
`    -f2 read2.fq.gz \ `  
`    -r hg38_plus_viral.fa \ `  
`    -o outDir/`  

input FASTQ files must be gzipped, and readnames must end with either "/1" or "/2". Currently Exogene-SR does not support single-end reads.

### Running Exogene-LR (with FASTQ input, e.g. PacBio HiFi reads)

`/home/Exogene-LR.sh \ `  
`    -f input.fq.gz \ `  
`    -r hg38_plus_viral.fa \ `  
`    -m hifi \ `  
`    -o outDir/`  

### Running Exogene-LR (with FASTA input, e.g. PacBio CLR reads)

`/home/Exogene-LR.sh \ `  
`    -f input.fa.gz \ `  
`    -r hg38_plus_viral.fa \ `  
`    -m clr \ `  
`    -o outDir/`  

### Running Exogene-LR (with BAM input)

`/home/Exogene-LR.sh \ `  
`    -b input.bam \ `  
`    -r hg38_plus_viral.fa \ `  
`    -m [hifi/clr] \ `  
`    -o outDir/`  

### Intersecting Exogene-SR and Exogene-LR results:

`python /home/scripts/combine_reports.py \`  
`    -s Viral_Reads_Report.tsv \ `  
`    -l Viral_Junctions_LongReads.tsv \ `  
`    -o combined_report_outDir/ \ `  
`    -ms minimum_number_of_softclipped_reads_per_site [1] \ `  
`    -md minimum_number_of_discordant_pairs_per_site [5]`  

Either -s or -l must be specified (or both, for a combined report). Viral_Reads_Report.tsv is created in the output directory of Exogene-SR, Viral_Junctions_LongReads.tsv is created in the output directory of Exogene-LR.

### Output Report Description:

**Exogene-SR:**

* TBD

**Exogene-LR:**

* TBD

**Combined Reports:**

* TBD
