----
## Exogene

A workflow for detecting viral integrations from both short read and long read sequencing data.

----
## usage

`docker pull zstephens/exogene:v15`

`docker run -it -v ${HOME}:${HOME} zstephens/exogene:v15`

### Create human + viral reference sequence:

(from inside the container)

`/home/init_ref.sh \ `  
`    -i /path/to/hg38.fa \ `  
`    -o /path/to/hg38_plus_viral.fa`  

#### Using custom viral references:

If you wish to use viral reference sequences different than what Exogene uses by default, you can use the `-v` input option to specify a fasta file of viral genomes. Exogene expects that contigs are named in the following format:

`>accession_id full_name` (space delimited)

For example: `>NC_009334 Human herpesvirus 4, complete genome.`  

Additionally, the viral fasta should be indexed using `bwa index`  

### Running Exogene-SR (with BAM input)

`/home/Exogene-SR.sh \ `  
`    -b input.bam \ `  
`    -r hg38_plus_viral.fa \ `  
`    -o outDir/`  

If custom viral sequences were used, the `-v` input option will be required.

### Running Exogene-SR (with FQ input)

`/home/Exogene-SR.sh \ `  
`    -f1 read1.fq.gz \ `  
`    -f2 read2.fq.gz \ `  
`    -r hg38_plus_viral.fa \ `  
`    -o outDir/`  

Input FASTQ files must be gzipped. Currently Exogene-SR does not support single-end reads. If custom viral sequences were used, the `-v` input option will be required.

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

`python /home/combine_reports.py \`  
`    -s Viral_Reads_Report.tsv \ `  
`    -l Viral_Junctions_LongReads.tsv \ `  
`    -o combined_report_outDir/ \ `  
`    -ms minimum_number_of_softclipped_reads_per_site [1] \ `  
`    -md minimum_number_of_discordant_pairs_per_site [5]`  

Either -s or -l must be specified (or both, for a combined report). Viral_Reads_Report.tsv is created in the output directory of Exogene-SR, Viral_Junctions_LongReads.tsv is created in the output directory of Exogene-LR.

### Test Data:

The Docker container contains a small quantity of test data which can be processed as follows:

`/home/Exogene-SR.sh \ `  
`    -f1 /home/test_data/SRR3104446_1.fq.gz \ `  
`    -f2 /home/test_data/SRR3104446_2.fq.gz \ `  
`    -r /path/to/hg38_and_viral.fa \ `  
`    -o /path/to/out_SR/ `  

`/home/Exogene-LR.sh \ `  
`    -f /home/test_data/a1el_ccs.fq.gz \ `  
`    -r /path/to/hg38_and_viral.fa \ `  
`    -m hifi \ `  
`    -o /path/to/out_LR/ `  

For the included hg38+viral reference, the bwa/pbmm2 alignment steps require ~32GB of memory.
