#!/bin/bash

USAGE_1="Usage (bam): `basename $0` -b input.bam -r ref.fa -o output_dir/"
USAGE_2="Usage (fq):  `basename $0` -f1 read1.fq.gz -f2 read2.fq.gz -r ref.fa -o output_dir/"
USAGE_3="optional arguments: -v viral.fa -k bwa_seed_size [30] -d duster_exclude_frac [70] -t transcriptome_exclude_frac [90]"

#input argument parsing
ARG_BAM=""
ARG_R1=""
ARG_R2=""
ARG_REF=""
ARG_OUT=""
ARG_HVR=""
DELETE_TEMP="false"
# bwa seed size
ARG_BWA_SEED=30
# toss out reads where >70% is flagged by duster as low complexity
ARG_DUSTER_FRAC=70
# toss out reads where >90% aligns to transcriptome reference
ARG_TRANSCRIPT_FRAC=90
# read input args
while [[ "$#" -gt 0 ]]; do case $1 in
  -b|--bam) ARG_BAM="$2"; shift;;
  -f1|--read1) ARG_R1="$2"; shift;;
  -f2|--read2) ARG_R2="$2"; shift;;
  -r|--reference) ARG_REF="$2"; shift;;
  -o|--outdir) ARG_OUT="$2"; shift;;
  -v|--viral) ARG_HVR="$2"; shift;;
  -k|--bwaseed) ARG_BWA_SEED="$2"; shift;;
  -d|--dusterfrac) ARG_DUSTER_FRAC="$2"; shift;;
  -t|--transcriptfrac) ARG_TRANSCRIPT_FRAC="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done
# ensure we're properly in BAM or FQ mode
if [ "$ARG_BAM" != "" ] && [ "$ARG_R1" == "" ] && [ "$ARG_R2" == "" ]; then
  INPUT_MODE='bam'
elif [ "$ARG_BAM" == "" ] && [ "$ARG_R1" != "" ] && [ "$ARG_R2" != "" ]; then
  INPUT_MODE='fq'
else
  echo
  echo "Must specify either -b OR -f1 -f2 inputs"
  echo $USAGE_1
  echo $USAGE_2
  echo $USAGE_3
  echo
  exit 1
fi
# make sure a ref was specified
if [ "$ARG_REF" == "" ]; then
  echo
  echo "-r input missing"
  echo $USAGE_1
  echo $USAGE_2
  echo $USAGE_3
  echo
  exit 1
fi
# make sure an output dir was specified
if [ "$ARG_OUT" == "" ]; then
  echo
  echo "-o input missing"
  echo $USAGE_1
  echo $USAGE_2
  echo $USAGE_3
  echo
  exit 1
fi
# which viral references are we using?
if [ "$ARG_HVR" == "" ]; then
  HVR=/home/exogene/dev/resources/HumanViral_Reference_02-07-2022.fa
  echo
  echo "Using default viral reference sequences:"
  echo $HVR
  echo
else
  HVR=$ARG_HVR
  echo
  echo "Using user-specified viral reference sequences:"
  echo $HVR
  echo
fi

# exes
samtools=/opt/conda/envs/samtools/bin/samtools
bwa=/usr/bin/bwa
perl=/usr/bin/perl
duster=/usr/bin/dustmasker
awk=/usr/bin/awk
closestbed=/opt/conda/envs/bedtools/bin/closestBed
intersectbed=/opt/conda/envs/bedtools/bin/intersectBed
clusterbed=/opt/conda/envs/bedtools/bin/clusterBed
mergeBed=/opt/conda/envs/bedtools/bin/mergeBed
sortBed=/opt/conda/envs/bedtools/bin/sortBed
jq=/opt/conda/bin/jq

# scripts
python_path=python
exogene_path=/home/exogene/dev
duster_filt="$python_path $exogene_path/duster_filter.py"
readlist_to_fq="$python_path $exogene_path/readlist_2_fq.py"
readlist_to_fq_bam="$python_path $exogene_path/readlist_2_fq_from_bam.py"
aln_match_filter="$python_path $exogene_path/aln_match_filter.py"
viralreads_to_report="$python_path $exogene_path/createViralReadsReport.py"
combine_reports="$python_path $exogene_path/combine_reports.py"

# refs
RNA=/home/refs/GRCh38_RNA_Reference_4-30-2019.fa
Decoy=/home/refs/GRCh38_Decoy_Reference_4-30-2019.fa

# other resources
ExcludeRegions=/home/exogene/dev/resources/Merged_ExcludeRegions.bed
Genes1KB=/home/exogene/dev/resources/HG38_ProteinCoding-1KB.bed
GenesGS=/home/exogene/dev/resources/HG38_ProteinCoding-GS.bed
GenesGSs=/home/exogene/dev/resources/HG38_ProteinCoding-GS_sort.bed
viral_json=${ARG_REF}.exogene.json

if [ "$INPUT_MODE" == "bam" ]; then
  # check input bam file path
  if [ ! -f $ARG_BAM ]; then
    echo "
  Input BAM file not found:"
    echo $ARG_BAM
    exit 1
  fi
elif [ "$INPUT_MODE" == "fq" ]; then
  # check input fq file paths
  if [ ! -f $ARG_R1 ]; then
    echo "
  Input read1 file not found:"
    echo $ARG_R1
    exit 1
  fi
  if [ ! -f $ARG_R2 ]; then
    echo "
  Input read2 file not found:"
    echo $ARG_R2
    exit 1
  fi
fi

# check input ref file path
if [ ! -f $ARG_REF ]; then
  echo "
Reference fasta not found:"
  echo $ARG_REF
  exit 1
fi
# check input exogene json
if [ ! -f $viral_json ]; then
  echo "
exogene.json file not found:"
  echo $viral_json
  exit 1
fi
# check specified output space path
mkdir -p $ARG_OUT
if [ ! -d $ARG_OUT ]; then
  echo "
Please check the path to the specified output directory."
  exit 1
fi

# get number of human ref contigs from exogene json
n_contigs=$(${jq} ".n_contigs" ${viral_json} | tr -d \")

# check samtools version (for samtools sort input syntax)
sam_major=$($samtools 2>&1 | grep "Version:" | cut -d " " -f 2 | cut -d "." -f 1)
sam_minor=$($samtools 2>&1 | grep "Version:" | cut -d " " -f 2 | cut -d "." -f 2)
sam_minor=$(echo "$sam_minor < 3" | bc)
if [ "$sam_major" == "0" ]; then
  sam_is_new="0"
elif [ "$sam_minor" == "1" ]; then
  sam_is_new="0"
else
  sam_is_new="1"
fi

# create result folder
name=exogeneSR
cd $ARG_OUT
date >> software.log

#
# STEP 1) INITIAL SINGLE-END ALIGNMENT OF ALL READS TO VIRAL REFERENCE DATABASE
#
if [ "$INPUT_MODE" == "bam" ]; then
  echo "input = $ARG_BAM" > software.log
  echo "Step 1) SE alignment from BAM">> software.log
  # run bwa against viral reference, don't align again unless we have to
  if [ ! -f viral_reads_se.reads ]; then
    $samtools view -F 2048 $ARG_BAM | $perl -lane 'print"\@$F[0]\n$F[9]\n+\n$F[10]";' | $bwa mem -k $ARG_BWA_SEED -t 4 $HVR - | egrep -v '^@' | $perl -lane 'if($F[2]ne"\*"){print"$_"};' > viral_reads_se.reads
  fi
elif [ "$INPUT_MODE" == "fq" ]; then
  echo "input = $ARG_R1 $ARG_R2" > software.log
  echo "Step 1) SE alignment from FASTQ" >> software.log
  # run bwa against viral reference, don't align again unless we have to
  if [ ! -f viral_reads_se.reads ]; then
    zcat $ARG_R1 $ARG_R2 | $bwa mem -k $ARG_BWA_SEED -t 4 $HVR - | egrep -v '^@' | $perl -lane 'if($F[2]ne"\*"){print"$_"};' > viral_reads_se.reads
  fi
fi
date >> software.log

#
# STEP 2) REMOVE VIRAL-ALIGNED READS THAT ARE LOW COMPLEXITY
#
echo "Step 2) Low complexity filter" >> software.log
if [ ! -f viral_reads_se.fa ] || [ ! -f viral_reads_se.ids ]; then
  cut -f6 viral_reads_se.reads | sed s/'[0-9]'//g > viral_reads_se.cigar
  #paste viral_reads_se.cigar viral_reads_se.reads | $perl -lane 'if(($F[0]eq"M")||($F[0]eq"SM")||($F[9]eq"MS")){print">$F[1]\n$F[10]"};' > viral_reads_se.fa
  paste viral_reads_se.cigar viral_reads_se.reads | $perl -lane 'print">$F[1]\n$F[10]";' > viral_reads_se.fa
  fgrep ">" viral_reads_se.fa | cut -b2- | sort | uniq > viral_reads_se.ids
fi

if [ ! -f viral_reads_se.keep ]; then
  $duster -in viral_reads_se.fa -outfmt fasta | $duster_filt $ARG_DUSTER_FRAC duster.out duster.retain duster.remove
  cat duster.remove viral_reads_se.ids | sort | uniq -u > viral_reads_se.ids.cleaned
  cat duster.retain viral_reads_se.ids.cleaned | sort | uniq > viral_reads_se.keep
fi
date >> software.log

#
# STEP 3) RETRIEVE MATES OF ALL VIRAL-ALIGNED READS, CONSTRUCT PAIRED FASTQ
#
echo "Step 3) Retrieving mates to reconstruct PE FASTQ" >> software.log
if [ "$INPUT_MODE" == "bam" ]; then
  if [ ! -f viral_1.fq ] || [ ! -f viral_2.fq ]; then
    $samtools view $ARG_BAM | $readlist_to_fq_bam viral_reads_se.keep viral_1.fq viral_2.fq totalReads.count
  fi
elif [ "$INPUT_MODE" == "fq" ]; then
  if [ ! -f viral_1.fq ] || [ ! -f viral_2.fq ]; then
    $readlist_to_fq viral_reads_se.keep $ARG_R1 $ARG_R2 viral_1.fq viral_2.fq totalReads.count
  fi
fi
date >> software.log

#
# STEP 4) PAIRED-END ALIGNMENT TO HUMAN+VIRUS COMBINED REFERENCE
#
echo "Step 4) Aligning PE reads to human + viral reference" >> software.log
if [ ! -f ${name}_viral.bam ] || [ ! -f bwa.log ]; then
  $bwa mem -Y -k $ARG_BWA_SEED -t 4 $ARG_REF viral_1.fq viral_2.fq 2>bwa.log | $samtools view -bS - > Viral.bam
  if [ "$sam_is_new" == "1" ]; then
    $samtools sort -@ 3 -T Viral.temp -o Viral.sort.bam Viral.bam
  elif [ "$sam_is_new" == "0" ]; then
    $samtools sort Viral.bam Viral.sort
  fi
  $samtools index Viral.sort.bam
  # discard reads which align very well to transcriptome reference
  $bwa mem -Y -k $ARG_BWA_SEED -t 4 $RNA viral_reads_se.fa | $aln_match_filter $ARG_TRANSCRIPT_FRAC > rna_hits.ids
  # discard reads which align very well to decoy reference
  $bwa mem -Y -k $ARG_BWA_SEED -t 4 $Decoy viral_reads_se.fa | $aln_match_filter $ARG_TRANSCRIPT_FRAC > decoy_hits.ids
  sort rna_hits.ids decoy_hits.ids | uniq > bad.list
  $samtools view -h Viral.sort.bam | fgrep -v -w -f bad.list > ${name}_viral.sam
  $samtools view -Sb ${name}_viral.sam > ${name}_viral.bam
  $samtools index ${name}_viral.bam
  rm Viral.bam ${name}_viral.sam bad.list
fi
date >> software.log

#
# STEP 5) DETECT INTEGRATION BREAKPOINTS FROM SOFT-CLIP AND DISCORDANT EVIDENCE
#
echo "Step 5) Identifying integration sites" >> software.log
sed -e "1,${n_contigs}d" $ARG_REF.fai | cut -f1 > VReads_0.accessions
$samtools view ${name}_viral.bam | cut -f1,3 > VReads_1.1
fgrep -f VReads_0.accessions VReads_1.1 > VReads_1.2
cut -f1 VReads_1.2 | sort | uniq > VReads_1.3
$samtools view ${name}_viral.bam | fgrep -f VReads_1.3 > VReads_1
$viralreads_to_report VReads_1 $viral_json Viral_Reads_Report.tsv
rm VReads_*
# run combine_reports.py with some lenient thresholds to spit out all integrations
$combine_reports -s Viral_Reads_Report.tsv -o results/ -ms 2 -md 5 -v1 $viral_json
echo ""
cat results/integrations.tsv
echo ""
date >> software.log

#
# CLEAN UP TEMP FILES
#
if [ "$DELETE_TEMP" == "false" ]; then
  mkdir -p temp_files
  mv viral* duster.* Viral* *_hits.ids temp_files/
  #mv temp_files/Viral_Presence_Report.tsv ./
  mv temp_files/Viral_Reads_Report.tsv ./
elif [ "$DELETE_TEMP" == "true" ]; then
  rm *.ids*
  rm *.bam*
  rm *.fq
  rm *.count
  rm duster.*
  rm viral_reads_se.*
fi
echo "Done!" >> software.log
date >> software.log
