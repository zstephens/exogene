#!/bin/bash

USAGE_1="Usage (fq.gz): `basename $0` -f input.fq.gz -r ref.fa -m mode -o output_dir/"

#input argument parsing
ARG_READS=""
ARG_REF=""
ARG_MODE=""
ARG_OUT=""
# read input args
while [[ "$#" -gt 0 ]]; do case $1 in
  -f|--reads) ARG_READS="$2"; shift;;
  -r|--reference) ARG_REF="$2"; shift;;
  -m|--readmode) ARG_MODE="$2"; shift;;
  -o|--outdir) ARG_OUT="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# make sure reads file was specified
if [ "$ARG_READS" == "" ]; then
  echo
  echo "-f input missing"
  echo $USAGE_1
  echo
  exit 1
fi
# make sure mode was specified
if [ "$ARG_MODE" == "" ]; then
  echo
  echo "-m input missing"
  echo $USAGE_1
  echo
  exit 1
fi
# make sure a ref was specified
if [ "$ARG_REF" == "" ]; then
  echo
  echo "-r input missing"
  echo $USAGE_1
  echo
  exit 1
fi
# make sure an output dir was specified
if [ "$ARG_OUT" == "" ]; then
  echo
  echo "-o input missing"
  echo $USAGE_1
  echo
  exit 1
fi

# check input reads file path
RCK=$(ls -lt $ARG_READS | cut -d' ' -f5 | $perl -lane '$math=$F[0]*0;print"$math";')
if [ "$RCK" != "0" ]; then
  echo "
Please check the path to the input reads file."
  exit 1
fi

# check input ref file path
RCK=$(ls -lt $ARG_REF | cut -d' ' -f5 | $perl -lane '$math=$F[0]*0;print"$math";')
if [ "$RCK" != "0" ]; then
  echo "
Please check the path to the input ref file."
  exit 1
fi

# check specified output space path
mkdir -p $ARG_OUT
OSCK=$(ls -lt $ARG_OUT | cut -d' ' -f5 | $perl -lane '$math=$F[0]*0;print"$math";' | head -1)
if [ "$OSCK" != "0" ]; then
  echo "
Please check the path to the specified output directory."
  exit 1
fi

if [ "$ARG_MODE" == "hifi" ]; then
  pbmm2_preset="--preset CCS"
  pbsv_disc_preset="-q 10"
  pbsv_call_preset="--ccs"
elif [ "$ARG_MODE" == "clr" ]; then
  pbmm2_preset="--preset SUBREAD"
  pbsv_disc_preset="-q 10"
  pbsv_call_preset=""
else
  echo
  echo "-m must be either hifi or clr"
  echo $USAGE_1
  echo
  exit 1
fi

# references and tools
samtools=/opt/conda/envs/samtools/bin/samtools
pbmm2=/opt/conda/envs/pbmm2/bin/pbmm2
pbsv=/opt/conda/envs/pbmm2/bin/pbsv
viral_db_json=/home/resources/HumanViral_Reference_12-12-2018_simpleNames.json

vcf_to_fa="python /home/scripts/readlist_2_fq.py"
grep_virus="python /home/scripts/grep_virus_from_sam.py"

mkdir -p $ARG_OUT
cd $ARG_OUT

# alignment
$pbmm2 align $ARG_REF $ARG_READS pbmm2_aln.bam $pbmm2_preset --sort --sample sample1 --rg '@RG\tID:movie1'
$samtools view pbmm2_aln.bam | $grep_virus $viral_db_json > pbmm2_viralReads.sam

# SV calling
$pbsv discover $pbsv_disc_preset pbmm2_aln.bam temp.svsig.gz
$pbsv call $pbsv_call_preset -j 4 -t INS,DEL,INV,DUP,BND $ARG_REF temp.svsig.gz pbsv_out.vcf

exit 0

# extract large insertions to check for viral sequence
$vcf_to_fa pbsv_out.vcf pbsv_ins.fa
$pbmm2 align $ARG_REF pbsv_ins.fa pbsv_ins.bam $pbmm2_preset --sort --sample sample1 --rg '@RG\tID:movie1'
$samtools view pbsv_ins.bam | $grep_virus $viral_db_json > pbsv_ins_virus.sam
