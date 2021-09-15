#!/bin/bash

USAGE_1="Usage (hifi): `basename $0` -f input.fq.gz -r ref.fa -m hifi -o output_dir/"
USAGE_2="Usage (cls):  `basename $0` -f input.fa.gz -r ref.fa -m clr -o output_dir/"
USAGE_3="Usage (bam):  `basename $0` -b input.bam -r ref.fa -m [hifi/clr] -o output_dir/"

#input argument parsing
ARG_READS=""
ARG_BAM=""
ARG_REF=""
ARG_MODE=""
ARG_OUT=""
# read input args
while [[ "$#" -gt 0 ]]; do case $1 in
  -f|--reads) ARG_READS="$2"; shift;;
  -b|--bam) ARG_BAM="$2"; shift;;
  -r|--reference) ARG_REF="$2"; shift;;
  -m|--readmode) ARG_MODE="$2"; shift;;
  -o|--outdir) ARG_OUT="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

# make sure reads/bam file was specified
if [ "$ARG_READS" == "" ] && [ "$ARG_BAM" == "" ]; then
  echo
  echo "must specify either -f or -b"
  echo $USAGE_1
  echo $USAGE_2
  echo $USAGE_3
  echo
  exit 1
fi
if [ "$ARG_READS" != "" ] && [ "$ARG_BAM" != "" ]; then
  echo
  echo "must specify either -f or -b"
  echo $USAGE_1
  echo $USAGE_2
  echo $USAGE_3
  echo
  exit 1
fi
# make sure mode was specified
if [ "$ARG_MODE" == "" ]; then
  echo
  echo "-m input missing"
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

# check input reads file path
if [ ! -f $ARG_READS ]; then
  echo "
Please check the path to the input reads file."
  exit 1
fi

# check input ref file path
if [ ! -f $ARG_REF ]; then
  echo "
Please check the path to the input ref file."
  exit 1
fi

# check specified output space path
mkdir -p $ARG_OUT
if [ ! -d $ARG_OUT ]; then
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
  echo $USAGE_2
  echo $USAGE_3
  echo
  exit 1
fi

# exes
samtools=/research/bsi/tools/biotools/samtools/1.10/bin/samtools
pbmm2=/research/bsi/tools/biotools/smrtlink/8.0/bin/smrtcmds/bin/pbmm2
pbsv=/research/bsi/tools/biotools/smrtlink/8.0/bin/smrtcmds/bin/pbsv

# scripts
grep_virus="/research/bsi/tools/biotools/smrtlink/8.0/bin/smrtcmds/bin/python /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/grep_virus_from_sam.py"
gen_report="/research/bsi/tools/biotools/smrtlink/8.0/bin/smrtcmds/bin/python /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/plot_viral_long_reads.py"
vcf_to_fa="/research/bsi/tools/biotools/smrtlink/8.0/bin/smrtcmds/bin/python /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/vcf_2_insfa.py"

# resources
viral_db_json=/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/resources/HumanViral_Reference_12-12-2018_simpleNames.json

mkdir -p $ARG_OUT
cd $ARG_OUT

# alignment
if [ ! -f pbmm2_viralReads.sam ]; then
  if [ "$ARG_BAM" == "" ]; then
    $pbmm2 align $ARG_REF $ARG_READS pbmm2_aln.bam $pbmm2_preset --sort --sample sample1 --rg '@RG\tID:movie1'
    $samtools view pbmm2_aln.bam | $grep_virus $viral_db_json > pbmm2_viralReads.sam
    MY_BAM="pbmm2_aln.bam"
  else
    $samtools view $ARG_BAM | $grep_virus $viral_db_json > pbmm2_viralReads.sam
    MY_BAM=$ARG_BAM
  fi
fi

# identify viral junctions
if [ ! -f Viral_Junctions_LongReads.tsv ]; then
  mkdir -p longread_plots
  $gen_report -s pbmm2_viralReads.sam -o Viral_Junctions_LongReads.tsv -p longread_plots/ -t t2t_v1.1_genes.bed.gz
fi

# SV calling
if [ ! -f temp.svsig.gz ]; then
  $pbsv discover $pbsv_disc_preset $MY_BAM temp.svsig.gz
fi
if [ ! -f pbsv_out.vcf ]; then
  $pbsv call $pbsv_call_preset -j 4 -t INS,DEL,INV,DUP,BND $ARG_REF temp.svsig.gz pbsv_out.vcf
fi

# extract large insertions to check for viral sequence
if [ ! -f pbsv_ins_virus.sam ]; then
  $vcf_to_fa -v pbsv_out.vcf -p PBSV -o pbsv_ins.fa
  $pbmm2 align $ARG_REF pbsv_ins.fa pbsv_ins.bam $pbmm2_preset --sort --sample sample1 --rg '@RG\tID:movie1'
  $samtools view pbsv_ins.bam | $grep_virus $viral_db_json > pbsv_ins_virus.sam
fi
