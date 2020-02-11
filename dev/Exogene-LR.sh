#!/bin/bash

USAGE_1="Usage (bam): `basename $0` -f input.fq -r ref.fa -m mode -o output_dir/"

# references and tools
samtools=/opt/conda/envs/samtools/bin/samtools
pbmm2=/opt/conda/envs/pbmm2/bin/pbmm2
pbsv=/opt/conda/envs/pbmm2/bin/pbsv
viral_db_json=/home/resources/HumanViral_Reference_12-12-2018_simpleNames.json

vcf_to_fa="python /home/scripts/readlist_2_fq.py"
grep_virus="python /home/scripts/grep_virus_from_sam.py"

if [ "$ARG_MODE" == "hifi" ]; then
  pbmm2_preset="--preset CCS"
  pbsv_disc_preset="-q 10"
  pbsv_call_preset="--ccs"
elif [ "$ARG_MODE" == "clr" ]; then
  pbmm2_preset="--preset SUBREAD"
  pbsv_disc_preset="-q 10"
  pbsv_call_preset=""
fi

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

