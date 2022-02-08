#!/bin/bash

USAGE_1="Usage: `basename $0` -i input.fa -v viral.fa -o output.fa"

ARG_REF=""
ARG_OUT=""
ARG_HVR=""
# read input args
while [[ "$#" -gt 0 ]]; do case $1 in
  -i|--inref) ARG_REF="$2"; shift;;
  -v|--viral) ARG_HVR="$2"; shift;;
  -o|--outref) ARG_OUT="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done
if [ "$ARG_REF" == "" ]; then
  echo
  echo "-i input missing"
  echo $USAGE_1
  echo
  exit 1
fi
if [ "$ARG_OUT" == "" ]; then
  echo
  echo "-o input missing"
  echo $USAGE_1
  echo
  exit 1
fi
if [ "$ARG_HVR" == "" ]; then
  HVR=/home/refs/HumanViral_Reference_02-07-2022.fa
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

# scripts
make_json="python /home/exogene/dev/make_exogene_json.py"

# exes
bwa=/usr/bin/bwa
samtools=/opt/conda/envs/samtools/bin/samtools

$samtools faidx $ARG_REF
n_contigs=$(wc -l < ${ARG_REF}.fai)

cat $ARG_REF $HVR > $ARG_OUT
$samtools faidx $ARG_OUT
$make_json $HVR $n_contigs ${ARG_OUT}.exogene.json
$bwa index $ARG_OUT
