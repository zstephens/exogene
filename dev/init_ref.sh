#!/bin/bash

USAGE_1="Usage: `basename $0` -i input.fa -o output.fa"

ARG_REF=""
ARG_OUT=""
# read input args
while [[ "$#" -gt 0 ]]; do case $1 in
  -i|--inref) ARG_REF="$2"; shift;;
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

bwa=/opt/conda/envs/bwa/bin/bwa
HVR=/home/refs/HumanViral_Reference_12-12-2018.fa

cat $ARG_REF $HVR > $ARG_OUT
$bwa index $ARG_OUT
