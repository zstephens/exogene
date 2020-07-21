#!/bin/bash

USAGE_1="Usage (bam): `basename $0` -b input.bam -r ref.fa -o output_dir/"
USAGE_2="Usage (fq):  `basename $0` -f1 read1.fq.gz -f2 read2.fq.gz -r ref.fa -o output_dir/"
USAGE_3="optional arguments: -k bwa_seed_size [30] -d duster_exclude_frac [70]"

#input argument parsing
ARG_BAM=""
ARG_R1=""
ARG_R2=""
ARG_REF=""
ARG_OUT=""
# bwa seed size
ARG_BWA_SEED=30
# toss out reads where >70% is flagged by duster as low complexity
ARG_DUSTER_FRAC=70
# read input args
while [[ "$#" -gt 0 ]]; do case $1 in
  -b|--bam) ARG_BAM="$2"; shift;;
  -f1|--read1) ARG_R1="$2"; shift;;
  -f2|--read2) ARG_R2="$2"; shift;;
  -r|--reference) ARG_REF="$2"; shift;;
  -o|--outdir) ARG_OUT="$2"; shift;;
  -k|--bwaseed) ARG_BWA_SEED="$2"; shift;;
  -d|--dusterfrac) ARG_DUSTER_FRAC="$2"; shift;;
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

# exes
samtools=/research/bsi/tools/biotools/samtools/0.1.18/samtools
bwa=/research/bsi/tools/biotools/bwa/0.7.17/bwa
perl=/research/bsi/tools/biotools/perl/5.16.3/bin/perl
duster=/research/bsi/tools/biotools/blast/2.6.0/bin/dustmasker
awk=/bin/awk
closestbed=/research/bsi/tools/biotools/bedtools/2.26.0/bin/closestBed
intersectbed=/research/bsi/tools/biotools/bedtools/2.26.0/bin/intersectBed
clusterbed=/research/bsi/tools/biotools/bedtools/2.26.0/bin/clusterBed
mergeBed=/research/bsi/tools/biotools/bedtools/2.26.0/bin/mergeBed
sortBed=/research/bsi/tools/biotools/bedtools/2.26.0/bin/sortBed

# scripts
duster_filt="/research/bsi/tools/biotools/smrtlink/5.1.0/smrtcmds/bin/python2.7 /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/duster_filter.py"
readlist_to_fq="/research/bsi/tools/biotools/smrtlink/5.1.0/smrtcmds/bin/python2.7 /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/readlist_2_fq.py"
readlist_to_fq_bam="/research/bsi/tools/biotools/smrtlink/5.1.0/smrtcmds/bin/python2.7 /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/readlist_2_fq_from_bam.py"
viralreads_to_report="/research/bsi/tools/biotools/smrtlink/5.1.0/smrtcmds/bin/python2.7 /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/createViralReadsReport.py"

# refs
HVR=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/HumanViral_Reference_12-12-2018.fa
RNA=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/GRCh38_RNA_Reference_4-30-2019.fa
Decoy=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/GRCh38_Decoy_Reference_4-30-2019.fa

# other resources
ExcludeRegions=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/Merged_ExcludeRegions.bed
Genes1KB=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/HG38_ProteinCoding-1KB.bed
GenesGS=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/HG38_ProteinCoding-GS.bed
GenesGSs=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/HG38_ProteinCoding-GS_sort.bed
ViralKey=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/HumanViral_Reference_12-12-2018.key
ViralAcc=/research/bsi/projects/staff_analysis/m112876/Viral_Workflow/references/HumanViral_Reference_12-12-2018.accessions

if [ "$INPUT_MODE" == "bam" ]; then
  # check input bam file path
  if [ ! -f $ARG_BAM ]; then
    echo "
  Please check the path to the input BAM file."
    exit 1
  fi
elif [ "$INPUT_MODE" == "fq" ]; then
  # check input fq file paths
  if [ ! -f $ARG_R1 ]; then
    echo "
  Please check the path to the input read_1 file."
    exit 1
  fi
  if [ ! -f $ARG_R2 ]; then
    echo "
  Please check the path to the input read_2 file."
    exit 1
  fi
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

# create result folder
name=exogeneSR
cd $ARG_OUT
date >> software.log

if [ "$INPUT_MODE" == "bam" ]; then
  echo "input = $ARG_BAM" > software.log
  echo "identifying potential viral read candidates">> software.log
  # run bwa against viral reference, don't align again unless we have to
  if [ ! -f viral_reads_se.reads ]; then
    $samtools view $ARG_BAM | $perl -lane 'print"\@$F[0]\n$F[9]\n+\n$F[10]";' | $bwa mem -k $ARG_BWA_SEED -t 4 $HVR - | egrep -v '^@' | $perl -lane 'if($F[2]ne"\*"){print"$_"};' > viral_reads_se.reads
  fi
elif [ "$INPUT_MODE" == "fq" ]; then
  echo "input = $ARG_R1 $ARG_R2" > software.log
  echo "identifying potential viral read candidates" >> software.log
  # run bwa against viral reference, don't align again unless we have to
  if [ ! -f viral_reads_se.reads ]; then
    zcat $ARG_R1 $ARG_R2 | $bwa mem -k $ARG_BWA_SEED -t 4 $HVR - | egrep -v '^@' | $perl -lane 'if($F[2]ne"\*"){print"$_"};' > viral_reads_se.reads
  fi
fi
cut -f6 viral_reads_se.reads | sed s/'[0-9]'//g > viral_reads_se.cigar
#paste viral_reads_se.cigar viral_reads_se.reads | $perl -lane 'if(($F[0]eq"M")||($F[0]eq"SM")||($F[9]eq"MS")){print">$F[1]\n$F[10]"};' > viral_reads_se.fa
paste viral_reads_se.cigar viral_reads_se.reads | $perl -lane 'print">$F[1]\n$F[10]";' > viral_reads_se.fa
fgrep ">" viral_reads_se.fa | cut -b2- | sort | uniq > viral_reads_se.ids
date >> software.log

# run blast duster
if [ ! -f viral_reads_se.ids.cleaned ]; then
  echo "evaluating viral reads for repeats" >> software.log
  $duster -in viral_reads_se.fa -outfmt fasta | $duster_filt $ARG_DUSTER_FRAC duster.out duster.retain duster.remove
  cat duster.remove viral_reads_se.ids | sort | uniq -u > viral_reads_se.ids.cleaned
  cat duster.retain viral_reads_se.ids.cleaned | sort | uniq > viral_reads_se.keep
  echo "repeats evaluated" >> software.log
  date >> software.log
fi

# create paired end fastq from the reliable viral reads
#
# TO SIMPLIFY MY LIFE: IF MODE IS FQ, READS MUST END IN /1 or /2
#
if [ "$INPUT_MODE" == "bam" ]; then
  if [ ! -f viral_1.fq ] || [ ! -f viral_2.fq ]; then
    $samtools view $ARG_BAM | $readlist_to_fq_bam viral_reads_se.keep viral_1.fq viral_2.fq totalReads.count
  fi
elif [ "$INPUT_MODE" == "fq" ]; then
  if [ ! -f viral_1.fq ] || [ ! -f viral_2.fq ]; then
    $readlist_to_fq viral_reads_se.keep $ARG_R1 $ARG_R2 viral_1.fq viral_2.fq
  fi
fi
echo "ran fastq conversion on reliable viral reads and their mates" >> software.log
date >> software.log

# run bwa against viral reference plus HG38 and removed alignments to viral/human sequence similarity regions and hard to align human regions
if [ ! -f ${name}_viral.bam ] || [ ! -f bwa.log ]; then
  $bwa mem -Y -k $ARG_BWA_SEED -t 4 $ARG_REF viral_1.fq viral_2.fq 2>bwa.log | $samtools view -bS - > Viral.bam
  $samtools sort Viral.bam Viral.sort
  $samtools index Viral.sort.bam
  # $samtools view Viral.sort.bam -L $ExcludeRegions | cut -f1 > bad.list-tmp
  touch bad.list-tmp
  $bwa mem -Y -k $ARG_BWA_SEED -t 4 $RNA viral_reads_se.fa | $perl -lane' print"$F[2]\t$F[0]\t$F[5]";' | egrep '^ENST' | cut -f2- | rev | cut -b2- | rev | $perl -lane 'if($F[1]>100){print"$F[0]"};' >> bad.list-tmp
  sort bad.list-tmp | uniq > bad.list
  # $samtools view Viral.sort.bam | cut -f1,6 |  $awk -v OFS="\t" '{if($2~/.*M.*S.*M.*/)print}' | cut -f1 >> bad.list
  $samtools view -h Viral.sort.bam | fgrep -v -w -f bad.list > ${name}_viral.sam
  $samtools view -Sb ${name}_viral.sam > ${name}_viral.bam
  $samtools index ${name}_viral.bam
  rm Viral.bam bad.lis* ${name}_viral.sam
  echo "evaluated viral integration sites into the human genome (hg38)" >> software.log
  date >> software.log
fi

# create gene disturbance report
$samtools view ${name}_viral.bam | $perl -lane 'print"$F[2]\t$_";' | egrep -v '^chr' | cut -f2 | sort | uniq > genes.virreads
$samtools view -H ${name}_viral.bam > genes.sam
$samtools view ${name}_viral.bam | $perl -lane 'print"$F[2]\t$_";' | egrep '^chr' | fgrep -f genes.virreads | cut -f2- >> genes.sam
$samtools view -Sb genes.sam > genes.bam
$intersectbed -loj -bed -abam genes.bam -b $Genes1KB > genes.bed
cut -f16 genes.bed | sort | uniq -c | $perl -lane 'print"$F[0]\t$F[1]";' | sort -k2,2 > genes.hits
printf "Gene\tSupporting_Reads_1KB\tChr\tGene_Start\tGene_End\n" > Gene_Disturbance_Report.tsv
join -1 2 -2 4 genes.hits $GenesGS | tr ' ' '\t' | sort -k2nr,2 >> Gene_Disturbance_Report.tsv
rm genes*
echo "evaluated viral integration support within gene regions" >> software.log
date >> software.log

# create viral presence summary report
$samtools view ${name}_viral.bam | cut -f1,3 | $perl -lane 'print"$F[1]\t$F[0]";' | egrep -v '^chr' | sort -k1,1 -k2,2 | uniq | cut -f1 | uniq -c | $perl -lane 'print"$F[1]\t$F[0]";' | sort -k1,1 > viralpresence.hits
printf "Accession\tSupporting_Reads\tVirus_Description\n" > Viral_Presence_Report.tsv
join -1 1 -2 1 viralpresence.hits $ViralKey | tr ' ' '\t' | tr '!' ' ' | sort -k2nr,2 >> Viral_Presence_Report.tsv
rm viralpresence*
echo "created viral presence summary report" >> software.log
date >> software.log

# create viral presence read detail report
$samtools view ${name}_viral.bam | cut -f1,3 > VReads_1.1
fgrep -f $ViralAcc VReads_1.1 > VReads_1.2
cut -f1 VReads_1.2 | sort | uniq > VReads_1.3
$samtools view ${name}_viral.bam | fgrep -f VReads_1.3 | cut -f1,2,3,4,5,6,7,8,10 > VReads_1
$viralreads_to_report VReads_1 $ViralKey Viral_Reads_Report.tsv
rm VReads_*
echo "created viral read details report" >> software.log
date >> software.log

# create gene integration report
#    creating the virpos tmp file
$samtools view ${name}_viral.bam | $perl -lane '$math=$F[7]+1;if(($F[2]!~/^chr/)&($F[6]=~/^chr/)){print"$F[6]\t$F[7]\t$math\t$F[2]"};' | sed s/'chr'//g | sed s/'^X'/'23'/g | egrep -v '^[a-z]|^[A-Z]|random' | sort -k1n,1 -k2n,2 | $perl -lane 'print"chr$_";' | sed s/'chr23'/'chrX'/g > tmp_1
$clusterbed -d 500 -i tmp_1 | $perl -lane 'print"$_\t$F[3]-$F[4]";' > tmp_2
cut -f6 tmp_2 | sort | uniq -c | sort -nr | $perl -lane 'print"$F[0]\t$F[1]";' > tmp_3
echo "" > tmp_4; sed -i '1d' tmp_4
cut -f2 tmp_3 | while read i; do fgrep -w $i tmp_2 | sort -R | head -1 >> tmp_4; done
paste tmp_3 tmp_4 | $perl -lane 'print"$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[0]";' | sort -k1,1 -k2n,2 > tmp_5
$closestbed -d -a tmp_5 -b $GenesGSs | cut -f1,2,4,5,9,10 > tmp_6
printf "chr\tintegration_site\tviral_accession\tsupporting_reads\tnearest_gene\tdistance\n" > tmp_7-virpos
sort -k4nr,4 -k5,5 tmp_6 >> tmp_7-virpos
rm tmp_4
#    creating the humpos tmp file
$samtools view ${name}_viral.bam | $perl -lane '$math=$F[3]+1;if(($F[2]=~/^chr/)&(($F[6]!~/^chr/)&($F[6]ne"="))){print"$F[2]\t$F[3]\t$math\t$F[6]"};' | sed s/'chr'//g | sed s/'^X'/'23'/g | egrep -v '^[a-z]|^[A-Z]|random' | sort -k1n,1 -k2n,2 | $perl -lane 'print"chr$_";' | sed s/'chr23'/'chrX'/g > tmp_1
$clusterbed -d 500 -i tmp_1 | $perl -lane 'print"$_\t$F[3]-$F[4]";' > tmp_2
cut -f6 tmp_2 | sort | uniq -c | sort -nr | $perl -lane 'print"$F[0]\t$F[1]";' > tmp_3
echo "" > tmp_4; sed -i '1d' tmp_4
cut -f2 tmp_3 | while read i; do fgrep -w $i tmp_2 | sort -R | head -1 >> tmp_4; done
paste tmp_3 tmp_4 | $perl -lane 'print"$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[0]";' | sort -k1,1 -k2n,2 > tmp_5
$closestbed -d -a tmp_5 -b $GenesGSs | cut -f1,2,4,5,9,10 > tmp_6
printf "chr\tintegration_site\tviral_accession\tsupporting_reads\tnearest_gene\tdistance\n" > tmp_7-humpos
sort -k4nr,4 -k5,5 tmp_6 >> tmp_7-humpos
#    creating the merged report
sed '1d' tmp_7-humpos | $perl -lane '$start=$F[1]-500;$end=$F[1]+500;print"$F[0]\t$start\t$end\t$F[2]\t$F[3]\t$F[4]\t$F[5]";' > tmp_8-humpos
sed '1d' tmp_7-virpos | $perl -lane '$start=$F[1]-500;$end=$F[1]+500;print"$F[0]\t$start\t$end\t$F[2]\t$F[3]\t$F[4]\t$F[5]";' > tmp_8-virpos
printf "Chr\tIntegration_Site\tViral_Accession\tHuman_Read_Support\tViral_Read_Support\tNearest_Gene\tDistance\n" > tmp_9-integrations
$intersectbed -wao -a tmp_8-humpos -b tmp_8-virpos | $perl -lane '$math1=($F[1]+500);$math2=($F[8]+500);$math3=(($math1+$math2)/2);$math4=(($F[6]+$F[13])/2);$math5=(abs($math2-$math1));if(($math5<500)&($F[3]eq$F[10])&($F[5]eq$F[12])){print"$F[0]\t$math3\t$F[3]\t$F[4]\t$F[11]\t$F[12]\t$math4"};' | sed s/'\.[0-9]'//g | sed s/'^chr'//g | sed s/'^X'/'23'/g | sort -k1n,1 -k2n,2 -k3,3 | $perl -lane 'print"chr$_";' | sed s/'chr23'/'chrX'/g >> tmp_9-integrations
echo "Virus_Description" > tmp_10
sed '1d' tmp_9-integrations | cut -f3 | while read i; do fgrep -w $i Viral_Presence_Report.tsv | head -1 | cut -f3-; done >> tmp_10
paste tmp_9-integrations tmp_10 > Gene_Integration_Full_Report.tsv
#    creating the truncated report
head -1 tmp_9-integrations > Gene_Integration_Truncated_Report.tsv
sed '1d' tmp_9-integrations | $perl -lane '$end=$F[1]+1;$math=$F[3]+$F[4];print"$F[0]\t$F[1]\t$end\t$F[2]\t$F[3]\t$F[4]\t$F[5]\t$F[6]\t$math";' > truc_1
$clusterbed -d 500 -i truc_1 | sort -k10n,10 -k9nr,9 -k8n,8 | uniq -f9 | cut -f1,2,4,5,6,7,8 >> Gene_Integration_Truncated_Report.tsv
echo "Virus_Description" > tmp_11
sed '1d' Gene_Integration_Truncated_Report.tsv | cut -f3 | while read i; do fgrep -w $i Viral_Presence_Report.tsv | head -1 | cut -f3-; done >> tmp_11
paste Gene_Integration_Truncated_Report.tsv tmp_11 > tmp_12
mv tmp_12 Gene_Integration_Truncated_Report.tsv
rm truc_1 tmp_*
echo "evaluated viral integration support genome wide and annotated nearby genes" >> software.log
date >> software.log

# create sample QC report
if [ "$INPUT_MODE" == "bam" ]; then
  treads=$(cat totalReads.count)
  rm totalReads.count
elif [ "$INPUT_MODE" == "fq" ]; then
  fq_lc=$(wc -l < $ARG_R1)
  treads=$(echo "${fq_lc}/2" | bc)
fi
vreads=$(cat viral_reads_se.ids | wc -l)
lcomplex=$(sort duster.remove | uniq | wc -l)
decoyc=$($bwa mem -Y -t 4 $Decoy viral_1.fq viral_2.fq | $perl -lane' print"$F[2]\t$F[0]";' | egrep '^chrUn' | cut -f2 | sort | uniq | wc -l)
decoy=$($perl -e '$math=('$decoyc'/'$vreads')*100;print"$math";'); 
hq=$(sed '1d' Viral_Reads_Report.tsv | cut -f1,2,8,9 | $perl -lane 'print"$F[1]\t$F[0]\n$F[3]\t$F[2]";' | egrep -v '^chr' | cut -f2 | sort | uniq | wc -l)
excludecount=$($samtools view Viral.sort.bam -L $ExcludeRegions | cut -f1 | sort | uniq | wc -l)
seqsim=$(echo "(${vreads}-(${lcomplex}+${hq}+${excludecount}))" | bc)
printf "TotalReads\tPaired_ViralReads\tViralDecoy_Percent\tLowComplexity_ViralReads\tExcludeRegion_ViralReads\tSeqSimilarity_ViralReads\tHighQuality_ViralReads\n" | tr '\t' '\n' > Sample_QC_Report-1
printf "$treads\t$vreads\t$decoy\t$lcomplex\t$excludecount\t$seqsim\t$hq\n" | tr '\t' '\n' > Sample_QC_Report-2
paste Sample_QC_Report-1 Sample_QC_Report-2 > Sample_QC_Report.tsv
rm Sample_QC_Report-*
echo "created a sample level QC report" >> software.log
date >> software.log

# cleaning up temp files
mkdir -p temp_files
mv viral* duster.* Viral* temp_files/
mv temp_files/Viral_Presence_Report.tsv ./
mv temp_files/Viral_Reads_Report.tsv ./
echo "software finished" >> software.log
date >> software.log
