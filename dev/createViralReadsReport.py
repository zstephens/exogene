import sys

# python createViralReadsReport.py viralReadDat viralKeys outReport.tsv

V_DAT = sys.argv[1]
V_KEY = sys.argv[2]
OUT_R = sys.argv[3]

viral_keyDict = {}
f = open(V_KEY,'r')
for line in f:
	splt = line.strip().split('\t')
	viral_keyDict[splt[0]] = splt[1].replace('!',' ')
f.close()

outDat_by_rName = {}
f = open(V_DAT,'r')
for line in f:
	splt = line.strip('\n').split('\t')
	myRefName = splt[1]
	if myRefName in viral_keyDict:
		myRefName = viral_keyDict[myRefName]
	if splt[0] not in outDat_by_rName:
		outDat_by_rName[splt[0]] = []
	outDat_by_rName[splt[0]].append([splt[0], splt[1], splt[2], splt[3], splt[4], splt[7], myRefName])
f.close()

# R1_ID    R1_Contig    R1_Pos    R1_MAPQ    R1_CIGAR    R1_Seq    R1_RefName
f = open(OUT_R,'w')
f.write('R1_ID\tR1_Contig\tR1_Pos\tR1_MAPQ\tR1_CIGAR\tR1_Seq\tR1_RefName\t')
f.write('R2_ID\tR2_Contig\tR2_Pos\tR2_MAPQ\tR2_CIGAR\tR2_Seq\tR2_RefName\n')
for k in sorted(outDat_by_rName.keys()):
	if len(outDat_by_rName[k]) == 2:
		f.write('\t'.join(outDat_by_rName[k][0]) + '\t' + '\t'.join(outDat_by_rName[k][1]) + '\n')
f.close()
