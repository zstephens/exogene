import re
import sys

# python createViralReadsReport.py viralReadDat viralKeys outReport.tsv

HUMAN_CHR  = [str(n) for n in range(1,23)] + ['X', 'Y', 'M', 'Mt']
HUMAN_CHR += ['chr'+n for n in HUMAN_CHR]

# returns reference span
REF_CHAR = 'MX=D'
def parse_cigar_for_match(cigar):
	if cigar == '*':
		return 0
	letters = re.split(r"\d+",cigar)[1:]
	numbers = [int(n) for n in re.findall(r"\d+",cigar)]
	adj = 0
	for i in xrange(len(letters)):
		if letters[i] in REF_CHAR:
			adj += numbers[i]
	return adj

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
	myRefName = splt[2]
	if myRefName in viral_keyDict:
		myRefName = viral_keyDict[myRefName]
	if splt[0] not in outDat_by_rName:
		outDat_by_rName[splt[0]] = [[], []]
	myFlag = int(splt[1])
	if myFlag&1:
		if myFlag&64:
			outDat_by_rName[splt[0]][0].append([splt[0], splt[2], splt[3], splt[4], splt[5], splt[8], myRefName])
		elif myFlag&128:
			outDat_by_rName[splt[0]][1].append([splt[0], splt[2], splt[3], splt[4], splt[5], splt[8], myRefName])
f.close()

# R1_ID    R1_Contig    R1_Pos    R1_MAPQ    R1_CIGAR    R1_Seq    R1_RefName
f  = open(OUT_R, 'w')
f2 = open(OUT_R+'.multimapped', 'w')
f.write('R1_ID\tR1_Contig\tR1_Pos\tR1_MAPQ\tR1_CIGAR\tR1_Seq\tR1_RefName\t')
f.write('R2_ID\tR2_Contig\tR2_Pos\tR2_MAPQ\tR2_CIGAR\tR2_Seq\tR2_RefName\n')
for k in sorted(outDat_by_rName.keys()):
	if len(outDat_by_rName[k][0]) and len(outDat_by_rName[k][1]):

		# if both reads are multi-mapped...
		if len(outDat_by_rName[k][0]) >= 2 and len(outDat_by_rName[k][1]) >= 2:
			# lets go the extra mile to try and salvage some human/viral junctions, if we can
			isHuman0 = [n[1] in HUMAN_CHR for n in outDat_by_rName[k][0]]
			isHuman1 = [n[1] in HUMAN_CHR for n in outDat_by_rName[k][1]]
			if (True in isHuman0 and False in isHuman1) or (True in isHuman0 and False in isHuman1):
				allRefSpans = []
				for i in xrange(len(isHuman0)):
					if isHuman0[i]:
						allRefSpans.append((parse_cigar_for_match(outDat_by_rName[k][0][i][4]), 0, i))
				for i in xrange(len(isHuman1)):
					if isHuman1[i]:
						allRefSpans.append((parse_cigar_for_match(outDat_by_rName[k][1][i][4]), 1, i))
				allRefSpans = sorted(allRefSpans,reverse=True)
				outDat_by_rName[k][allRefSpans[0][1]] = [outDat_by_rName[k][allRefSpans[0][1]][allRefSpans[0][2]]]
			else:
				for n in outDat_by_rName[k][0]:
					f2.write('1\t' + '\t'.join(n) + '\n')
				for n in outDat_by_rName[k][1]:
					f2.write('2\t' + '\t'.join(n) + '\n')
				continue

		if len(outDat_by_rName[k][0]) == 1 and len(outDat_by_rName[k][1]) >= 2:
			indM = 1
		elif len(outDat_by_rName[k][0]) >= 2 and len(outDat_by_rName[k][1]) == 1:
			indM = 0
		else:
			indM = -1

		if indM >= 0:
			# if a read is multi-mapped, prioritize the following:
			# --- if mate is mapped to virus, choose human
			# --- if mate is mapped to human, choose virus
			# ------ this should cover 99.99% of cases, in the event of a tie, throw out the read
			goodMate  = abs(indM-1)
			mateHuman = outDat_by_rName[k][goodMate][0][1] in HUMAN_CHR
			candidates = []
			if mateHuman == False:
				candidates = [n for n in outDat_by_rName[k][indM] if n[1] in HUMAN_CHR]
			else:
				human_cig  = [n for n in outDat_by_rName[k][indM] if n[1] in HUMAN_CHR]
				print 'mate:', outDat_by_rName[k][goodMate]
				print human_cig
				candidates = [n for n in outDat_by_rName[k][indM] if not(n[1] in HUMAN_CHR)]
			if len(candidates) == 1:
				outDat_by_rName[k][indM] = [candidates[0]]
			else:
				for n in outDat_by_rName[k][0]:
					f2.write('1\t' + '\t'.join(n) + '\n')
				for n in outDat_by_rName[k][1]:
					f2.write('2\t' + '\t'.join(n) + '\n')
				continue

		outDat_by_rName[k] = [outDat_by_rName[k][0][0], outDat_by_rName[k][1][0]]
		f.write('\t'.join(outDat_by_rName[k][0]) + '\t' + '\t'.join(outDat_by_rName[k][1]) + '\n')
f2.close()
f.close()
