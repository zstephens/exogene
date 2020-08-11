import sys
import os

#python append_viral_ins.py pbsv_ins_virus.sam Viral_Junctions_LongReads.tsv Viral_Ins_LongReads.tsv

IN_SAM  = sys.argv[1]
IN_TSV  = sys.argv[2]
OUT_TSV = sys.argv[3]

BUFFER  = 100

existing_junctions = []
f = open(IN_TSV, 'r')
for line in f:
	splt = line.strip().split('\t')
	splt2 = splt[1].split(':')
	existing_junctions.append((splt2[0], int(splt2[1])))
f.close()

f = open(IN_SAM, 'r')
for line in f:
	splt = line.strip().split('\t')
	flag = int(splt[1])
	ref  = splt[2]
	if not(flag&2048) and ref[:5] == 'virus':		# only use primary alignments to virus
		splt2 = splt[0].split('_')
		r_ref = splt2[2]
		r_pos = int(splt2[3])
		any_hit = False
		for n in existing_junctions:
			if n[0] == r_ref and abs(r_pos - n[1]) <= BUFFER:
				any_hit = True
				break
		if any_hit == False:
			print splt
f.close()
