import sys
import json

# dustmasker -in reads.fa -outfmt acclist | python duster_filter.py max_thresh duster.out duster.retain duster.remove

THRESH   = int(float(sys.argv[1]))
F_OUT    = sys.argv[2]
F_RETAIN = sys.argv[3]
F_REMOVE = sys.argv[4]

if not sys.stdin.isatty():
	input_stream = sys.stdin
else:
	print 'No input.'
	exit(1)

total_byReadName = {}
for line in input_stream:
	splt  = line.strip().split('\t')
	rname = splt[0][1:]
	sInd  = sorted([int(splt[1]), int(splt[2])])
	if rname not in total_byReadName:
		total_byReadName[rname] = 0
	total_byReadName[rname] += sInd[1] - sInd[0]

f1 = open(F_OUT,'w')
f2 = open(F_RETAIN,'w')
f3 = open(F_REMOVE,'w')
for k in sorted(total_byReadName.keys()):
	f1.write(str(total_byReadName[k])+'\t'+k+'\n')
	if total_byReadName[k] <= THRESH:
		f2.write(k+'\n')
	else:
		f3.write(k+'\n')
f3.close()
f2.close()
f1.close()
