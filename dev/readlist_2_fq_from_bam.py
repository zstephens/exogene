import sys
import gzip

# samtools view input.bam | python readlist_2_fq_from_bam.py readlist out_r1.fq out_r2.fq

RLIST  = sys.argv[1]
OUT_R1 = sys.argv[2]
OUT_R2 = sys.argv[3]

def get_file_handle(fn, rw='r'):
	if fn[-6:].lower() == '.fq.gz' or fn[-9:].lower() == '.fastq.gz':
		return gzip.open(fn, rw)
	elif fn[-3:].lower() == '.fq' or fn[-6:].lower() == '.fastq':
		return open(fn, rw)
	else:
		print 'Input must be .fq, .fastq, .fq.gz or .fastq.gz'
		exit(1)

f = open(RLIST,'r')
rDict = {}
for line in f:
	#rDict[line.strip()] = ['','','','']
	rDict[line.strip()] = [[],[]]
f.close()

if not sys.stdin.isatty():
	input_stream = sys.stdin
else:
	print 'No input.'
	exit(1)

for line in input_stream:
	if len(line) and line[0] != '#':
		splt  = line.strip().split('\t')
		rnm   = splt[0]
		if rnm not in rDict:
			continue
		flag  = int(splt[1])
		rdat  = splt[9]
		qdat  = splt[10]

		if flag&1:
			if flag&64:
				rDict[rnm][0].append([rdat, qdat])
			elif flag&128:
				rDict[rnm][1].append([rdat, qdat])

# only include reads where we can find both members of a pair!
for k in sorted(rDict.keys()):
	if len(rDict[k][0]) and len(rDict[k][1]):
		rDict[k] = [rDict[k][0][0][0], rDict[k][0][0][1], rDict[k][1][0][0], rDict[k][1][0][1]]
	else:
		rDict[k] = ['','','','']

fo_1 = get_file_handle(OUT_R1, 'w')
fo_2 = get_file_handle(OUT_R2, 'w')
for k in sorted(rDict.keys()):
	if rDict[k][0] == '' or rDict[k][1] == '' or rDict[k][2] == '' or rDict[k][3] == '':
		continue
	fo_1.write('@' + k + '/1\n')
	fo_1.write(rDict[k][0] + '\n')
	fo_1.write('+\n')
	fo_1.write(rDict[k][1] + '\n')
	fo_2.write('@' + k + '/2\n')
	fo_2.write(rDict[k][2] + '\n')
	fo_2.write('+\n')
	fo_2.write(rDict[k][3] + '\n')
fo_2.close()
fo_1.close()

