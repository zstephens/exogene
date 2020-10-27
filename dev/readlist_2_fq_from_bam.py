import sys
import gzip
import time

# samtools view input.bam | python readlist_2_fq_from_bam.py readlist out_r1.fq out_r2.fq outReadCount

RLIST   = sys.argv[1]
OUT_R1  = sys.argv[2]
OUT_R2  = sys.argv[3]
OUT_CNT = sys.argv[4]

def get_file_handle(fn, rw='r'):
	if fn[-6:].lower() == '.fq.gz' or fn[-9:].lower() == '.fastq.gz':
		return gzip.open(fn, rw)
	elif fn[-3:].lower() == '.fq' or fn[-6:].lower() == '.fastq':
		return open(fn, rw)
	else:
		print 'Input must be .fq, .fastq, .fq.gz or .fastq.gz'
		exit(1)

# return the reverse complement of a string
RC_DICT = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
def RC(s):
	return ''.join(RC_DICT[n] for n in s[::-1])

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

fo_1 = get_file_handle(OUT_R1, 'w')
fo_2 = get_file_handle(OUT_R2, 'w')

tt = time.time()
nWritten = 0
for line in input_stream:
	if len(line) and line[0] != '#':
		splt  = line.strip().split('\t')
		rnm   = splt[0]
		flag  = int(splt[1])

		if rnm[-2:] == '/1' or rnm[-2:] == '/2':
			rnm = rnm[:-2]
		if rnm not in rDict:
			continue

		if flag&2048:	# skip supplementary alignments
			continue
		
		if flag&16:
			rdat = RC(splt[9])
			qdat = splt[10][::-1]
		else:
			rdat = splt[9]
			qdat = splt[10]

		if flag&1:
			if flag&64:
				rDict[rnm][0] = [rdat, qdat]
			elif flag&128:
				rDict[rnm][1] = [rdat, qdat]

		if len(rDict[rnm][0]) and len(rDict[rnm][1]):
			fo_1.write('@' + rnm + '/1\n')
			fo_1.write(rDict[rnm][0][0] + '\n')
			fo_1.write('+\n')
			fo_1.write(rDict[rnm][0][1] + '\n')
			fo_2.write('@' + rnm + '/2\n')
			fo_2.write(rDict[rnm][1][0] + '\n')
			fo_2.write('+\n')
			fo_2.write(rDict[rnm][1][1] + '\n')
			del rDict[rnm]
			nWritten += 1
			if nWritten%1000000 == 0:
				print nWritten, 'read pairs written', '('+str(int(time.time()-tt))+' sec)', splt[2] + ':' + splt[3]
print nWritten, 'read pairs written', '('+str(int(time.time()-tt))+' sec)'

fo_2.close()
fo_1.close()

fo_0 = open(OUT_CNT, 'w')
fo_0.write(str(nWritten))
fo_0.close()
