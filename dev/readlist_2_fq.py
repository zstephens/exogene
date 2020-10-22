import sys
import gzip

# python readlist_2_fq.py readlist r1.fq r2.fq out_r1.fq out_r2.fq outReadCount

RLIST   = sys.argv[1]
IN_R1   = sys.argv[2]
IN_R2   = sys.argv[3]
OUT_R1  = sys.argv[4]
OUT_R2  = sys.argv[5]
OUT_CNT = sys.argv[6]

def get_file_handle(fn, rw='r'):
	if fn[-6:].lower() == '.fq.gz' or fn[-9:].lower() == '.fastq.gz':
		return gzip.open(fn, rw)
	elif fn[-3:].lower() == '.fq' or fn[-6:].lower() == '.fastq':
		return open(fn, rw)
	else:
		print 'Input must be .fq, .fastq, .fq.gz or .fastq.gz'
		exit(1)

def read_fq_entry(fq_file):
	myName = fq_file.readline().strip()[1:]
	myName = myName.split(' ')[0]
	if not myName:
		return ('','','')
	myRDat = fq_file.readline().strip()
	skip   = fq_file.readline().strip()
	myQDat = fq_file.readline().strip()
	return (myRDat, myQDat, myName)

f = open(RLIST,'r')
rDict = {}
for line in f:
	rDict[line.strip()] = ['','','','']
f.close()

fi_1 = get_file_handle(IN_R1, 'r')
fi_2 = get_file_handle(IN_R2, 'r')
while True:
	(myRDat, myQDat, myName) = read_fq_entry(fi_1)
	if myRDat == '':
		break
	if myName[-2:] == '/1' or myName[-2:] == '/2':
		myName = myName[:-2]
	if myName in rDict:
		rDict[myName][0] = myRDat
		rDict[myName][1] = myQDat
while True:
	(myRDat, myQDat, myName) = read_fq_entry(fi_2)
	if myRDat == '':
		break
	if myName[-2:] == '/1' or myName[-2:] == '/2':
		myName = myName[:-2]
	if myName in rDict:
		rDict[myName][2] = myRDat
		rDict[myName][3] = myQDat
fi_2.close()
fi_1.close()

fo_1 = get_file_handle(OUT_R1, 'w')
fo_2 = get_file_handle(OUT_R2, 'w')
nWritten = 0
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
	nWritten += 1
fo_2.close()
fo_1.close()

fo_0 = open(OUT_CNT, 'w')
fo_0.write(str(nWritten))
fo_0.close()
