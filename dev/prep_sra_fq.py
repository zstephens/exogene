import sys
import gzip

# python prep_sra_fq.py input.fq output.fq readnum

# example command for getting input fq from SRA:
#
# ./fastq-dump --split-3 SRR3105120

def read_fq_entry(fq_file):
	myName = fq_file.readline().strip()[1:]
	myName = myName.split(' ')[0]
	if not myName:
		return ('','','')
	myRDat = fq_file.readline().strip()
	skip   = fq_file.readline().strip()
	myQDat = fq_file.readline().strip()
	return (myRDat, myQDat, myName)

IN_FQ  = sys.argv[1]
OUT_FQ = sys.argv[2]
R_NUM  = sys.argv[3]

if R_NUM != '1' and R_NUM != '2':
	print 'Readnumber must be 1 or 2'
	exit(1)

if IN_FQ[-6:].lower() == '.fq.gz' or IN_FQ[-9:].lower() == '.fastq.gz':
	f_in  = gzip.open(IN_FQ,'r')
elif IN_FQ[-3:].lower() == '.fq' or IN_FQ[-6:].lower() == '.fastq':
	f_in  = open(IN_FQ,'r')
else:
	print 'Input must be .fq, .fastq, .fq.gz or .fastq.gz'
	exit(1)

if OUT_FQ[-6:].lower() == '.fq.gz' or OUT_FQ[-9:].lower() == '.fastq.gz':
	f_out = gzip.open(OUT_FQ,'w')
elif OUT_FQ[-3:].lower() == '.fq' or OUT_FQ[-6:].lower() == '.fastq':
	f_out = open(OUT_FQ,'w')
else:
	print 'Output must be .fq, .fastq, .fq.gz or .fastq.gz'
	exit(1)

while True:
	(myRDat, myQDat, myName) = read_fq_entry(f_in)
	if myRDat == '':
		break
	f_out.write('@' + myName + '/' + R_NUM + '\n')
	f_out.write(myRDat + '\n')
	f_out.write('+\n')
	f_out.write(myQDat + '\n')
f_out.close()
f_in.close()
