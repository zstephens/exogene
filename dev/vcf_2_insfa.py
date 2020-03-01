import os
import argparse

parser = argparse.ArgumentParser(description='vcf_2_insfa.py')
parser.add_argument('-v', type=str, required=True,  metavar='<str>', help="* input.vcf")
parser.add_argument('-p', type=str, required=True,  metavar='<str>', help="* program [pbsv/sniffles]")
parser.add_argument('-o', type=str, required=True,  metavar='<str>', help="* output.fa")
parser.add_argument('-l', type=int, required=False, metavar='<int>', help="minimum ins length [50]", default=50)
args = parser.parse_args()

IN_VCF    = args.v
PROGRAM   = args.p
OUT_READS = args.o
MIN_LEN   = args.l

f = open(IN_VCF, 'r')
INTERESTING_INSERTIONS = []
for line in f:
	if line[0] != '#':
		splt  = line.strip().split('\t')
		myDat = [splt[0], int(splt[1]), splt[3], splt[4]]
		if PROGRAM == 'PBSV' and 'INS' in  splt[2] and splt[4] != '<DUP>':
			INTERESTING_INSERTIONS.append([n for n in myDat])
		if PROGRAM == 'SNIFFLES' and 'SVTYPE=INS;' in splt[7]+';':
			INTERESTING_INSERTIONS.append([n for n in myDat])
f.close()

f = open(OUT_READS, 'w')
i = 0
for n in INTERESTING_INSERTIONS:
	if len(n[3]) >= MIN_LEN:
		f.write('>insertion_'+str(i)+'_'+str(n[0])+'_'+str(n[1])+'\n')
		f.write(n[3] + '\n')
		i += 1
f.close()
