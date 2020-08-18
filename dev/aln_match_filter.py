import re
import sys

# bwa mem -Y ref.fa reads.fa | python aln_match_filter.py thresh_percent > bad.list

thresh_percent = float(sys.argv[1])	# e.g. 90 = discard if >=90% of readlength is matched in a single alignment

if not sys.stdin.isatty():
	input_stream = sys.stdin
else:
	print 'No input.'
	exit(1)

bad_rname = {}

for line in input_stream:
	if len(line) and line[0] != '#' and line[0] != '@':
		splt  = line.strip().split('\t')
		cigar = splt[5]
		letters = re.split(r"\d+",cigar)[1:]
		numbers = [int(n) for n in re.findall(r"\d+",cigar)]
		readlen = len(splt[9])
		matchlen = sum([numbers[n] for n in xrange(len(numbers)) if letters[n] in 'M='])
		matchPercent = 100.*(float(matchlen)/float(readlen))
		if matchPercent >= thresh_percent:
			bad_rname[splt[0]] = True
		sys.stderr.write(' '.join([str(n) for n in [cigar, matchlen, readlen, matchPercent]])+'\n')

for k in sorted(bad_rname.keys()):
	sys.stdout.write(k + '\n')
