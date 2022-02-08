import sys
import json

# samtools view input.bam | python grep_virus_from_sam.py viral_ref_names.json

if not sys.stdin.isatty():
	input_stream = sys.stdin
else:
	print 'No input.'
	exit(1)

# read in input readname:data dict
f = open(sys.argv[1],'r')
READNAME_DICT = json.load(f)
f.close()

REV_DICT      = {READNAME_DICT[k]:k for k in READNAME_DICT.keys() if k[:5] == 'virus'}
REV_DICT_ABRV = {READNAME_DICT[k].split(' ')[0]:k for k in READNAME_DICT.keys()}

ks1 = sorted(REV_DICT.keys())
ks2 = sorted(REV_DICT_ABRV.keys())

#
# replace all instances of annoyingly long viral ref names with shorthand
#
for line in input_stream:
	if len(line) and line[0] != '#' and line[0] != '@':
		splt = line.strip().split('\t')
		# replace refs in chr field
		if splt[2] in REV_DICT:
			splt[2] = REV_DICT[splt[2]]
		elif splt[2] in REV_DICT_ABRV:
			splt[2] = REV_DICT_ABRV[splt[2]]
		# replace refs in supplementary alignments
		for i in xrange(11,len(splt)):
			if splt[i][:3] == 'SA:' or splt[i][:3] == 'XA:':
				splt2 = [n.split(',') for n in splt[i].split(':')[-1].split(';') if len(n)]
				for j in xrange(len(splt2)):
					if splt2[j][0] in REV_DICT:
						splt2[j][0] = REV_DICT[splt2[j][0]]
					elif splt2[j][0] in REV_DICT_ABRV:
						splt2[j][0] = REV_DICT_ABRV[splt2[j][0]]
				splt[i] = splt[i][:5] + ';'.join([','.join(n) for n in splt2])
		outStr = '\t'.join(splt)+'\n'
		# if we got some virus, share it with the world!
		if 'virus' in outStr:
			sys.stdout.write(outStr)
