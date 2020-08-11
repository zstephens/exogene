import sys
import os

#python append_viral_ins.py pbsv_ins_virus.sam duster.retain Viral_Junctions_LongReads.tsv Viral_Ins_LongReads.tsv [hifi/clr]

IN_SAM  = sys.argv[1]
IN_DUST = sys.argv[2]
IN_TSV  = sys.argv[3]
OUT_TSV = sys.argv[4]
MODE    = sys.argv[5]

BUFFER  = 100

if MODE == 'hifi':
	read_suffix = '/ccs'
elif MODE == 'clr':
	read_suffix = '/clr'
else:
	print 'Unknown mode.'
	exit(1)

existing_junctions = []
f = open(IN_TSV, 'r')
for line in f:
	splt = line.strip().split('\t')
	splt2 = splt[1].split(':')
	existing_junctions.append((splt2[0], int(splt2[1])))
f.close()

read_whitelist = {}
f = open(IN_DUST, 'r')
for line in f:
	read_whitelist[line.strip()] = True
f.close()

out_dat = []
f = open(IN_SAM, 'r')
for line in f:
	splt = line.strip().split('\t')
	flag = int(splt[1])
	ref  = splt[2]
	if not(flag&2048) and ref[:5] == 'virus':		# only use primary alignments to virus
		if splt[0] in read_whitelist:				# and we passed complexity filters
			splt2 = splt[0].split('_')
			r_ref = splt2[2]
			r_pos = int(splt2[3])
			any_hit = False
			for n in existing_junctions:
				if n[0] == r_ref and abs(r_pos - n[1]) <= BUFFER:
					any_hit = True
					break
			if any_hit == False:
				ref_len    = str(100)			# pbsv default threshold
				ref_mapq   = str(20)			# pbsv default threshold
				viral_len  = str(len(splt[9]))
				viral_mapq = splt[4]
				qc_gap     = str(0)
				out_dat.append(['ins_'+str(len(out_dat)+1)+read_suffix,
					            r_ref+':'+str(r_pos),
					            splt[2]+':'+splt[3],
					            ref_len+','+ref_mapq,
					            viral_len+','+viral_mapq,
					            qc_gap])
f.close()

f = open(OUT_TSV, 'w')
for n in out_dat:
	f.write('\t'.join(n)+'\n')
f.close()
