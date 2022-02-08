import json
import sys

# python make_exogene_json.py viral.fa num-human-contigs out.json

VIRAL_FA  = sys.argv[1]
N_CONTIGS = int(sys.argv[2])
OUT_JSON  = sys.argv[3]

out_dict   = {}
virus_num  = 1
lines_read = 0

f = open(VIRAL_FA,'r')
for line in f:
	if line[0] == '>':
		out_dict['virus'+str(virus_num)] = line[1:-1]
		virus_num += 1
f.close()

out_dict['n_contigs'] = str(N_CONTIGS)

# write out read data dict
f = open(OUT_JSON,'w')
json.dump(out_dict,f,sort_keys=True,indent=4)
f.close()
