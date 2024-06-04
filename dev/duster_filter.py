import sys

# dustmasker -in reads.fa -outfmt fasta | python duster_filter.py max_lowcomplex_frac duster.out duster.retain duster.remove

THRESH   = float(sys.argv[1])/100.  # input is given as integer percent, e.g. "70"
F_OUT    = sys.argv[2]
F_RETAIN = sys.argv[3]
F_REMOVE = sys.argv[4]

if not sys.stdin.isatty():
    input_stream = sys.stdin
else:
    print('No input.')
    exit(1)

f1 = open(F_OUT,'w')
f2 = open(F_RETAIN,'w')
f3 = open(F_REMOVE,'w')
current_readname = ''
total_len        = 0
low_complex      = 0
for line in input_stream:
    if line[0] == '>':
        if current_readname != '':
            f1.write(str(current_readname)+'\t'+str(low_complex)+'\t'+str(total_len)+'\n')
            if float(low_complex)/float(total_len) <= THRESH:
                f2.write(current_readname+'\n')
            else:
                f3.write(current_readname+'\n')
        current_readname = line.strip()[1:]
        low_complex      = 0
        total_len        = 0
    else:
        rDat = line.strip()
        low_complex += rDat.count('a')
        low_complex += rDat.count('c')
        low_complex += rDat.count('g')
        low_complex += rDat.count('t')
        total_len   += len(rDat)

if total_len > 0:
    f1.write(str(current_readname)+'\t'+str(low_complex)+'\t'+str(total_len)+'\n')
    if float(low_complex)/float(total_len) <= THRESH:
        f2.write(current_readname+'\n')
    else:
        f3.write(current_readname+'\n')

f3.close()
f2.close()
f1.close()
