import os
import re
import gzip
import argparse

import numpy as np

import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def getColor(i,N,colormap='jet'):
    cm = mpl.get_cmap(colormap)
    cNorm  = colors.Normalize(vmin=0, vmax=N+1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    colorVal = scalarMap.to_rgba(i)
    return colorVal


HUMAN_CHR  = [str(n) for n in range(1,22+1)] + ['X', 'Y']
HUMAN_CHR += ['chr'+n for n in HUMAN_CHR]
HUMAN_CHR  = {n:True for n in HUMAN_CHR}


# returns reference span
REF_CHAR  = 'MX=D'
READ_CHAR = 'MX=I'
def parse_cigar(cigar):
    letters = re.split(r"\d+",cigar)[1:]
    numbers = [int(n) for n in re.findall(r"\d+",cigar)]
    startPos = 0
    if letters[0] == 'S':
        startPos = numbers[0]
    endClip = 0
    if len(letters) > 1 and letters[-1] == 'S':
        endClip = numbers[-1]
    adj  = 0
    radj = 0
    for i in range(len(letters)):
        if letters[i] in REF_CHAR:
            adj += numbers[i]
        if letters[i] in READ_CHAR:
            radj += numbers[i]
    return (startPos, adj, radj, endClip)


def exists_and_is_nonZero(fn):
    if os.path.isfile(fn):
        if os.path.getsize(fn) > 0:
            return True
    return False


def makedir(d):
    if not os.path.isdir(d):
        os.system('mkdir '+d)


def get_nearest_transcript(myChr, pos, bedDat, max_dist=20000):
    if myChr in bedDat:
        # lazy and slow, but it gets the job done!
        closest_dist = 99999999999
        closest_meta = ''
        for n in bedDat[myChr]:
            if pos >= n[0] and pos <= n[1]:
                closest_dist = 0
                closest_meta = [n[2], n[3]]
                break
            my_dist = min([abs(pos-n[0]), abs(pos-n[1])])
            if my_dist < closest_dist:
                closest_dist = my_dist
                closest_meta = [n[2], n[3]]
        if closest_dist <= max_dist:
            return (closest_dist, closest_meta)
    return None


"""//////////////////////////////////////////////////
////////////    PARSE INPUT ARGUMENTS    ////////////
//////////////////////////////////////////////////"""

parser = argparse.ArgumentParser(description='plot_viral_long_reads.py')
parser.add_argument('-s', type=str, required=True,  metavar='<str>', help="* input.sam")
parser.add_argument('-o', type=str, required=True,  metavar='<str>', help="* output_report.tsv")
parser.add_argument('-p', type=str, required=True,  metavar='<str>', help="* output_pics/")
parser.add_argument('-t', type=str, required=False, metavar='<str>', help="transcripts.bed", default='transcripts_hg38.bed.gz')
args = parser.parse_args()

OUT_REPORT = args.o
OUT_DIR    = args.p
if OUT_DIR[-1] != '/':
    OUT_DIR += '/'
makedir(OUT_DIR)

# track for finding nearest-gene
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
BED_DIR = SIM_PATH + 'resources/'
TRANSCRIPT_TRACK = {}
#
IN_TRANSCRIPTS_BED = args.t
if IN_TRANSCRIPTS_BED[-3:] == '.gz':
    f = gzip.open(BED_DIR + IN_TRANSCRIPTS_BED, 'rt')
else:
    f = open(BED_DIR + IN_TRANSCRIPTS_BED, 'r')
for line in f:
    splt = line.strip().split('\t')
    if splt[0] not in TRANSCRIPT_TRACK:
        TRANSCRIPT_TRACK[splt[0]] = []
    TRANSCRIPT_TRACK[splt[0]].append([int(splt[1]), int(splt[2]), splt[3], splt[4]])
f.close()

# [readpos_start, readpos_end, ref, pos_start, pos_end, orientation, mapq]
ALIGNMENTS_BY_RNAME = {}
READDAT_BY_RNAME = {}
READLEN_BY_RNAME = {}

f = open(args.s, 'r')
for line in f:
    splt = line.strip().split('\t')
    cigar = splt[5]
    flag  = int(splt[1])
    ref   = splt[2]
    pos   = int(splt[3])
    rdat  = splt[9]
    rnm   = splt[0]
    mapq  = int(splt[4])

    orientation = 'FWD'
    if flag & 16:
        orientation = 'REV'

    cigDat = parse_cigar(cigar)

    readPos1 = cigDat[0]
    readPos2 = cigDat[0] + cigDat[2]
    readLen  = cigDat[0] + cigDat[2] + cigDat[3]
    if orientation == 'REV':
        [readPos1, readPos2] = [readLen - readPos2, readLen - readPos1]

    if rnm not in ALIGNMENTS_BY_RNAME:
        ALIGNMENTS_BY_RNAME[rnm] = []
        READLEN_BY_RNAME[rnm]    = 0
    if orientation == 'FWD':
        pos1, pos2 = pos, pos + cigDat[1]
    elif orientation == 'REV':
        pos1, pos2 = pos + cigDat[1], pos
    ALIGNMENTS_BY_RNAME[rnm].append([readPos1, readPos2, ref, pos1, pos2, orientation, mapq])

    READDAT_BY_RNAME[rnm] = rdat
    READLEN_BY_RNAME[rnm] = max([READLEN_BY_RNAME[rnm], readLen])
f.close()

# get rid of pathological cases where grep "virus" found a single alignment where "virus"
# just so happened to be in the quality score string. WHAT ARE THE ODDS???
for k in sorted(ALIGNMENTS_BY_RNAME.keys()):
    if len(ALIGNMENTS_BY_RNAME) <= 1:
        del ALIGNMENTS_BY_RNAME[k]

#
#
#
f_out  = open(OUT_REPORT, 'w')
vseq_out = {}
DELIM_MAJOR = '='
DELIM_MINOR = '-'
nPlot  = 1
for k in sorted(ALIGNMENTS_BY_RNAME.keys()):
    #print(k)
    abns_k = sorted(ALIGNMENTS_BY_RNAME[k])
    for i in range(len(abns_k)):
        n = sorted(abns_k)[i]
        #print(n)
        if i >= 1:
            if n[2] in HUMAN_CHR and abns_k[i-1][2] not in HUMAN_CHR:
                qc1 = str(n[1]-n[0]) + ',' + str(n[6])
                qc2 = str(abns_k[i-1][1]-abns_k[i-1][0]) + ',' + str(abns_k[i-1][6])
                qc3 = str(n[0]-abns_k[i-1][1])
                outStr = k + '\t' + n[2] + ':' + str(n[3]) + '\t' + abns_k[i-1][2] + ':' + str(abns_k[i-1][4])
                f_out.write(outStr + '\t' + qc1 + '\t' + qc2 + '\t' + qc3 + '\n')
        if i < len(abns_k)-1:
            if n[2] in HUMAN_CHR and abns_k[i+1][2] not in HUMAN_CHR:
                qc1 = str(n[1]-n[0]) + ',' + str(n[6])
                qc2 = str(abns_k[i+1][1]-abns_k[i+1][0]) + ',' + str(abns_k[i+1][6])
                qc3 = str(abns_k[i+1][0]-n[1])
                outStr = k + '\t' + n[2] + ':' + str(n[4]) + '\t' + abns_k[i+1][2] + ':' + str(abns_k[i+1][3])
                f_out.write(outStr + '\t' + qc1 + '\t' + qc2 + '\t' + qc3 + '\n')
    #print('')

    viral_spans = []
    for i in range(len(abns_k)):
        n = abns_k[i]
        if n[2] not in HUMAN_CHR:
            if len(viral_spans) == 0:
                viral_spans.append([i])
            else:
                if i == viral_spans[-1][-1]+1:
                    viral_spans[-1].append(i)
                else:
                    viral_spans.append([i])

    #print(abns_k)
    #print(viral_spans)
    #print('')
    for n in viral_spans:

        (anchored_left, anchored_right)  = (False, False)
        (anchorseq_left, anchorseq_right) = ('', '')
        affected_genes = []
        if min(n) > 0:
            anchored_left = True
            # anchor seq format: "chr:pos:orr:gap:nearest_gene"
            a = abns_k[min(n)-1]
            out_nge = ''
            nearest_hit = get_nearest_transcript(a[2], a[4], TRANSCRIPT_TRACK)
            if nearest_hit != None:
                out_nge = nearest_hit[1][1]
                affected_genes.append(out_nge)
            anchorseq_left = a[2] + ':' + str(a[4]) + ':' + a[5] + ':' + str(abns_k[min(n)][0] - a[1]) + ':' + out_nge
        if max(n) < len(abns_k)-1:
            anchored_right = True
            # anchor seq format: "chr:pos:orr:gap:nearest_gene"
            a = abns_k[max(n)+1]
            out_nge = ''
            nearest_hit = get_nearest_transcript(a[2], a[3], TRANSCRIPT_TRACK)
            if nearest_hit != None:
                out_nge = nearest_hit[1][1]
                affected_genes.append(out_nge)
            anchorseq_right = a[2] + ':' + str(a[3]) + ':' + a[5] + ':' + str(a[0] - abns_k[max(n)][1]) + ':' + out_nge

        myAnchor = 'l'*anchored_left+'r'*anchored_right

        seq_name  = ''
        seq_name += k + DELIM_MAJOR                                             # read name
        seq_name += DELIM_MINOR.join([abns_k[m][2]+':'+abns_k[m][5] for m in n]) + DELIM_MAJOR  # viral refs this read spans
        #seq_name += DELIM_MINOR.join([str(m) for m in n]) + DELIM_MAJOR        # aln num of viral alignments (out of all alignments in the read)
        seq_name += myAnchor                                        # anchored on left, right, or both?
        if len(anchorseq_left):
            seq_name += DELIM_MAJOR + anchorseq_left
        if len(anchorseq_right):
            seq_name += DELIM_MAJOR + anchorseq_right

        seq_dat  = READDAT_BY_RNAME[k][abns_k[n[0]][0]:abns_k[n[-1]][1]]
        if len(seq_dat) >= 400:
            gk = ','.join(affected_genes)
            if gk not in vseq_out:
                vseq_out[gk] = []
            vseq_out[gk].append('>'+seq_name+'\n')
            vseq_out[gk].append(seq_dat+'\n')

    BN = -0.05 * READLEN_BY_RNAME[k]
    BP =  0.05 * READLEN_BY_RNAME[k]

    y_offsets = [0.7*BP]

    polygons = []
    p_alpha  = []
    p_color  = []
    p_text   = []
    for i in range(len(ALIGNMENTS_BY_RNAME[k])):
        n = sorted(ALIGNMENTS_BY_RNAME[k])[i]
        p_alpha.append(float(n[6]+15.)/(60.+15.))

        if 'virus' in n[2]:
            p_color.append('yellow')
        else:
            p_color.append('purple')

        delta_pointy     = 0.8*(n[1]-n[0])
        delta_pointy_rev = 0.2*(n[1]-n[0])

        if n[5] == 'FWD':
            polygons.append(Polygon(np.array([[n[0],BN], [n[0],BP], [n[0]+delta_pointy,BP], [n[1],0.], [n[0]+delta_pointy,BN]]), closed=True))
        else:
            polygons.append(Polygon(np.array([[n[0]+delta_pointy_rev,BN], [n[0],0.], [n[0]+delta_pointy_rev,BP], [n[1],BP], [n[1],BN]]), closed=True))

        p_text.append((n[0], y_offsets[-1], n[2]+' : '+str(n[3])+' - '+str(n[4])))
        y_offsets.append(y_offsets[-1] - 0.40*BP)

    mpl.rcParams.update({'font.size': 18, 'font.weight':'bold'})

    fig = mpl.figure(0,figsize=(12,6))
    # <plot stuff here>
    ax = mpl.gca()
    for i in range(len(polygons)):
        ax.add_collection(PatchCollection([polygons[i]], alpha=p_alpha[i], color=p_color[i]))

    for i in range(len(p_text)):
        mpl.text(p_text[i][0], p_text[i][1], p_text[i][2], ha='left', fontsize=9)

    mpl.axis([0, READLEN_BY_RNAME[k]+1, 4.0*BN, 4.0*BP])
    mpl.yticks([],[])
    mpl.grid(linestyle='--', alpha=0.5)
    mpl.xlabel('read coordinates')
    mpl.title(k)
    mpl.tight_layout()

    mpl.savefig(OUT_DIR+'read_'+str(nPlot)+'.png')
    nPlot += 1
    mpl.close(fig)
f_out.close()

f_out = open(OUT_DIR+'viral_sequences.fa', 'w')
for k in sorted(vseq_out.keys()):
    for n in vseq_out[k]:
        f_out.write(n)
f_out.close()
