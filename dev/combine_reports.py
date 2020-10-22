#!/usr/bin/env python
# encoding: utf-8

from __future__ import print_function

import os
import re
import sys
import gzip
import json
import argparse
import numpy as np

import matplotlib.pyplot as mpl
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# absolute path to this script
SIM_PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
sys.path.append(SIM_PATH)

from mappability_corgi import MappabilityTrack
from sra_comparison import COMPARE

HUMAN_CHR  = [str(n) for n in range(1,22+1)] + ['X', 'Y']
HUMAN_CHR += ['chr'+n for n in HUMAN_CHR]
HUMAN_CHR  = {n:True for n in HUMAN_CHR}

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
                'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
                'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18, 'chrX' :23, 'chrY' :24}

TELOMERE_HG38 = [('chr1',0,10000), ('chr1',248946422,248956422), ('chr2',0,10000), ('chr2',242183529,242193529), ('chr3',0,10000), ('chr3',198285559,198295559), ('chr4',0,10000), ('chr4',190204555,190214555), ('chr5',0,10000), ('chr5',181528259,181538259), ('chr6',0,10000), ('chr6',170795979,170805979), ('chr7',0,10000), ('chr7',159335973,159345973), ('chr8',0,10000), ('chr8',145128636,145138636), ('chr9',0,10000), ('chr9',138384717,138394717), ('chrX',0,10000), ('chrX',156030895,156040895), ('chrY',0,10000), ('chrY',57217415,57227415), ('chr10',0,10000), ('chr10',133787422,133797422), ('chr11',0,10000), ('chr11',135076622,135086622), ('chr12',0,10000), ('chr12',133265309,133275309), ('chr13',0,10000), ('chr13',114354328,114364328), ('chr14',0,10000), ('chr14',107033718,107043718), ('chr15',0,10000), ('chr15',101981189,101991189), ('chr16',0,10000), ('chr16',90328345,90338345), ('chr17',0,10000), ('chr17',83247441,83257441), ('chr18',0,10000), ('chr18',80363285,80373285), ('chr19',0,10000), ('chr19',58607616,58617616), ('chr20',0,10000), ('chr20',64434167,64444167), ('chr21',0,10000), ('chr21',46699983,46709983), ('chr22',0,10000), ('chr22',50808468,50818468)]
TELOMERE_BUFF = 50000
CENTROMERE_HG38 = [('chr1',122503247,124785432), ('chr1',122026459,122224535), ('chr1',122224635,122503147), ('chr1',124849229,124932724), ('chr1',124785532,124849129), ('chr2',92188145,94090557), ('chr3',91553419,93655574), ('chr3',90772458,91233586), ('chr3',91233686,91247622), ('chr4',49712061,51743951), ('chr5',47153439,47296069), ('chr5',47309184,49591369), ('chr5',46485900,46569062), ('chr5',46569162,46796725), ('chr5',46796825,47061288), ('chr5',49667531,49721203), ('chr5',49721303,50059807), ('chr5',47106994,47153339), ('chr6',58553888,59829934), ('chr7',58169653,60828234), ('chr7',61377788,61528020), ('chr8',44033744,45877265), ('chr9',43389635,45518558), ('chrX',58605579,62412542), ('chrY',10316944,10544039), ('chr10',39686682,39935900), ('chr10',39936000,41497440), ('chr10',41545820,41593521), ('chr10',41497540,41545720), ('chr11',51090417,54342399), ('chr11',54342499,54425074), ('chr11',51078348,51090317), ('chr12',34835295,37185252), ('chr12',34769407,34816611), ('chr13',16282173,17416384), ('chr13',17418662,18051248), ('chr13',16110759,16164892), ('chr13',16249397,16256067), ('chr13',16000000,16022537), ('chr13',16022637,16110659), ('chr13',16164992,16228527), ('chr13',16228627,16249297), ('chr13',16256167,16259412), ('chr13',16259512,16282073), ('chr13',17416484,17416824), ('chr13',17416924,17417264), ('chr13',17417364,17418562), ('chr14',16404448,17538659), ('chr14',17540937,18173523), ('chr14',16228749,16282882), ('chr14',16377502,16400063), ('chr14',16000000,16022537), ('chr14',16140627,16228649), ('chr14',16282982,16346517), ('chr14',16346617,16367287), ('chr14',16367387,16374057), ('chr14',16374157,16377402), ('chr14',17538759,17539099), ('chr14',17539199,17539539), ('chr14',17539639,17540837), ('chr15',17499051,18355008), ('chr15',18355108,19725254), ('chr15',17083673,17498951), ('chr16',36337666,38265669), ('chr16',36311158,36334460), ('chr17',23195018,26566633), ('chr17',22813679,23194918), ('chr17',26566733,26616164), ('chr18',15797855,20561439), ('chr18',15460899,15780377), ('chr18',20696389,20736025), ('chr18',20839797,20861206), ('chr18',20603247,20696289), ('chr18',20736125,20813083), ('chr19',24908689,27190874), ('chr19',24498980,24552652), ('chr19',24552752,24891256), ('chr20',26608145,28494539), ('chr20',26436232,26586955), ('chr20',28648108,28728874), ('chr20',29917404,30038348), ('chr20',28508997,28556953), ('chr20',29125793,29204668), ('chr21',11146733,12280944), ('chr21',12283222,12915808), ('chr21',10864560,10887097), ('chr21',10975319,11029452), ('chr21',11124072,11146633), ('chr21',10887197,10975219), ('chr21',11029552,11093087), ('chr21',11093187,11113857), ('chr21',11113957,11120627), ('chr21',11120727,11123972), ('chr21',12281044,12281384), ('chr21',12281484,12281824), ('chr21',12281924,12283122), ('chr22',13285243,14419454), ('chr22',14421732,15054318), ('chr22',12954788,12977325), ('chr22',13021422,13109444), ('chr22',13227412,13248082), ('chr22',13109544,13163677), ('chr22',13163777,13227312), ('chr22',13248182,13254852), ('chr22',13254952,13258197), ('chr22',13258297,13280858), ('chr22',14419554,14419894), ('chr22',14419994,14420334), ('chr22',14420434,14421632)]
CENTROMERE_BUFF = 50000

#
#
#

def get_compare(c,p):
	#closest = None
	#minDist = -1
	#for n in COMPARE:
	#	if n[0] == c and abs(p-n[1]) < 1000:
	#		if minDist < 0:
	#			minDist = abs(p-n[1])
	#			closest = n
	#		elif abs(p-n[1]) < minDist:
	#			minDist = abs(p-n[1])
	#			closest = n
	#return closest
	outList = []
	for n in COMPARE:
		if n[0] == c and abs(p-n[1]) < 1000:
			outList.append(n)
	return outList

def isInBadRange(c,p):
	for n in TELOMERE_HG38:
		if n[0] == c and p > n[1]-TELOMERE_BUFF and p < n[2]+TELOMERE_BUFF:
			return True
	for n in CENTROMERE_HG38:
		if n[0] == c and p > n[1]-CENTROMERE_BUFF and p < n[2]+CENTROMERE_BUFF:
			return True
	return False

# returns (softclip_size, ref_position_clipping_begins)
REF_CHAR  = 'MX=D'
READ_CHAR = 'MX=I'
def parse_cigar_for_softclip(readDat):
	cigar = readDat[2]
	if 'S' not in cigar:
		return (0,None)
	letters = re.split(r"\d+",cigar)[1:]
	numbers = [int(n) for n in re.findall(r"\d+",cigar)]
	[left_clip, right_clip] = [(0,None), (0,None)]
	[adj, radj] = [0,0]
	for i in xrange(len(letters)):
		if i == 0 and letters[i] == 'S':
			left_clip = (numbers[i], readDat[1])
		elif i == len(letters)-1 and letters[i] == 'S':
			right_clip = (numbers[i], readDat[1] + adj)
		if letters[i] in REF_CHAR:
			adj += numbers[i]
		if letters[i] in READ_CHAR:
			radj += numbers[i]
	# if we have both left and right clipping, return the largest
	if left_clip[0] > right_clip[0]:
		return left_clip
	return right_clip

# returns reference span
def parse_cigar_for_match(readDat):
	cigar   = readDat[2]
	if cigar == '*':
		return 0
	letters = re.split(r"\d+",cigar)[1:]
	numbers = [int(n) for n in re.findall(r"\d+",cigar)]
	adj = 0
	for i in xrange(len(letters)):
		if letters[i] in REF_CHAR:
			adj += numbers[i]
	return adj

def close_enough(j1, j2, maxDist=100):
	# junction format: [chr_human, viral_name, human_pos_start, human_pos_end]
	if j1[0] != j2[0] or j1[1] != j2[1]:
		return False
	span_bounds = []
	span_bounds.append(abs(j1[2] - j2[2]))
	span_bounds.append(abs(j1[2] - j2[3]))
	span_bounds.append(abs(j1[3] - j2[2]))
	span_bounds.append(abs(j1[3] - j2[3]))
	if min(span_bounds) > maxDist:
		return False
	return True

def sorted_coord_list_to_ranges(l):
	if len(l) == 0:
		return []
	if len(l) == 1:
		return [(l[0],l[0])]
	ind = l[0]
	prev = l[0]
	outRanges = []
	for i in xrange(1,len(l)):
		if l[i] == prev+1:
			pass
		else:
			outRanges.append((ind,l[i-1]))
			ind = l[i]
		prev = l[i]
	if len(outRanges) == 0 or outRanges[-1][1] != l[-1]:
		outRanges.append((ind,l[-1]))
	return outRanges

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

def makedir(d):
	if not os.path.isdir(d):
		os.system('mkdir '+d)

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def rm(fn):
	if exists_and_is_nonZero(fn):
		os.system('rm '+fn)


"""//////////////////////////////////////////////////
////////////    PARSE INPUT ARGUMENTS    ////////////
//////////////////////////////////////////////////"""


parser = argparse.ArgumentParser(description='combine_reports.py')
parser.add_argument('-o',  type=str, required=True,  metavar='<str>', help="* output_dir/")
parser.add_argument('-s',  type=str, required=False, metavar='<str>', help="short_reads_report.tsv", default='')
parser.add_argument('-l',  type=str, required=False, metavar='<str>', help="long_reads_report.tsv", default='')
parser.add_argument('-v',  type=str, required=False, metavar='<str>', help="virus_of_interest", default='')
parser.add_argument('-v1', type=str, required=False, metavar='<str>', help="virusNames.json", default='')
parser.add_argument('-v2', type=str, required=False, metavar='<str>', help="virusAccessionToCommonName.nbr", default='')
parser.add_argument('-c',  type=str, required=False, metavar='<str>', help="SRR id to compare against", default='')
parser.add_argument('-a',  type=str, required=False, metavar='<str>', help="path/to/bed/annotations/", default='')
parser.add_argument('-b',  type=str, required=False, metavar='<str>', help="original_input.bam (for proximal softclip)", default='')
parser.add_argument('-ms', type=int, required=False, metavar='<int>', help="min number of SC reads per event", default=2)
parser.add_argument('-md', type=int, required=False, metavar='<int>', help="min number of disc pairs per event (if no sc)", default=5)
parser.add_argument('-ml', type=int, required=False, metavar='<int>', help="min number of long reads per event", default=0)
parser.add_argument('-mq', type=int, required=False, metavar='<int>', help="discard long read alns below this mapq", default=3)
parser.add_argument('-sn', type=str, required=False, metavar='<str>', help="sample name", default='')
parser.add_argument('--report-exclude', required=False, action='store_true', default=False, help='output bed-excluded integration sites')
args = parser.parse_args()

# basic parameters
(IN_SHORT, IN_LONG, OUT_DIR) = (args.s, args.l, args.o)

VOI = args.v
VOI = 'Hepatitis B virus'  ### FOR TESTING ONLY, REMOVE ME LATER

if IN_SHORT == '' and IN_LONG == '':
	print('Must specify either -s or -l')
	exit(1)

# look in same directory as output report for bwa.log
IN_BWALOG = None
if IN_SHORT != '':
	IN_BWALOG = IN_SHORT.split('/')
	IN_BWALOG[-1] = 'bwa.log'
	IN_BWALOG = '/'.join(IN_BWALOG)
	if exists_and_is_nonZero(IN_BWALOG) == False:
		print('Warning: bwa.log not found. Will be using default values.')

if OUT_DIR[-1] != '/':
	OUT_DIR += '/'
makedir(OUT_DIR)

COMP_SAMPLE    = args.c.upper()
COMPARE        = [(n.split('\t')[1], int(n.split('\t')[2]), n.split('\t')[3]) for n in COMPARE if n.split('\t')[0].upper() == COMP_SAMPLE]
COMPARE_OUT    = {n:[] for n in COMPARE}
COMPARE_OUT_FP = []

BED_DIR = args.a
if BED_DIR == '':
	BED_DIR = SIM_PATH + 'resources/'
BED_TRACKS = [['centromere',  MappabilityTrack(BED_DIR + 'hg38_centromere.bed',          bed_buffer=50000)],
              ['gap',         MappabilityTrack(BED_DIR + 'hg38_gap.bed',                 bed_buffer=1000)],
              ['repeats',     MappabilityTrack(BED_DIR + 'hg38_simpleRepeats.bed',       bed_buffer=50)],
              ['mappability', MappabilityTrack(BED_DIR + 'hg38_e2_l400_mappability.bed', bed_buffer=50)],
              ['exclude',     MappabilityTrack(BED_DIR + 'Merged_ExcludeRegions.bed',    bed_buffer=500)]]
              #['satellites',  MappabilityTrack(BED_DIR + 'hg38_microsatellites.bed',     bed_buffer=10)]]

# track for finding nearest-gene
TRANSCRIPT_TRACK = {}
f = gzip.open(BED_DIR + 'transcripts_hg38.bed.gz', 'r')
for line in f:
	splt = line.strip().split('\t')
	if splt[0] not in TRANSCRIPT_TRACK:
		TRANSCRIPT_TRACK[splt[0]] = []
	TRANSCRIPT_TRACK[splt[0]].append([int(splt[1]), int(splt[2]), splt[3], splt[4]])
f.close()

# read in viral short-hand that I used for the long read workflow
VIRAL_JSON = args.v1
if VIRAL_JSON == '':
	VIRAL_JSON = SIM_PATH + 'resources/HumanViral_Reference_12-12-2018_simpleNames.json'
if exists_and_is_nonZero(VIRAL_JSON) == False:
	print('Error: Viral simple-names json not found (-v1)')
	exit(1)
f = open(VIRAL_JSON, 'r')
VIRAL_NAME_DICT = json.load(f)
f.close()

# get viral accession ids and whatnot
VIRAL_ACCESSIONS = args.v2
if VIRAL_ACCESSIONS == '':
	VIRAL_ACCESSIONS = SIM_PATH + 'resources/taxid10239.nbr'
if exists_and_is_nonZero(VIRAL_ACCESSIONS) == False:
	print('Error: Viral accessions-to-name file not found (-v2)')
	exit(1)
f = open(VIRAL_ACCESSIONS, 'r')
ID_TO_ACCESSION       = {}
ACCESSION_TO_TAXONOMY = {}
for line in f:
	if line[0] != '#':
		splt = line.strip().split('\t')
		ID_TO_ACCESSION[splt[1]] = splt[0]
		ACCESSION_TO_TAXONOMY[splt[0]] = splt[4]
f.close()

SAMTOOLS    = '/opt/conda/envs/samtools/bin/samtools'
INPUT_BAM   = args.b
SAMPLE_NAME = args.sn

#
#
#

MIN_LONG_READ_MAPQ = args.mq	# skip a junction if both ends are below this mapq
PLOT_BUFF          = 100
MAX_LONG_READ_GAP  = 1000		# skip human --> viral junctions that have more than this much unexplained read sequence between them
MIN_SOFTCLIP       = args.ms
MIN_DISC_ONLY      = args.md	# if discordant reads are our only source of evidence, demand we have at least this many
MIN_SOFTCLIP_SIZE  = 10			# minimum softclipped size for reads that aren't anchored in virus (proximal softclips)
MIN_SOFTCLIP_MULT  = 2			# must have at least this many unanchored softclips at a particular coordinate
MIN_LONG_PER_EVENT = args.ml

REPORT_ALL_BEDHITS = args.report_exclude

if MIN_SOFTCLIP < 0:
	print('Error: -ms must be > 0')
	exit(1)
if MIN_DISC_ONLY < 0:
	print('Error: -md must be > 0')
	exit(1)
if MIN_LONG_PER_EVENT < 0:
	print('Error: -ml must be > 0')
	exit(1)
if MIN_SOFTCLIP <= 0 and MIN_DISC_ONLY <= 0:
	print('Error: either -ms or -md must be > 0')
	exit(1)

clustered_events = []
evidence_sc      = []	# softclip
evidence_pe      = []	# discordant paired-end reads
evidence_pb      = []	# pacbio

#
#	PARSE SHORT READS REPORT
#
tlen_mean = None
tlen_std  = None
if len(IN_SHORT):
	f = open(IN_SHORT,'r')
	isFirst = True
	data_byReadName = {}
	altc_byReadName = {}
	readLens = []
	for line in f:
		if isFirst:
			isFirst = False
			splt = line.strip().split('\t')
			ind_name  = splt.index('R1_ID')
			ind_ref   = splt.index('R1_Contig')
			ind_pos   = splt.index('R1_Pos')
			ind_cigar = splt.index('R1_CIGAR')
			ind_seq   = splt.index('R1_Seq')
			ind_mapq  = splt.index('R1_MAPQ')
			r2_offset = 7
			ind_altc  = None
			if 'alt_cigar' in splt:
				ind_altc = splt.index('alt_cigar')
			continue

		splt   = line.strip().split('\t')
		myName = splt[ind_name]
		myRef  = splt[ind_ref]
		myPos  = int(splt[ind_pos])
		myCig  = splt[ind_cigar]
		myMapQ = int(splt[ind_mapq])
		readLens.append(len(splt[ind_seq]))

		# read 1
		if myName not in data_byReadName:
			data_byReadName[myName] = []
		data_byReadName[myName].append([myRef, myPos, myCig, myMapQ])
		# read 2
		myRef  = splt[ind_ref+r2_offset]
		myPos  = int(splt[ind_pos+r2_offset])
		myCig  = splt[ind_cigar+r2_offset]
		myMapQ = int(splt[ind_mapq+r2_offset])
		data_byReadName[myName].append([myRef, myPos, myCig, myMapQ])
		# alternate cigar
		if ind_altc != None and splt[ind_altc] != '-':
			splt2 = splt[ind_altc].split('_')
			altc_byReadName[myName] = [splt2[0], int(splt2[1]), splt2[3], int(splt2[2])]

	f.close()

	# tlen stats (for use in estimating breakpoint from discordant paired-end evidence)
	if exists_and_is_nonZero(IN_BWALOG):
		f = open(IN_BWALOG, 'r')
		for line in f:
			if '[M::mem_pestat] mean and std.dev:' in line:
				splt = line.strip().split('(')[-1][:-1].split(',')
				tlen_mean = int(float(splt[0]) + 2.*np.mean(readLens))
				tlen_std  = int(float(splt[1]))
				break
		f.close()
	if tlen_mean == None:
		print('We were unable to grab template length stats from bwa.log, making a complete guess...')
		tlen_mean = 350
		tlen_std  = 50
	print('estimated template length:', tlen_mean, tlen_std)

	# cluster reads by event
	for k in data_byReadName.keys():
		# only read pairs where both passed filters...
		if len(data_byReadName[k]) == 2:
			# human - virus only
			if sum([1*(n[0] in HUMAN_CHR) for n in data_byReadName[k]]) == 1:

				# sort pairs such that read 1 is mapped to human
				if data_byReadName[k][0][0] in HUMAN_CHR:
					[r1, r2] = data_byReadName[k]
				else:
					[r2, r1] = data_byReadName[k]

				# collapse viral reference ids to common name
				if r2[0] in ID_TO_ACCESSION:
					r2[0] = ID_TO_ACCESSION[r2[0]]
				if r2[0] in ACCESSION_TO_TAXONOMY:
					r2[0] = ACCESSION_TO_TAXONOMY[r2[0]]
				#print('--',[r1,r2])

				# speed things up, only extract reads pertaining to virus-of-interest
				if len(VOI) and r2[0] != VOI:
					continue

				# softclip evidence
				sc1 = parse_cigar_for_softclip(r1)
				sc2 = parse_cigar_for_softclip(r2)
				scDat = None
				if sc1[0] > 0:
					scDat = [r1[0], sc1[1], sc1[0], r1[3]]
				# prefer alt cigar if present
				if k in altc_byReadName:
					if altc_byReadName[k][0] == r1[0] and abs(altc_byReadName[k][1]-r1[1]) <= tlen_mean+10*tlen_std:
						sc3 = parse_cigar_for_softclip(altc_byReadName[k])
						if sc3[0] > 0:
							scDat = [r1[0], sc3[1], sc3[0], max([r2[3], altc_byReadName[k][3]])]	# use mapping quality of viral alignment (suppl guaranteed to be 0?)

				# paired-end evidence
				span1 = parse_cigar_for_match(r1)
				span2 = parse_cigar_for_match(r2)
				bp_range1_left  = [max([0, r1[1] - tlen_mean - tlen_std + span1 + span2]), r1[1]]
				bp_range1_right = [r1[1] + span1, r1[1] + tlen_mean + tlen_std - span2]
				#print('--',[bp_range1_left, bp_range1_right])

				# closest-point hierarchal clustering
				myJunction = [r1[0], r2[0], r1[1], r1[1] + span1]
				found_a_hit = False
				for i in xrange(len(clustered_events)):
					for j in xrange(len(clustered_events[i])):
						found_a_hit = close_enough(myJunction, clustered_events[i][j])
						if found_a_hit:
							clustered_events[i].append([n for n in myJunction])
							if scDat == None:
								evidence_sc[i].append(None)
							else:
								evidence_sc[i].append(tuple([n for n in scDat]))
							evidence_pe[i].append(tuple(bp_range1_left + bp_range1_right + [r1[3]]))
							break
					if found_a_hit:
						break
				if not found_a_hit:
					clustered_events.append([[n for n in myJunction]])
					if scDat == None:
						evidence_sc.append([None])
					else:
						evidence_sc.append([tuple([n for n in scDat])])
					evidence_pe.append([tuple(bp_range1_left + bp_range1_right + [r1[3]])])
					evidence_pb.append([])
				#print('SUM:', len(clustered_events), sum([len(n) for n in clustered_events]))

#
#	PARSE LONG READ EVENTS
#
if len(IN_LONG):
	f = open(IN_LONG, 'r')
	for line in f:
		splt = line.strip().split('\t')
		bp1  = splt[1].split(':')
		bp2  = splt[2].split(':')
		qc1  = splt[3].split(',')
		qc2  = splt[4].split(',')
		qc3  = splt[5]
	
		if '/ccs' in splt[0]:
			myType = 'CCS'
		else:
			myType = 'CLR'
	
		viral_ref = VIRAL_NAME_DICT[bp2[0]].split(' ')[0]
		if viral_ref in ID_TO_ACCESSION:
			viral_ref = ID_TO_ACCESSION[viral_ref]
		if viral_ref in ACCESSION_TO_TAXONOMY:
			viral_ref = ACCESSION_TO_TAXONOMY[viral_ref]
	
		if int(qc1[1]) < MIN_LONG_READ_MAPQ and int(qc2[1]) < MIN_LONG_READ_MAPQ:
			continue
	
		if int(qc3) > MAX_LONG_READ_GAP:
			continue
	
		myJunction = [bp1[0], viral_ref, int(bp1[1]), int(bp1[1])]
	
		# I'm copy-pasting this code and I don't even care!
		found_a_hit = False
		for i in xrange(len(clustered_events)):
			for j in xrange(len(clustered_events[i])):
				found_a_hit = close_enough(myJunction, clustered_events[i][j])
				if found_a_hit:
					clustered_events[i].append([n for n in myJunction])
					evidence_pb[i].append((int(bp1[1]), myType))
					break
			if found_a_hit:
				break
		if not found_a_hit:
			clustered_events.append([[n for n in myJunction]])
			evidence_pb.append([(int(bp1[1]), myType)])
	f.close()

# sort!
order_to_process_clusters = []
for i in xrange(len(clustered_events)):
	n = clustered_events[i]
	r1, r2 = n[0][0], n[0][1]
	pmin = min([m[2] for m in n])
	pmax = max([m[3] for m in n])
	if isInBadRange(r1,pmin) == False and isInBadRange(r1,pmax) == False:
		order_to_process_clusters.append((LEXICO_2_IND[r1],r2,pmin,pmax,i))
order_to_process_clusters = [n[4] for n in sorted(order_to_process_clusters)]

#
# if we're hunting for proximal softclips, grab all the relevant reads
#
proximal_softclips = {}
if len(INPUT_BAM) and exists_and_is_nonZero(INPUT_BAM) and len(IN_SHORT):
	sc_regions = []
	for i in order_to_process_clusters:
		if len(VOI) and clustered_events[i][0][1] != VOI:
			continue
		(min_coord, max_coord) = (999999999999,-1)
		if i < len(evidence_sc):
			for n in evidence_sc[i]:
				if n != None:
					if n[1] > max_coord: max_coord = n[1]
					if n[1] < min_coord: min_coord = n[1]
			for n in evidence_pe[i]:
				min_coord = min([min_coord, n[0], n[1], n[2], n[3]])
				max_coord = max([max_coord, n[0], n[1], n[2], n[3]])
		sc_regions.append([clustered_events[i][0][0], max([0, min_coord-tlen_mean]), max_coord+tlen_mean])
	sc_f = open(OUT_DIR+'sc_temp.bed', 'w')
	for n in sc_regions:
		sc_f.write(str(n[0])+'\t'+str(n[1])+'\t'+str(n[2])+'\n')
	sc_f.close()
	sc_track = MappabilityTrack(OUT_DIR+'sc_temp.bed', read_pickle=False)
	rm(OUT_DIR+'sc_temp.sam')
	os.system('touch ' + OUT_DIR+'sc_temp.sam')
	for k in sc_track.all_tracks.keys():
		for i in xrange(1, len(sc_track.all_tracks[k])-2, 2):
			myLoc = k+':'+str(sc_track.all_tracks[k][i:i+2][0])+'-'+str(sc_track.all_tracks[k][i:i+2][1])
			os.system(SAMTOOLS + ' view ' + INPUT_BAM + ' ' + myLoc + ' >> ' + OUT_DIR+'sc_temp.sam')
	sc_f = open(OUT_DIR+'sc_temp.sam', 'r')
	for line in sc_f:
		splt   = line.strip().split('\t')
		rname  = splt[0]
		myRef  = splt[2]
		myPos  = int(splt[3])
		myMapQ = int(splt[4])
		myCig  = splt[5]
		if rname not in data_byReadName and 'S' in myCig:
			myR = [myRef, myPos, myCig, myMapQ]
			myS = parse_cigar_for_softclip(myR)
			if myS[0] >= MIN_SOFTCLIP_SIZE:
				#print('rawr:', myRef, myPos, myCig, myS)
				myK = (myRef, myS[1])
				if myK not in proximal_softclips:
					proximal_softclips[myK] = 0
				proximal_softclips[myK] += 1
	sc_f.close()
	for k in proximal_softclips.keys():
		if proximal_softclips[k] < MIN_SOFTCLIP_MULT:
			del proximal_softclips[k]
	rm(OUT_DIR+'sc_temp.bed')
	rm(OUT_DIR+'sc_temp.sam')
#for k in sorted(proximal_softclips.keys()):
#	print(k, proximal_softclips[k])
#exit(1)

zfill_num = len(str(len(clustered_events)))

POLY_PE_ALPHA = [0.1, 0.2, 0.3, 0.5, 0.8]
POLY_PE_STEPS = len(POLY_PE_ALPHA)
POLY_PE_ALPHA = [0.0] + POLY_PE_ALPHA

#
#	LET'S TAKE A LOOK
#
nPlot = 0
class_counts = {}
bp_dist_sc_ccs = []
bp_dist_sc_clr = []
bp_dist_min_sc_ccs = []
bp_dist_min_sc_clr = []
bp_dev_sc  = []
bp_dev_ccs = []
bp_dev_clr = []
sc_count_by_class = {}
pe_count_by_class = {}
out_report_data = []
for i in order_to_process_clusters:

	#print(clustered_events[i][0][0] + ' --> ' + clustered_events[i][0][1])

	if len(VOI) and clustered_events[i][0][1] != VOI:
		continue

	sc_str  = ''
	ccs_str = ''
	clr_str = ''
	
	# collapse all sc/pe evidence down by choosing the mode
	sc_count = {}
	sc_mapq_byPos = {}
	pe_count = {}
	pe_mapq_byPos = {}
	pe_to_report_stratified = []
	if i < len(evidence_sc):
		for n in evidence_sc[i]:
			if n != None:
				if n[1] not in sc_count:
					sc_count[n[1]] = 0
					sc_mapq_byPos[n[1]] = []
				sc_count[n[1]] += 1
				sc_mapq_byPos[n[1]].append(n[3])
		for n in evidence_pe[i]:
			for j in range(n[0],n[1]+1)+range(n[2],n[3]+1):
				if j not in pe_count:
					pe_count[j] = 0
					pe_mapq_byPos[j] = []
				pe_count[j] += 1
				pe_mapq_byPos[j].append(n[4])

	mySCCount = 0
	myPECount = 0
	if not len(sc_count):
		sc_to_report = None
		max_sc       = 0
		sc_str       = ''

	else:
		sc_coord_list = [[k]*sc_count[k] for k in sc_count.keys()]
		sc_coord_list = [item for sublist in sc_coord_list for item in sublist]
		if len(sc_coord_list) > 1:
			scm = np.median(sc_coord_list)
			scl = [int(abs(n-scm)+0.5) for n in sc_coord_list]
			sc_str = 'type: shortread count: ' + str(len(scl)) + ' avg_pos: ' + str(int(scm+0.5)) + ' std: {0:.2f}'.format(np.mean(scl))
			bp_dev_sc.extend(scl)
			scm = int(scm+0.5)
		elif len(sc_coord_list) == 1:
			sc_str = 'type: shortread count: ' + str(1) + ' avg_pos: ' + str(sc_coord_list[0]) + ' std: 0.00'
			scm = sc_coord_list[0]
			scl = [-1]

		scm_closest   = sorted([(abs(scm-k_sc), k_sc) for k_sc in sc_mapq_byPos.keys()])[0][1]
		mapq0_frac    = float(sc_mapq_byPos[scm_closest].count(0))/len(sc_mapq_byPos[scm_closest])
		mapq0_percent = '{0:0.2f}%'.format(100.*mapq0_frac)
		sc_str += ' MAPQ0: ' + mapq0_percent

		if (len(args.c) or len(COMPARE)) and len(scl) >= MIN_SOFTCLIP:
			#sc_str += ' closest: ' + str(get_compare(clustered_events[i][0][0], int(scm+0.5)))
			sc_coords_to_compare = sorted([k for k in sc_count.keys() if sc_count[k] >= MIN_SOFTCLIP])
			sc_fps = []
			sc_anyHits = False
			for scc_to_use in sc_coords_to_compare:
				gcl = get_compare(clustered_events[i][0][0], scc_to_use)
				if len(gcl) > 0:
					sc_anyHits = True
					for gc in gcl:
						COMPARE_OUT[gc].append((clustered_events[i][0][0], scc_to_use, str(sc_count[scc_to_use])+'/'+str(len(scl)), mapq0_percent, abs(int(scc_to_use+0.5)-gc[1]), 'SOFTCLIP'))
				else:
					sc_fps.append((clustered_events[i][0][0], scc_to_use, str(sc_count[scc_to_use])+'/'+str(len(scl)), mapq0_percent))
			if sc_anyHits == False and len(sc_fps):
				#COMPARE_OUT_FP.append([n for n in sc_fps])
				COMPARE_OUT_FP.append((clustered_events[i][0][0], int(scm+0.5), len(scl), mapq0_percent))
		mySCCount = len(sc_coord_list)
		max_sc = max(sc_count.values())

		# if 1 read, just report it, if multiple, report all coords with >1 read
		if max_sc == 1:
			sc_to_report = sorted([k for k in sc_count.keys()])
			sc_to_report = sorted_coord_list_to_ranges(sc_to_report)
		else:
			sc_to_report = sorted([k for k in sc_count.keys() if sc_count[k] > 1])
			sc_to_report = sorted_coord_list_to_ranges(sc_to_report)

	if not len(pe_count):
		pe_to_report = None
		max_pe       = 0
		pe_str       = ''

	else:
		max_pe = max(pe_count.values())
		threshes = [n*(max_pe/POLY_PE_STEPS) for n in xrange(POLY_PE_STEPS)] + [max_pe-1]
		myPECount = max_pe
		#print('DISC_THRESHES:',threshes)
		#pe_to_report = sorted([k for k in pe_count.keys() if pe_count[k] == max_pe])
		pe_to_report = sorted([k for k in pe_count.keys()])
		pe_density   = [pe_count[k] for k in pe_to_report]
		for j in xrange(len(pe_density)):
			myPEDind = 0
			while pe_density[j] > threshes[myPEDind]:
				myPEDind += 1
				if myPEDind == len(threshes):
					break
			pe_density[j] = myPEDind-1
		pe_to_report_stratified = [[pe_to_report[n] for n in xrange(len(pe_to_report)) if pe_density[n] == m] for m in xrange(0,POLY_PE_STEPS+1)]
		pe_to_report_stratified = [sorted_coord_list_to_ranges(n) for n in pe_to_report_stratified]
		#mpl.figure(6)
		#mpl.plot(pe_to_report, pe_density)
		#mpl.savefig('debug.png')
		#exit(1)

		if (len(args.c) or len(COMPARE)) and myPECount >= MIN_DISC_ONLY:
			pe_coords_to_compare = sorted_coord_list_to_ranges([pe_to_report[n] for n in xrange(len(pe_to_report)) if pe_density[n] >= POLY_PE_STEPS-1])
			pe_coords_to_compare = [int(np.mean(n)+0.5) for n in pe_coords_to_compare]
			for pem in pe_coords_to_compare:
				mapq0_percent = 'N/A'
				if pem in pe_mapq_byPos:
					mapq0_frac    = float(pe_mapq_byPos[pem].count(0))/len(pe_mapq_byPos[pem])
					mapq0_percent = '{0:0.2f}%'.format(100.*mapq0_frac)
				myPECount_2 = myPECount
				if pem in pe_count:
					myPECount_2 = pe_count[pem]

				gcl = get_compare(clustered_events[i][0][0], pem)
				if len(gcl) > 0:
					for gc in gcl:
						COMPARE_OUT[gc].append((clustered_events[i][0][0], pem, str(myPECount_2)+'/'+str(myPECount), mapq0_percent, abs(pem-gc[1]), 'DISCORDANT'))
				#print('PEM:', pem, mapq0_percent)

		pe_to_report = sorted_coord_list_to_ranges(pe_to_report)

	# read support filters
	if len(evidence_pb[i]) == 0:
		if MIN_SOFTCLIP > 0 and max_sc < MIN_SOFTCLIP:
			not_enough_softclip = True
		else:
			not_enough_softclip = False
		if MIN_DISC_ONLY > 0 and max_pe < MIN_DISC_ONLY:
			not_enough_disc = True
		else:
			not_enough_disc = False
		if not_enough_softclip:
			continue
		if MIN_SOFTCLIP <= 0 and not_enough_disc:
			continue
		if MIN_LONG_PER_EVENT > 0:	# we have sufficient short read evidence, but no pacbio (and we want pacbio!)
			continue

	elif len(evidence_pb[i]) >= MIN_LONG_PER_EVENT:
		ccs_coord_list = [n[0] for n in evidence_pb[i] if n[1] == 'CCS']
		if len(ccs_coord_list) > 1:
			scm = np.median(ccs_coord_list)
			scl = [int(abs(n-scm)+0.5) for n in ccs_coord_list]
			ccs_str = 'type: longread-ccs count: ' + str(len(scl)) + ' avg_pos: ' + str(int(scm+0.5)) + ' std: ' + str(np.mean(scl))
			bp_dev_ccs.extend(scl)
		elif len(ccs_coord_list) == 1:
			ccs_str = 'type: longread-ccs count: ' + str(1) + ' avg_pos: ' + str(ccs_coord_list[0]) + ' std: 0.00'

		clr_coord_list = [n[0] for n in evidence_pb[i] if n[1] == 'CLR']
		if len(clr_coord_list) > 1:
			scm = np.median(clr_coord_list)
			scl = [int(abs(n-scm)+0.5) for n in clr_coord_list]
			clr_str = 'type: longread-clr count: ' + str(len(scl)) + ' avg_pos: ' + str(int(scm+0.5)) + ' std: ' + str(np.mean(scl))
			bp_dev_clr.extend(scl)
		elif len(clr_coord_list) == 1:
			clr_str = 'type: longread-clr count: ' + str(1) + ' avg_pos: ' + str(clr_coord_list[0]) + ' std: 0.00'

	#
	#	PLOTTING
	#
	PWIDTH = 0.3
	PWIDTH_SC = 0.4
	poly_pe = []
	poly_pe_stratified = [[] for n in pe_to_report_stratified]
	poly_sc = []
	poly_pb1 = []
	poly_pb2 = []
	if max_pe:
		for n in pe_to_report:
			poly_pe.append(Polygon(np.array([[n[0]-0.5,-PWIDTH+2], [n[0]-0.5,PWIDTH+2], [n[1]+0.5,PWIDTH+2], [n[1]+0.5,-PWIDTH+2]]), closed=True, linewidth=0))
		for j in xrange(len(poly_pe_stratified)):
			for n in pe_to_report_stratified[j]:
				poly_pe_stratified[j].append(Polygon(np.array([[n[0]-0.5,-PWIDTH+2], [n[0]-0.5,PWIDTH+2], [n[1]+0.5,PWIDTH+2], [n[1]+0.5,-PWIDTH+2]]), closed=True, linewidth=0))
	if max_sc:
		for n in sc_to_report:
			poly_sc.append(Polygon(np.array([[n[0]-0.5,-PWIDTH_SC+2], [n[0]-0.5,PWIDTH_SC+2], [n[1]+0.5,PWIDTH_SC+2], [n[1]+0.5,-PWIDTH_SC+2]]), closed=True, linewidth=2))
	if len(evidence_pb[i]):
		for n in evidence_pb[i]:
			if n[1] == 'CCS':
				poly_pb1.append(Polygon(np.array([[n[0]-0.5,-PWIDTH+1], [n[0]-0.5,PWIDTH+1], [n[0]+0.5,PWIDTH+1], [n[0]+0.5,-PWIDTH+1]]), closed=True, linewidth=2))
			else:
				poly_pb2.append(Polygon(np.array([[n[0]-0.5,-PWIDTH+1], [n[0]-0.5,PWIDTH+1], [n[0]+0.5,PWIDTH+1], [n[0]+0.5,-PWIDTH+1]]), closed=True, linewidth=2))

	mpl.rcParams.update({'font.size': 16, 'font.weight':'bold'})

	fig = mpl.figure(0,figsize=(12,6))
	ax = mpl.gca()
	#ax.add_collection(PatchCollection(poly_pe, alpha=0.6, color='tab:blue', linewidth=0))
	for j in xrange(len(poly_pe_stratified)):
		ax.add_collection(PatchCollection(poly_pe_stratified[j], alpha=POLY_PE_ALPHA[j], color='blue', linewidth=0))
	ax.add_collection(PatchCollection(poly_sc, alpha=0.6, color='black', linewidth=2))
	ax.add_collection(PatchCollection(poly_pb2, alpha=0.6, color='purple', linewidth=2))
	ax.add_collection(PatchCollection(poly_pb1, alpha=0.6, color='red', linewidth=2))

	allv  = []
	legText = []
	if max_pe:
		allv += [n[0] for n in pe_to_report] + [n[1] for n in pe_to_report]
		mpl.plot([-2,-1],[0,0],color='blue')
		legText.append('Discordant pairs')
	if max_sc:
		allv += [n[0] for n in sc_to_report] + [n[1] for n in sc_to_report]
		mpl.plot([-2,-1],[0,0],color='black')
		legText.append('Softclipped reads')
	if len(evidence_pb[i]):
		allv += [n[0] for n in evidence_pb[i]]
		mpl.plot([-2,-1],[0,0],color='purple')
		legText.append('PacBio CLR')
		mpl.plot([-2,-1],[0,0],color='red')
		legText.append('PacBio CCS')

	leg = mpl.legend(legText,loc=4)
	for l in leg.legendHandles:            
		l.set_linewidth(10)

	mpl.axis([min(allv)-PLOT_BUFF, max(allv)+PLOT_BUFF, -1, 3])
	mpl.grid(linestyle='--', alpha=0.5)
	mpl.yticks([0,1,2],['', 'Exogene-LR', 'Exogene-SR'])
	mpl.xlabel('Reference coordinates')
	mpl.title(clustered_events[i][0][0] + ' --> ' + clustered_events[i][0][1])
	mpl.tight_layout()

	igv_pos = clustered_events[i][0][0] + ':' + str(min(allv)-PLOT_BUFF) + '-' + str(max(allv)+PLOT_BUFF)

	myClass = (1*(len(poly_sc) > 0), 1*(len(poly_pe) > 0), 1*(len(poly_pb1) > 0), 1*(len(poly_pb2) > 0))
	if myClass not in class_counts:
		class_counts[myClass] = 0
	class_counts[myClass] += 1
	if myClass not in sc_count_by_class:
		sc_count_by_class[myClass] = []
		pe_count_by_class[myClass] = []
	sc_count_by_class[myClass].append(mySCCount)
	pe_count_by_class[myClass].append(myPECount)
	#exit(1)

	#
	# bed intersection with average breakpoint position
	#
	# - this is weird and inefficient. not sure why I did it this way.
	#
	avg_sc, avg_pe, avg_pb = 0,0,0
	if max_sc:
		avg_sc = int(sum([n[0]+n[1] for n in sc_to_report])/(2*len(sc_to_report)))
	if max_pe:
		avg_pe = int(sum([n[0]+n[1] for n in pe_to_report])/(2*len(pe_to_report)))
	if len(evidence_pb[i]):
		avg_pb = int(sum([n[0] for n in evidence_pb[i]])/len(evidence_pb[i]))
	avg_all = [n for n in [avg_sc, avg_pe, avg_pb] if n]
	avg_all = int(sum(avg_all)/len(avg_all))
	bed_chr = clustered_events[i][0][0]
	bed_pos = avg_all

	# annotate!
	bed_out = []
	for bed_i in xrange(len(BED_TRACKS)):
		if BED_TRACKS[bed_i][1].query(bed_chr, bed_pos):
			#print('--',BED_TRACKS[bed_i][0])
			bed_out.append(BED_TRACKS[bed_i][0])

	#
	# PRINT BREAKPOINT SUMMARY
	#
	if len(bed_out) == 0 or REPORT_ALL_BEDHITS:
		print('CLUSTER', str(nPlot)+':', igv_pos)
		print('== VIRUS:     ', clustered_events[i][0][1])
		if max_sc > 0: print('== SOFTCLIP:  ', sc_to_report, max_sc)
		if len(sc_str): print('==== ' + sc_str)
		if max_pe > 0: print('== DISCORDANT:', pe_to_report, max_pe)
		if len(evidence_pb[i]) > 0: print('== PACBIO:    ', evidence_pb[i])
		if len(ccs_str): print('==== ' + ccs_str)
		if len(clr_str): print('==== ' + clr_str)
		if len(bed_out): print('== ANNOTATION:', ' '.join(bed_out))
		print('')

	#
	# DATA FOR OUTPUT REPORT
	#
	# (CHR, INTEGRATION_POS, #READS, VIRUS, ANNOTATION, SOFTCLIP_POS, #SOFTCLIP, DISCORDANT_POS, #DISCORDANT, LONGREAD_POS, #LONGREAD, NEAREST_GENE)
	if len(bed_out) == 0 or REPORT_ALL_BEDHITS:
		out_chr = clustered_events[i][0][0]
		out_pos = '-'
		out_num = str(max_pe + len(evidence_pb[i]))
		out_vir = clustered_events[i][0][1]
		out_ann = '-'
		if len(bed_out):
			out_ann = ','.join(bed_out)
		#
		(out_scp, out_scc) = ('-', '0')
		if len(sc_str):
			out_scp = sc_str.split(' ')[5]
			out_scc = sc_str.split(' ')[3]
			out_pos = sc_str.split(' ')[5]
		#
		out_dip = '-'
		if max_pe > 0:
			out_dip = ','.join([str(n[0])+'-'+str(n[1]) for n in pe_to_report])
		out_dic = str(max_pe)
		#
		out_lrp = '-'
		if len(clr_str):
			out_lrp = clr_str.split(' ')[5]
			if out_pos == '-':
				out_pos = clr_str.split(' ')[5]
		if len(ccs_str):	# prio ccs result over clr
			out_lrp = ccs_str.split(' ')[5]
			out_pos = ccs_str.split(' ')[5]
		out_lrc = str(len(evidence_pb[i]))
		#
		out_nge = '-'
		if out_pos == '-':	# disc only
			backup_pos  = np.mean([np.mean([n[0],n[1]]) for n in pe_to_report])
			nearest_hit = get_nearest_transcript(out_chr, backup_pos, TRANSCRIPT_TRACK)
		else:
			nearest_hit = get_nearest_transcript(out_chr, int(out_pos), TRANSCRIPT_TRACK)
		if nearest_hit != None:
			out_nge = nearest_hit[1][1] + ' (' + nearest_hit[1][0] + ')'
		#
		out_report_data.append((out_chr, out_pos, out_num, out_vir, out_ann, out_scp, out_scc, out_dip, out_dic, out_lrp, out_lrc, out_nge))

		#mpl.show()
		mpl.savefig(OUT_DIR+'site_'+str(nPlot).zfill(zfill_num)+'_'+igv_pos.replace(':','_')+'.png')
		nPlot += 1
	mpl.close(fig)

	BIG_VAL = 9999999
	min_distance_between_sc_and_pb = BIG_VAL
	if max_sc and len(evidence_pb[i]):
		first_clr = True
		first_ccs = True
		for j in xrange(len(sc_to_report)):
			for k in xrange(len(evidence_pb[i])):
				pDist = abs(sc_to_report[j][0] - evidence_pb[i][k][0])
				if evidence_pb[i][k][1] == 'CLR':
					bp_dist_sc_clr.append(pDist)
					if first_clr:
						bp_dist_min_sc_clr.append(pDist)
						first_clr = False
					else:
						bp_dist_min_sc_clr[-1] = min([bp_dist_min_sc_clr[-1], pDist])
				elif evidence_pb[i][k][1] == 'CCS':
					bp_dist_sc_ccs.append(pDist)
					if first_ccs:
						bp_dist_min_sc_ccs.append(pDist)
						first_ccs = False
					else:
						bp_dist_min_sc_ccs[-1] = min([bp_dist_min_sc_ccs[-1], pDist])
				if pDist < min_distance_between_sc_and_pb:
					min_distance_between_sc_and_pb = pDist
	#if min_distance_between_sc_and_pb < BIG_VAL:
	#	print(igv_pos,'SC_PB_DIST =',min_distance_between_sc_and_pb)

#
# WRITE OUTPUT REPORT
#
header = ('CHR', 'INTEGRATION_POS', '#READS', 'VIRUS', 'ANNOTATION',
          'SOFTCLIP_POS', '#SOFTCLIP',
          'DISCORDANT_POS', '#DISCORDANT',
          'LONGREAD_POS', '#LONGREAD',
          'NEAREST_GENE')
f = open(OUT_DIR+'integrations.txt', 'w')
f.write('\t'.join(header)+'\n')
for n in out_report_data:
	f.write('\t'.join(n)+'\n')
f.close()

#
# WRITE OUTPUT BED
#
f = open(OUT_DIR+'integrations.bed', 'w')
for n in out_report_data:
	f.write(n[0] + '\t' + n[1] + '\t' + str(int(n[1])+1) + '\t')
	if len(SAMPLE_NAME):
		f.write(SAMPLE_NAME + ' ')
	f.write(n[3] + '\n')
f.close()

#
# BREAKPOINT STATS
#
bp_dev_sc  = [n for n in bp_dev_sc if n <= 300]
bp_dev_ccs = [n for n in bp_dev_ccs if n <= 300]
bp_dev_clr = [n for n in bp_dev_clr if n <= 300]
print('=== BREAKPOINT DEVIATIONS:')
if len(IN_SHORT):
	if len(bp_dev_sc):
		print('short:', np.mean(bp_dev_sc))
if len(IN_LONG):
	if len(bp_dev_ccs):
		print('ccs:  ', np.mean(bp_dev_ccs))
	if len(bp_dev_clr):
		print('clr:  ', np.mean(bp_dev_clr))

#kv = sorted([(class_counts[k],k) for k in class_counts.keys()], reverse=True)
#for n in kv:
#	print(n[0], n[1], 'sc_counts:', sc_count_by_class[n[1]], 'pe_counts:', pe_count_by_class[n[1]])

#
# LONG READ BREAKPOINT STATS
#
if len(IN_SHORT) and len(IN_LONG):
	print('')
	if len(bp_dist_sc_clr):
		print('avg sc vs. clr:', bp_dist_sc_clr, np.mean(bp_dist_sc_clr))
	if len(bp_dist_sc_ccs):
		print('avg sc vs. ccs:', bp_dist_sc_ccs, np.mean(bp_dist_sc_ccs))
	if len(bp_dist_min_sc_clr):
		print('min sc vs. clr:', bp_dist_min_sc_clr, np.mean(bp_dist_min_sc_clr))
	if len(bp_dist_min_sc_ccs):
		print('min sc vs. ccs:', bp_dist_min_sc_ccs, np.mean(bp_dist_min_sc_ccs))
	print('')

#
# SRA COMPARISON
#
if len(args.c) or len(COMPARE):
	print('')
	print('### START_COMPARE ###')
	for n in COMPARE:
		if len(COMPARE_OUT[n]) == 0:
			anno = [isInBadRange(n[0],n[1])*'gap', BED_TRACKS[2][1].query(n[0],n[1])*'repeats', BED_TRACKS[3][1].query(n[0],n[1])*'mappability', BED_TRACKS[4][1].query(n[0],n[1])*'exclude']
			print('MISS\t' + str(n) + '\t' + ', '.join([m for m in anno if len(m)]))
		else:
			min_distance_co_item = sorted([(COMPARE_OUT[n][i][4], i) for i in xrange(len(COMPARE_OUT[n]))])[0][1]
			s1 = ' '.join([str(m) for m in COMPARE_OUT[n][min_distance_co_item][0:6]])
			s2 = str(n)
			#s3 = str(COMPARE_OUT[n])
			print('HIT\t' + s1 + '\t' + s2)
	print('\n### NOT_IN_COMPARE ###\n')
	for n in COMPARE_OUT_FP:
		print(n[0], n[1], n[2], n[3])
	print('')
