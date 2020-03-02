import os

WORKING_DIR = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_sra/'
SRA_FQ_DIR  = '/research/bsi/data/controlled_access/dbGAP/PRJNA298941/fastqs/'

FQ_DIR  = WORKING_DIR + 'fq/'
OUT_DIR = WORKING_DIR + 'out/'
GIT_DIR = WORKING_DIR + 'git/exogene/dev/'

OUT_QSH  = WORKING_DIR + 'qsh/'
OUT_QLOG = WORKING_DIR + 'qlog/'
EMAIL    = 'stephens.zachary@mayo.edu'

QUEUES  = ['1-hour', '1-day', '4-day', '7-day', '30-day']
myQueue = '1-day'

PYTHON      = '/research/bsi/tools/biotools/smrtlink/5.1.0/smrtcmds/bin/python2.7'
PREP_SRA_FQ = PYTHON + ' ' + GIT_DIR + 'prep_sra_fq.py'
EXOGENE_SR  = GIT_DIR + 'Exogene-SR-mforge.sh'
COMBINE_REP = PYTHON + ' ' + GIT_DIR + 'combine_reports.py'

REF_HVR38   = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_sra/ref/hg38_and_viral.fa'

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

listing = [n for n in os.listdir(SRA_FQ_DIR) if n[-9:] == '.fastq.gz']
fq_dict = {n.split('_')[0]:[None,None] for n in listing}
for n in listing:
	if n.split('_')[-1][0] == '1':
		fq_dict[n.split('_')[0]][0] = SRA_FQ_DIR + n
	elif n.split('_')[-1][0] == '2':
		fq_dict[n.split('_')[0]][1] = SRA_FQ_DIR + n

for sampleName in sorted(fq_dict.keys()):
	r1 = fq_dict[sampleName][0]
	r2 = fq_dict[sampleName][1]
	if r1 == None or r2 == None:
		continue
	jobName = 'exogene_' + sampleName

	HEADER  = ''
	HEADER += '#!/bin/bash' + '\n'
	HEADER += '#$ -N ' + jobName + '\n'
	HEADER += '#$ -q ' + myQueue + '\n'
	HEADER += '#$ -o ' + OUT_QLOG + jobName + '.o' + '\n'
	HEADER += '#$ -e ' + OUT_QLOG + jobName + '.e' + '\n'
	#HEADER += '#$ -m ae' + '\n'
	#HEADER += '#$ -M ' + EMAIL + '\n'
	HEADER += '#$ -l h_vmem=32G ' + '\n'
	HEADER += '#$ -notify ' + '\n'

	r1_clean = FQ_DIR + r1.split('/')[-1]
	r2_clean = FQ_DIR + r2.split('/')[-1]
	exo_out  = OUT_DIR + jobName + '/'

	reads_report = exo_out + 'Viral_Reads_Report.tsv'
	plots_out    = exo_out + 'plots/'
	report_out   = exo_out + 'Integration_Summary.txt'
	
	CMD  = ''
	if exists_and_is_nonZero(r1_clean) == False or exists_and_is_nonZero(r2_clean) == False:
		CMD += PREP_SRA_FQ + ' ' + r1 + ' ' + r1_clean + ' 1 &' + '\n'
		CMD += PREP_SRA_FQ + ' ' + r2 + ' ' + r2_clean + ' 2 &' + '\n'
		CMD += 'wait' + '\n'
	#
	haveData = False
	if exists_and_is_nonZero(reads_report):
			f = open(reads_report,'r')
			fr = f.read()
			f.close()
			nlc = fr.count('\n')
			if nlc > 2:
				haveData = True
	if haveData == False:
		CMD += EXOGENE_SR + ' -f1 ' + r1_clean + ' -f2 ' + r2_clean + ' -r ' + REF_HVR38 + ' -o ' + exo_out + '\n'
	#
	if exists_and_is_nonZero(report_out) == False:
		CMD += COMBINE_REP + ' -ms 20 -s ' + reads_report + ' -o ' + plots_out + ' -c ' + sampleName + ' > ' + report_out + '\n'

	f = open(OUT_QSH+jobName+'.sh', 'w')
	f.write(HEADER+'\n'+CMD+'\n')
	f.close()

	break
