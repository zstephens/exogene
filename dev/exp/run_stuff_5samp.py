import os

WORKING_DIR = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_5new/'
SRA_FQ_DIR  = '/research/bsi/data/controlled_access/tcga/LIHC/'

BAM_DIR = '/research/bsi/data/controlled_access/tcga/LIHC/'
OUT_DIR = WORKING_DIR + 'out/'
GIT_DIR = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_sra/git/exogene/dev/'

OUT_QSH  = WORKING_DIR + 'qsh/'
OUT_QLOG = WORKING_DIR + 'qlog/'
EMAIL    = 'stephens.zachary@mayo.edu'

QUEUES  = ['1-hour', '1-day', '4-day', '7-day', '30-day']
myQueue = '1-day'

SAMPLES = [['9b605411-e103-4b26-b271-faa681ae3829', 'TCGA-DD-AACV-01A-11D-A40R-10'],
           ['faecbc9e-693a-45a3-b041-79b17c75c880', 'TCGA-DD-AACV-01A-11R-A41C-07'],
           ['f9764955-0347-43dc-9e7c-72501eec5c88', 'TCGA-DD-AACV-01A-11R-A41H-13'],
           ['5f15cedd-18b8-4082-8c30-a4db8551b17e', 'TCGA-DD-AACV-10A-01D-A40U-10'],
           ['958b0375-51ca-47e1-9bfc-424afba91cd6', 'TCGA-DD-AAD0-01A-11D-A40R-10'],
           ['ee1cbcf2-706f-4a6d-9cd8-d25666d3f093', 'TCGA-DD-AAD0-01A-11R-A41C-07'],
           ['260520b6-080c-4917-9001-35427b21dceb', 'TCGA-DD-AAD0-01A-11R-A41H-13'],
           ['56108e1f-849c-4d4d-9bff-0205180e5c94', 'TCGA-DD-AAD0-10A-01D-A40U-10'],
           ['1397d066-01c2-4cd6-9cda-6fbb95887be0', 'TCGA-DD-AADL-01A-11D-A40R-10'],
           ['9c516342-1859-45db-9202-23a9469e43f9', 'TCGA-DD-AADL-01A-11R-A41C-07'],
           ['388893f8-1fc6-4180-b0f6-6deb3b20bd48', 'TCGA-DD-AADL-01A-11R-A41H-13'],
           ['9c404bcc-59a0-4116-ab02-eebfefb8677b', 'TCGA-DD-AADL-10A-01D-A40U-10'],
           ['7552eb82-ff2a-4cc0-a0eb-45dbadc6136f', 'TCGA-DD-AADU-01A-11D-A40R-10'],
           ['802206f5-b7c3-4262-88d7-a7a373c70202', 'TCGA-DD-AADU-01A-11R-A41C-07'],
           ['92baa5a4-c31b-4aae-9e2e-c171e6430feb', 'TCGA-DD-AADU-01A-11R-A41H-13'],
           ['8238e3ae-6936-40be-8639-a32492b2e8fe', 'TCGA-DD-AADU-10A-01D-A40U-10'],
           ['ba001ff0-4895-45e6-85c6-58eb3b016198', 'TCGA-DD-AADV-01A-11D-A38X-10'],
           ['a376ce6f-6d63-45cc-83d0-8910f0c6161e', 'TCGA-DD-AADV-01A-11R-A39D-07'],
           ['fb880503-4b4c-4d7f-8c94-61954cd12259', 'TCGA-DD-AADV-01A-11R-A39J-13'],
           ['d2a4c924-f4c2-42d4-a3df-660748dcd440', 'TCGA-DD-AADV-10A-01D-A38X-10']]

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

for i in xrange(len(SAMPLES)):
	listing = [n for n in os.listdir(BAM_DIR+SAMPLES[i]+'/') if n[-4:] == '.bam']
	if len(listing) and exists_and_is_nonZero(listing[0]):
		myBam = BAM_DIR+SAMPLES[i][0]+'/'+listing[0]
	else:
		continue
	jobName = 'exogene_' + str(i)
	exo_out = OUT_DIR + SAMPLES[i][1] + '/'

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

	reads_report = exo_out + 'Viral_Reads_Report.tsv'
	plots_out    = exo_out + 'plots/'
	report_out   = exo_out + 'Integration_Summary.txt'

	CMD  = ''
	CMD += EXOGENE_SR + ' -b ' + myBam + ' -r ' + REF_HVR38 + ' -o ' + exo_out + '\n'
	CMD += COMBINE_REP + ' -ms 2 -s ' + reads_report + ' -o ' + plots_out + ' > ' + report_out + '\n'

	if runCMD:
		f = open(OUT_QSH+jobName+'.sh', 'w')
		f.write(HEADER+'\n'+CMD+'\n')
		f.close()
		print OUT_QSH+jobName+'.sh'
		#os.system('qsub '+OUT_QSH+jobName+'.sh')



####listing = [n for n in os.listdir(SRA_FQ_DIR) if n[-9:] == '.fastq.gz']
####fq_dict = {n.split('_')[0]:[None,None] for n in listing}
####for n in listing:
####	if n.split('_')[-1][0] == '1':
####		fq_dict[n.split('_')[0]][0] = SRA_FQ_DIR + n
####	elif n.split('_')[-1][0] == '2':
####		fq_dict[n.split('_')[0]][1] = SRA_FQ_DIR + n
####
####nProcessed = 0
####for sampleName in sorted(fq_dict.keys()):
####	r1 = fq_dict[sampleName][0]
####	r2 = fq_dict[sampleName][1]
####	if r1 == None or r2 == None:
####		continue
####	jobName = 'exogene_' + sampleName
####	exo_out = OUT_DIR + jobName + '/'
####
####	#if sampleName not in REPROCESS:
####	#	continue
####
####	HEADER  = ''
####	HEADER += '#!/bin/bash' + '\n'
####	HEADER += '#$ -N ' + jobName + '\n'
####	HEADER += '#$ -q ' + myQueue + '\n'
####	HEADER += '#$ -o ' + OUT_QLOG + jobName + '.o' + '\n'
####	HEADER += '#$ -e ' + OUT_QLOG + jobName + '.e' + '\n'
####	#HEADER += '#$ -m ae' + '\n'
####	#HEADER += '#$ -M ' + EMAIL + '\n'
####	HEADER += '#$ -l h_vmem=32G ' + '\n'
####	HEADER += '#$ -notify ' + '\n'
####
####	r1_clean = FQ_DIR + r1.split('/')[-1]
####	r2_clean = FQ_DIR + r2.split('/')[-1]
####	exo_out  = OUT_DIR + jobName + '/'
####
####	reads_report = exo_out + 'Viral_Reads_Report.tsv'
####	plots_out    = exo_out + 'plots/'
####	report_out   = exo_out + 'Integration_Summary.txt'
####
####	runCMD = False
####	
####	CMD  = ''
####	if exists_and_is_nonZero(r1_clean) == False or exists_and_is_nonZero(r2_clean) == False:
####		CMD += PREP_SRA_FQ + ' ' + r1 + ' ' + r1_clean + ' 1 &' + '\n'
####		CMD += PREP_SRA_FQ + ' ' + r2 + ' ' + r2_clean + ' 2 &' + '\n'
####		CMD += 'wait' + '\n'
####	#
####	haveData = False
####	if exists_and_is_nonZero(reads_report):
####			f = open(reads_report,'r')
####			fr = f.read()
####			f.close()
####			nlc = fr.count('\n')
####			if nlc > 2:
####				haveData = True
####	if haveData == False:
####		CMD += EXOGENE_SR + ' -f1 ' + r1_clean + ' -f2 ' + r2_clean + ' -r ' + REF_HVR38 + ' -o ' + exo_out + '\n'
####	#
####	if True or exists_and_is_nonZero(report_out) == False:
####		CMD += COMBINE_REP + ' -ms 20 -s ' + reads_report + ' -o ' + plots_out + ' -c ' + sampleName + ' > ' + report_out + '\n'
####		runCMD = True
####
####	if runCMD:
####		f = open(OUT_QSH+jobName+'.sh', 'w')
####		f.write(HEADER+'\n'+CMD+'\n')
####		f.close()
####		os.system('qsub '+OUT_QSH+jobName+'.sh')
####	nProcessed += 1
####
####	if nProcessed >= 10000:
####		break
