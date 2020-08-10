import os

WORKING_DIR = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_5new/'
#SRA_FQ_DIR  = '/research/bsi/data/controlled_access/tcga/LIHC/'

BAM_DIR = '/research/bsi/data/controlled_access/tcga/LIHC/'
OUT_DIR = WORKING_DIR + 'out/'
GIT_DIR = WORKING_DIR + 'exogene/dev/'
HGT_OUT = WORKING_DIR + 'hgtid_result/'
REALIGN_DIR = WORKING_DIR + 'bam_realign_to_hg38/'

OUT_QSH  = WORKING_DIR + 'qsh/'
OUT_QLOG = WORKING_DIR + 'qlog/'
EMAIL    = 'stephens.zachary@mayo.edu'
EMAIL    = 'zstephe2@illinois.edu'

QUEUES  = ['1-hour', '1-day', '4-day', '7-day', '30-day']
myQueue = '4-day'

SAMPLES = [['9b605411-e103-4b26-b271-faa681ae3829', 'TCGA-DD-AACV-01A-11D-A40R-10', 'WXS',   'TUMOR'],
           #['82f4bdfc-6861-11e4-9b49-60d7fc308db6', 'TCGA-DD-AACV-01A-11R-A41C-07', 'RNA',   'TUMOR'],
           #['f9764955-0347-43dc-9e7c-72501eec5c88', 'TCGA-DD-AACV-01A-11R-A41H-13', 'miRNA', 'TUMOR'],
           ['5f15cedd-18b8-4082-8c30-a4db8551b17e', 'TCGA-DD-AACV-10A-01D-A40U-10', 'WXS',   'NORMAL'],

           ['958b0375-51ca-47e1-9bfc-424afba91cd6', 'TCGA-DD-AAD0-01A-11D-A40R-10', 'WXS',   'TUMOR'],
           #['a463e1e4-7f21-11e4-bc03-c46bfc308db6', 'TCGA-DD-AAD0-01A-11R-A41C-07', 'RNA',   'TUMOR'],
           #['260520b6-080c-4917-9001-35427b21dceb', 'TCGA-DD-AAD0-01A-11R-A41H-13', 'miRNA', 'TUMOR'],
           ['56108e1f-849c-4d4d-9bff-0205180e5c94', 'TCGA-DD-AAD0-10A-01D-A40U-10', 'WXS',   'NORMAL'],

           ['1397d066-01c2-4cd6-9cda-6fbb95887be0', 'TCGA-DD-AADL-01A-11D-A40R-10', 'WXS',   'TUMOR'],
           #['a4b0985e-7f21-11e4-bc03-c46bfc308db6', 'TCGA-DD-AADL-01A-11R-A41C-07', 'RNA',   'TUMOR'],
           #['388893f8-1fc6-4180-b0f6-6deb3b20bd48', 'TCGA-DD-AADL-01A-11R-A41H-13', 'miRNA', 'TUMOR'],
           ['9c404bcc-59a0-4116-ab02-eebfefb8677b', 'TCGA-DD-AADL-10A-01D-A40U-10', 'WXS',   'NORMAL'],

           ['7552eb82-ff2a-4cc0-a0eb-45dbadc6136f', 'TCGA-DD-AADU-01A-11D-A40R-10', 'WXS',   'TUMOR'],
           #['802206f5-b7c3-4262-88d7-a7a373c70202', 'TCGA-DD-AADU-01A-11R-A41C-07', 'RNA',   'TUMOR'],
           #['92baa5a4-c31b-4aae-9e2e-c171e6430feb', 'TCGA-DD-AADU-01A-11R-A41H-13', 'miRNA', 'TUMOR'],
           ['8238e3ae-6936-40be-8639-a32492b2e8fe', 'TCGA-DD-AADU-10A-01D-A40U-10', 'WXS',   'NORMAL'],

           ['ba001ff0-4895-45e6-85c6-58eb3b016198', 'TCGA-DD-AADV-01A-11D-A38X-10', 'WXS',   'TUMOR'],
           #['a376ce6f-6d63-45cc-83d0-8910f0c6161e', 'TCGA-DD-AADV-01A-11R-A39D-07', 'RNA',   'TUMOR'],
           #['fb880503-4b4c-4d7f-8c94-61954cd12259', 'TCGA-DD-AADV-01A-11R-A39J-13', 'miRNA', 'TUMOR'],
           ['d2a4c924-f4c2-42d4-a3df-660748dcd440', 'TCGA-DD-AADV-10A-01D-A38X-10', 'WXS',   'NORMAL'],

           #['b7fcfb70-5390-4b50-b5cb-3308fb098226', 'TCGA-DD-A1EL-01A-11D-A152-10', 'WGS',   'TUMOR'],
           #['70ddc5d6-43d5-43e3-bbd4-33c905d7d5d8', 'TCGA-DD-A1EL-10A-01D-A152-10', 'WGS',   'NORMAL'],
           ['53668b25-f770-4bfd-b3e9-8fb5cb19054e', 'TCGA-DD-A1EL-01A-11D-A152-10', 'WXS',   'TUMOR'],
           ['71d753c8-2aa1-46b3-bb6a-edab9a90de6d', 'TCGA-DD-A1EL-10A-01D-A152-10', 'WXS',   'NORMAL']]

PYTHON      = '/research/bsi/tools/biotools/smrtlink/8.0/bin/smrtcmds/bin/python'
PREP_SRA_FQ = PYTHON + ' ' + GIT_DIR + 'prep_sra_fq.py'
EXOGENE_SR  = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_new/git/exogene/dev/Exogene-SR-mforge.sh'
COMBINE_REP = PYTHON + ' ' + GIT_DIR + 'combine_reports.py'

SAMTOOLS    = '/research/bsi/tools/biotools/samtools/1.10/bin/samtools'
BWA_MEM     = '/research/bsi/tools/biotools/bwa/0.7.10/bwa mem'
PERL        = '/research/bsi/tools/biotools/perl/5.30.0/bin/perl'

REF_HVR38   = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_sra/ref/hg38_and_viral.fa'
REF_HG38    = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/hgt-id/temp/hg38_clean.fa'

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

for i in xrange(len(SAMPLES)):
	listing = [n for n in os.listdir(BAM_DIR+SAMPLES[i][0]+'/') if n[-4:] == '.bam']
	#print i, listing
	if len(listing):
		myBam = BAM_DIR+SAMPLES[i][0]+'/'+listing[0]
	else:
		print SAMPLES[i][0], 'does not exist, sadly.'
		continue
	jobName = 'exogene_' + str(i)
	exo_out = OUT_DIR + SAMPLES[i][1] + '/'

	HEADER  = ''
	HEADER += '#!/bin/bash' + '\n'
	HEADER += '#$ -N ' + jobName + '\n'
	HEADER += '#$ -q ' + myQueue + '\n'
	HEADER += '#$ -o ' + OUT_QLOG + jobName + '.o' + '\n'
	HEADER += '#$ -e ' + OUT_QLOG + jobName + '.e' + '\n'
	HEADER += '#$ -m abe' + '\n'
	HEADER += '#$ -M ' + EMAIL + '\n'
	HEADER += '#$ -l h_vmem=64G ' + '\n'
	HEADER += '#$ -notify ' + '\n'

	reads_report = exo_out + 'Viral_Reads_Report.tsv'
	plots_out    = exo_out + 'plots/'
	report_out   = exo_out + 'Integration_Summary.txt'

	CMD  = ''
	CMD += EXOGENE_SR + ' -b ' + myBam + ' -r ' + REF_HVR38 + ' -o ' + exo_out + '\n'
	CMD += COMBINE_REP + ' -ms 2 -s ' + reads_report + ' -o ' + plots_out + ' > ' + report_out + '\n'

	f = open(OUT_QSH+jobName+'.sh', 'w')
	f.write(HEADER+'\n'+CMD+'\n')
	f.close()
	print OUT_QSH+jobName+'.sh'
	#os.system('qsub '+OUT_QSH+jobName+'.sh')





	jobName = 'realign_' + str(i)

	HEADER  = ''
	HEADER += '#!/bin/bash' + '\n'
	HEADER += '#$ -N ' + jobName + '\n'
	HEADER += '#$ -q ' + '30-day' + '\n'
	HEADER += '#$ -o ' + OUT_QLOG + jobName + '.o' + '\n'
	HEADER += '#$ -e ' + OUT_QLOG + jobName + '.e' + '\n'
	HEADER += '#$ -m abe' + '\n'
	HEADER += '#$ -M ' + EMAIL + '\n'
	HEADER += '#$ -l h_vmem=256G ' + '\n'
	HEADER += '#$ -notify ' + '\n'

	realign_bam  = REALIGN_DIR + SAMPLES[i][1] + '.bam'
	realign_temp = REALIGN_DIR + SAMPLES[i][1] + '-temp'
	realign_sort = REALIGN_DIR + SAMPLES[i][1] + '-sort.bam'

	CMD  = ''
	CMD += SAMTOOLS + ' view -F 2048 $ARG_BAM | '
	CMD += PERL + ' -lane \'print"\\@$F[0]\\n$F[9]\\n+\\n$F[10]";\' | '
	CMD += BWA_MEM + ' -t 4 ' + REF_HG38 + ' - | '
	CMD += SAMTOOLS + ' view -bS - > ' + realign_bam + ';\n\n'

	CMD += SAMTOOLS + ' sort -@ 3 -o ' + realign_sort + ' -T ' + realign_temp + ' ' + realign_bam + ';\n\n'

	CMD += SAMTOOLS + ' index ' + realign_sort + ';\n\n'

	f = open(OUT_QSH+jobName+'.sh', 'w')
	f.write(HEADER+'\n'+CMD+'\n')
	f.close()
	print OUT_QSH+jobName+'.sh'
	#os.system('qsub '+OUT_QSH+jobName+'.sh')






	jobName = 'hgt_' + str(i)

	HEADER  = ''
	HEADER += '#!/bin/bash' + '\n'
	HEADER += '#$ -N ' + jobName + '\n'
	HEADER += '#$ -q ' + '30-day' + '\n'
	HEADER += '#$ -o ' + OUT_QLOG + jobName + '.o' + '\n'
	HEADER += '#$ -e ' + OUT_QLOG + jobName + '.e' + '\n'
	HEADER += '#$ -m abe' + '\n'
	HEADER += '#$ -M ' + EMAIL + '\n'
	HEADER += '#$ -l h_vmem=256G ' + '\n'
	HEADER += '#$ -notify ' + '\n'

	hgtid = HGT_OUT + SAMPLES[i][1] + '/'

	CMD  = ''
	CMD += 'perl /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/hgt-id/src/hgt.pl -d '
	CMD += '-c /research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/hgt-id/src/config-mforge.txt '
	CMD += '-o ' + hgtid + ' -b ' + realign_sort + ';\n\n'

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
