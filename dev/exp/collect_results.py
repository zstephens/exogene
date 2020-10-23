import os

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

WORKING_DIR = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_sra/'

OUT_DIR = WORKING_DIR + 'out_redo/'

dat_out = {}

PREFIX = 'exo_'
PLEN   = len(PREFIX)

listing = [n for n in os.listdir(OUT_DIR) if n[:PLEN] == PREFIX]
for lDir in listing:
	listing2  = [n for n in os.listdir(OUT_DIR+lDir) if n == 'Integration_Summary.txt']
	if len(listing2):
		mySummary = OUT_DIR+lDir+'/'+listing2[0]
		#print mySummary, exists_and_is_nonZero(mySummary)
		if exists_and_is_nonZero(mySummary):
			if lDir[PLEN:] not in dat_out:
				dat_out[lDir[PLEN:]] = ''
			startReading = False
			f = open(mySummary,'r')
			for line in f:
				if '### START_COMPARE ###' in line:
					startReading = True
					line = line.strip()+'\t'+lDir[PLEN:]+'\n'
				if startReading:
					dat_out[lDir[PLEN:]] += line
			f.close()
			if startReading == False:
				print 'Warning:', lDir, 'has incomplete Integration_Summary.txt'
		else:
			print 'Warning:', lDir, 'has empty Integration_Summary.txt'
	else:
		print 'Warning:', lDir, 'does not have Integration_Summary.txt'

collected_results = WORKING_DIR + 'collected_results_redo.txt'
f = open(collected_results,'w')
for k in sorted(dat_out.keys()):
	f.write(dat_out[k])
	f.write('\n')
f.close()
