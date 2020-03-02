import os

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

WORKING_DIR = '/research/bsi/projects/PI/tertiary/Kocher_Jean-Pierre_m026645/s205842.Viral_Integration/processing/exogene_sra/'

OUT_DIR = WORKING_DIR + 'out/'

dat_out = {}

listing = [n for n in os.listdir(OUT_DIR) if n[:8] == 'exogene_']
for lDir in listing:
	listing2  = [n for n in os.listdir(OUT_DIR+lDir) if n == 'Integration_Summary.txt']
	if len(listing2):
		mySummary = OUT_DIR+lDir+'/'+listing2[0]
		#print mySummary, exists_and_is_nonZero(mySummary)
		if exists_and_is_nonZero(mySummary):
			if lDir[8:] not in dat_out:
				dat_out[lDir[8:]] = ''
			startReading = False
			f = open(mySummary,'r')
			for line in f:
				if '### START_COMPARE ###' in line:
					startReading = True
					line = line.strip()+'\t'+lDir[8:]+'\n'
				if startReading:
					dat_out[lDir[8:]] += line
			f.write('\n')
			f.close()
	else:
		print 'Warning:', lDir, 'does not have Integration_Summary.txt'

collected_results = WORKING_DIR + 'collected_results.txt'
f = open(collected_results,'w')
for k in sorted(dat_out.keys()):
	f.write(dat_out[k])
f.close()
