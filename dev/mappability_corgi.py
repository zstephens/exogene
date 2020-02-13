import os
import gzip
import bisect
import cPickle as pickle

LARGE_NUMBER = (2**31) - 1

SKIP_XY = ['X','Y','chrX','chrY']

def exists_and_is_nonZero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def condenseListOfRegions(l):
	delList = [False]
	prevR   = l[0]
	for i in xrange(1,len(l)):
		if l[i][1] <= prevR[1]:
			delList.append(True)
		else:
			delList.append(False)
			prevR = l[i]
	condensedInput = [l[i] for i in xrange(len(l)) if delList[i] == False]
	prevR = condensedInput[-1]
	delList = [False]
	for i in xrange(len(condensedInput)-2,-1,-1):
		if condensedInput[i][0] == prevR[0] and condensedInput[i][1] <= prevR[1]:
			delList.append(True)
		else:
			delList.append(False)
			prevR = condensedInput[i]
	delList.reverse()
	return [condensedInput[i] for i in xrange(len(condensedInput)) if delList[i] == False]

# mappability track
class MappabilityTrack:
	def __init__(self, bedfile, bed_buffer=0):

		#print '-', bedfile

		readFromPickle = False
		if exists_and_is_nonZero(bedfile+'.p'):
			#print 'reading data from pickle...'
			self.all_tracks = pickle.load(open(bedfile+'.p','rb'))
			readFromPickle = True

		if readFromPickle == False:

			if bedfile[-3:] == '.gz':
				f = gzip.open(bedfile,'r')
			else:
				f = open(bedfile,'r')

			# preprocess lines (remove redundancies and join overlaps)
			lines         = {}
			for line in f:
				splt = line.strip().split('\t')
				##### skip certain refs
				####if splt[0] in SKIP_XY:
				####	continue
				####print splt
				myChr = splt[0]
				myPos = max([0,int(splt[1])-bed_buffer])
				myEnd = int(splt[2])+bed_buffer
				if myEnd <= myPos:
					print 'skipping invalid bed region:', [myChr,myPos,myEnd]
					continue
				if myChr not in lines:
					lines[myChr] = []
				lines[myChr].append((myPos,myEnd))
			f.close()
			for k in sorted(lines.keys()):
				print k, len(lines[k]), '-->',
				s_dat = condenseListOfRegions(sorted(lines[k]))
				print len(s_dat), '-->',
				can_I_stop_yet = False
				i_start = len(s_dat)-2
				while can_I_stop_yet == False:
					can_I_stop_yet = True
					for i in xrange(i_start+1,0,-1):
						if s_dat[i-1][1] >= s_dat[i][0]:
							s_dat[i-1] = (s_dat[i-1][0],s_dat[i][1])
							del s_dat[i]
							i_start = min([i,len(s_dat)-2])
							can_I_stop_yet = False
							break
				lines[k] = condenseListOfRegions(s_dat)
				print len(lines[k])

			# convert to a bisect-friendly list
			self.all_tracks = {k:[-1] for k in lines.keys()}
			endSoFar        = {k:-1 for k in lines.keys()}
			for k in lines.keys():
				for dat in lines[k]:
					[myChr,myPos,myEnd] = [k,dat[0],dat[1]]
					if myPos < endSoFar[myChr]:
						print 'skipping unsorted bed region:', [myChr,myPos,myEnd], endSoFar[myChr]
						continue
					endSoFar[myChr] = myEnd
					self.all_tracks[myChr].append(myPos)
					self.all_tracks[myChr].append(myEnd)
				# append something huge to prevent bisect from going past the length of the list
				if LARGE_NUMBER <= self.all_tracks[k][-1]:
					print "\nYour large number isn't large enough!\n"
					exit(1)
				self.all_tracks[k].append(LARGE_NUMBER)

			print 'saving bed data as pickle...'
			pickle.dump(self.all_tracks,open(bedfile+'.p','wb'))

	# return True if index is in bedfile region, False if outside
	def query(self, myChr, coord, endPointInclusive=True):
		if myChr in self.all_tracks:
			bcoord = bisect.bisect(self.all_tracks[myChr],coord)
			if bcoord%2 == 0:
				return True
			if endPointInclusive and coord == self.all_tracks[myChr][bcoord-1]:
				return True
		return False

	# returns number of positions of a range that lie in bedfile regions
	def query_range(self, myChr, coord1, coord2, query_endPointInclusive=False):
		count = 0
		if myChr in self.all_tracks:
			for i in xrange(coord1,coord2+1*query_endPointInclusive):
				if self.query(myChr,i):
					count += 1
		return count

	# more efficient version of the above...
	def query_range_faster(self, myChr, coord1, coord2, query_endPointInclusive=False):
		count = 0
		if myChr in self.all_tracks:
			lb = bisect.bisect(self.all_tracks[myChr],coord1)
			ub = bisect.bisect(self.all_tracks[myChr],coord2+1*query_endPointInclusive)
			
			# entirely in a single region
			if lb == ub and lb%2 == 0:
				return coord2 - coord1 + 1*query_endPointInclusive
			if lb == ub and lb%2 == 1:
				return 0

			#print [coord1,coord2],(lb,ub)
			size_list   = [self.all_tracks[myChr][lb]-coord1]
			in_out_list = [(lb%2 == 0)]
			for i in xrange(lb+1,ub):
				size_list.append(self.all_tracks[myChr][i]-self.all_tracks[myChr][i-1])
				in_out_list.append((i%2 == 0))
			size_list.append(coord2-self.all_tracks[myChr][ub-1])
			in_out_list.append((ub%2 == 0))

			for i in xrange(len(size_list)):
				if in_out_list[i]:
					count += size_list[i] + 1*query_endPointInclusive

		return count