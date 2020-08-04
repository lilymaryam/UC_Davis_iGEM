#!/usr/bin/env python3

def readmemeout(memeoutxml): #'meme_out/meme.xml'
	memeoutxml = 'meme_out/meme.xml' #delete later
	motifname = ''
	consensusseq = ''
	w = 0
	sites = 0 
	e_value = 0
	seq = ''
	numsites = 0
	motifid = ''
	strand = ''
	posi = 0
	background = 0
	result = []
	with open(memeoutxml) as mm:# how to do this ?
		for line in mm.readlines():
			if line.startswith('<motif id'):
				line = line.split(' ')
				for i in range(len(line)):
					item = line[i].split('="')
					if i == 1: motifname = item[1][:-1]
					if i == 2: consensusseq = item[1][:-1]
					if i == 4: w = int(item[1][:-1])
					if i == 5: sites = int(item[1][:-1])
					if i == 10: e_value = float(item[1][:-1])
			elif line.startswith('<scanned_sites '): 
				line = line.split(' ')
				for i in range(len(line)):
					item = line[i].split('="')
					if i == 1: seq = item[1][:-1]
					if i == 3: numsites = int(item[1][0]) #how should i fix this?
					if numsites == 0:
						motifid = 'NA'
						strand =  'NA'
						posi = 'NA'
					elif numsites > 0:
						if i == 4: motifid = item[1][:-1]
						if i == 5: strand = item[1][:-1]
						if i == 6: posi = int(item[1][:-1])		
				if motifid == 'NA': result.append((seq, numsites, motifid, strand, posi))
				else:               result.append((seq, numsites, motifid, strand, posi,\
				consensusseq, w, sites, e_value))
	return result


def memepwm(memeouttxt): #'meme_out/meme.txt'
	mememotif = []
	nt = ['A','C','G','T']
	with open(memeouttxt) as mt:
		line = mt.readline()
		while line:
			if line.startswith('letter-probability matrix'): 
				line = line.split()
				wid = int(line[5])
				for i in range(wid):
					mememotif.append({})
					line = mt.readline()
					line = line.split()
					for j in range(len(nt)):
						mememotif[i][nt[j]] = float(line[j])
			line = mt.readline()
	return mememotif

def read_testmotif(motif_file):
	jpositions = []
	with open(motif_file) as tm:
		for line in tm.readlines():
			if line.startswith('>'): 
				line = line.split()
				jposi = line[1].strip('[]')
				jposi = jposi.strip("''")
				jposi = jposi.strip('+')
				jpositions.append(jposi)
	return jpositions
	

	
#if __name__ == '__main__':
#	inf = readtest			
#print(memepwm('meme_out/meme.txt'))						

#result = (readmemeout('meme_out/meme.xml'))		
#for i in range(0,len(result)):
#	print(result[i])