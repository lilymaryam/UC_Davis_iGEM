#!/usr/bin/env python3


def read_memetxt(memetxt):
	sites = []
	motif_stats = []
	with open(memetxt) as mt:
		while True:
			line = mt.readline()
			if line == '': break
			if line.startswith('MOTIF'):
				l = line.split()
				motifid = l[2]
				#skip to sequence information
				while True:
					line = mt.readline()
					if line.startswith('Sequence name'):
						line = mt.readline()
						break
				#parse through sequence information
				while True:
					line = mt.readline()
					if line.startswith('---'): break
					l = line.split()
					seq = l[0]
					if l[1] == '+' or l[1] == '-':
						strand = l[1]
						beg = int(l[2])
						pval = float(l[3])	
					else:
						strand = '+'
						beg = int(l[1])
						pval = float(l[2])
					sites.append((motifid,strand,seq,beg,pval))
				#could potentially read out info about the motifs
				while True:
					line = mt.readline()
					if line.startswith('letter-probability matrix'):
						l = line.split()
						wid = int(l[5])
						e_val = float(l[9])
						break						
		return sites

				
				
'''
#obsolete
def readmemeout(memeoutxml): #'meme_out/meme.xml'
	#memeoutxml = 'meme_out/meme.xml' #delete later
	seq = ''
	numsites = 0
	motifid = ''
	strand = ''
	posi = 0
	background = 0
	result = []
	motifs = []
	with open(memeoutxml) as mm:# how to do this ?
		for line in mm.readlines():
			if line.startswith('<motif id'):
				motifname = ''
				consensusseq = ''
				w = 0
				sites = 0 
				e_value = 0
				line = line.split(' ')
				for i in range(len(line)):
					item = line[i].split('="')
					if i == 1: motifname = item[1][:-1]
					if i == 2: consensusseq = item[1][:-1]
					if i == 4: w = int(item[1][:-1])
					if i == 5: sites = int(item[1][:-1])
					if i == 10: e_value = float(item[1][:-1])
				motifs.append((motifname, consensusseq, w, sites, e_value))
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
		print(motifs)
	return result
'''

def memepwm(memeouttxt): #'meme_out/meme.txt'
	motifs = []
	motif_stats = []
	nummot = 1
	nt = ['A','C','G','T']
	with open(memeouttxt) as mt:
		#iterates through the whole file
		while True:
			line = mt.readline()
			if line == '': break
			#finds the beginning of the pwms
			while True:
				if line == '': break
				line = mt.readline()
				if line.startswith('letter-probability matrix'):
					m = []
					l = line.split()
					wid = int(l[5])
					e_val = float(l[9])
					for i in range(wid):
						m.append({})
						#iterates through the pwms
					for i in range(len(m)):
						line = mt.readline()
						l = line.split()
						for n in range(len(nt)):
							m[i][nt[n]] = float(l[n])
					motifs.append(m)	
					motif_stats.append((f'MEME-{nummot}',wid,e_val))	
					nummot += 1
	return motifs, motif_stats

#are variable names too specific?
def read_testmotif(motif_file):
	jpositions = []
	seqnum = 0
	with open(motif_file) as tm:
		for line in tm.readlines():
			if line.startswith('>'): 
				line = line.split()
				jposi = line[1].strip("['")
				if len(line) == 3:
					strand = line[2].strip("]'")
					jposi = int(jposi)
				else:
					jposi = jposi.strip("]")
					strand = ''
				jpositions.append((f'seq-{seqnum}',jposi,strand))
				seqnum += 1
	return jpositions
	

	
if __name__ == '__main__':
	print(read_testmotif('testmotif0265.fa'))
	#print(read_memetxt('meme_out/meme.txt'))
	#print(memepwm('meme_out/meme.txt'))
	#print(read_testmotif('testmotif0265.fa'))
	