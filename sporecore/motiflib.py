#!/usr/bin/env python3
import random
import math



'''
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

def read_memetxt(memetxt):
	'''Opens text file output of meme and extracts the sites discovered by meme'''
	sites = []
	motif_stats = []
	with open(memetxt) as mt:
		while True:
			line = mt.readline()
			if line == '': break
			if line.startswith('MOTIF'):
				l = line.split()
				motifid = l[2]
				nsites = l[8]
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
					sites.append((motifid,nsites,strand,seq,beg,pval))
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
	'''Opens text file output of meme and extracts the position weight matrix
	 and it's corresponding information for every motif discovered by meme'''
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
					nsites = int(l[7])
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
					bits,p_bits = score_motifbit(m)	
					motif_stats.append((f'MEME-{nummot}',wid,nsites,e_val,bits,p_bits))	
					nummot += 1
	return motifs, motif_stats

#are variable names too specific?
def read_testmotif(motif_file):
	'''Opens the fasta file generated  '''
	jpositions = []
	numjsites = 0
	#numfsites = 0
	with open(motif_file) as tm:
		for line in tm.readlines():
			if line.startswith('>'):
				positions = []
				#numpos = 0 
				line = line.split()
				#print(line)
				numjsites += (int((len(line[1:(len(line))]))/2))
				seq  = line[0].strip('>')
				#jpositions.append(seq)
				for i in range(1,len(line)):
					#print(i)
					#print(i,line)
					#print('i',i,'len(line)',len(line))
					#print(line[i])
					if i == 1: 
						if line[1] == '[]':
							line[1] = line[1].strip("]")
						positions.append(line[1].strip("['"))
					#print(positions)
					elif i == len(line)-1:
						#print('yes!!!!!!!!!!!!!!!!!')
						#print('i+1',i,'len line', len(line))
						positions.append(line[i].strip("]'"))
					else:
						#jposi = int(jposi)
						#print(positions
						positions.append(line[i].strip("',"))
						#print('strand',line[i+1].strip("',"))
						#print('strand',line[i+1].strip("',"))
					#strand = ''
				#print(seq,positions)
				jpositions.append((seq,positions))
		#print(jpositions)
		#print(numjsites)
		#seqnum += 1
		#print(positions)
		#print(positions)
		#jpositions.append(positions)
	#print(len(jpositions))
	#print(jpositions)
	return jpositions, numjsites
	'''
				for i in range(1,len(line),2):
					print('i',i,'len(line)',len(line))
					print(line[i])
					if i == 1: #positions.append
						print('first',line[1].strip("['"))
					#print(positions)
					elif i+1 == len(line)-1:
						#print('yes!!!!!!!!!!!!!!!!!')
						print('i+1',i,'len line', len(line))
						#positions.append
						print('last',line[i+1].strip("]'"))
					#jposi = int(jposi)
					#print(positions
					print('pos',line[i].strip("',"))
					print('strand',line[i+1].strip("',"))
						#print('strand',line[i+1].strip("',"))
					#strand = ''
				#jpositions.append((f'seq-{seqnum}',jposi,strand))
				#seqnum += 1
	#return jpositions
	'''
	


def read_JASPAR(jasparfile):
	'''Opens a jaspar file containing the binding motif of interest and creates 
	a position weight matrix representing the binding motif'''
	motif = []
	linenum = 0
	nt = ['A', 'C', 'G', 'T']
	with open(jasparfile) as jf:
		for line in jf.readlines():
			if line.startswith('>'): continue
			line = line.split()	
			for i in range(2,len(line)-1):
				if len(motif) < (len(line)-3): motif.append({})
				motif[i-2][nt[linenum]] = int(line[i])
			linenum+= 1 
	for i in range(len(motif)):
		total = 0 
		for nt in motif[i]: total += motif[i][nt]
		for nt in motif[i]: motif[i][nt] /= total
	return motif 

def generate_site(motif):
	'''Generates a single DNA sequence from the position weight matrix '''
	motifseq = []
	total = 0
	for i in range(0,len(motif)):
		pA = motif[i]['A']
		pC = motif[i]['C']
		pG = motif[i]['G']
		pT = motif[i]['T']
		r = random.random()
		if r < pA:               motifseq += 'A'
		elif r < (pA + pC):      motifseq += 'C'
		elif r < (pA + pC + pG): motifseq += 'G' 
		else:                    motifseq += 'T'
	return motifseq
	
def new_pwm(sites):
	'''Accepts an array of binding site sequences and produces a position 
	weight matrix from them '''
	newmot =[]
	for i in range(len(sites[0])):
		A = 0
		C = 0
		G = 0
		T = 0
		newmot.append({})
		for j in range(len(sites)):
			s = sites[j][i]
			if s == 'A': A += 1
			elif s == 'C': C += 1
			elif s == 'G': G += 1
			elif s == 'T': T += 1
		newmot[i]['A'] = A/len(sites)
		newmot[i]['C'] = C/len(sites)
		newmot[i]['G'] = G/len(sites)
		newmot[i]['T'] = T/len(sites)
	return newmot
	
#meme.txt files do not find false negatives 

def pos_accuracy(mpos,jpos,mw,jw):
	'''Accepts and compares the starting positions of the meme motifs and the
	 embedded jaspar motifs to determine how far apart they are and how much
	  overlap there is   '''
	fp = 0
	fl = 0
	posdis = 0
	if jpos == '':
		fp += 1
		fl += 1
		posdis = 100
		overlap = 0
		overlap_p = 0
	elif jpos != '':
		mpos = int(mpos)
		jpos = int(jpos)
		mend = mpos + mw
		jend = jpos + jw
		if mpos >= jpos and mpos <= jend:
			posdis = mpos-jpos
			if mend >= jend:
				overlap = jend - mpos
			else:
				overlap = mend - mpos
			overlap_p = overlap/jw
			if overlap_p < 0.2:
				fp += 1
				fl += 1		
		elif mpos <= jpos and jpos <= mend:
			posdis = jpos-mpos
			if mend <= jend:
				overlap = mend - jpos
			else:
				overlap = jend - jpos
			overlap_p = overlap/jw
			if overlap_p < 0.2:
				fp += 1
				fl += 1
		else:
			fl += 1	
			posdis = 98
			overlap = 0
			overlap_p = 0
	return fp, posdis, fl, overlap, overlap_p





'''	
#uses manhattan distance (edit distance) to compare two pwms
#delete this function?
#this function is superfluous
def compare_motifs(motif1, motif2):
	if len(motif1) == len(motif2):
		d = 0
		for i in range(len(motif1)):
			for nt in motif1[i]:
				d += abs(motif1[i][nt]-motif2[i][nt])
	#if the motifs aren't the same length, will assess windows of the larger motif
	#against the shorter motif
	else:
		d = 0
		if len(motif1) > len(motif2):
			max = motif1
			min = motif2
		if len(motif2) > len(motif1):
			max = motif2
			min = motif1
		distances = []
		for i in range(0,len(max)-len(min)):#need +1?
			d = 0
			window = max[i:i+len(min)]
			for j in range(len(window)):
				for nt in min[j]:
					d += abs(min[j][nt]-window[j][nt])
			distances.append(d)
		d = 0
		winpos = 0
		for w in range(len(distances)):
			if w == 0: d = distances[w]
			elif distances[w] <= d: 
				d = distances[w]
				winpos = w
	return d
'''

def global_motcompare(motif1,motif2,background): #background info needs to be a dictionary
	'''Finds global similarity score between the meme motif and the jaspar 
	motif by finding the best matching region between the two motifs and 
	filling unmatched positions with background nucleotide distribution 
	frequency  '''
	distances = []
	score = 0 
	if len(motif1)==0 or len(motif2)== 0:
		score = 0
		distances.append(score)
		max_score = 1
	elif len(motif1) == len(motif2):
		max_score = 2*len(motif1)
		for i in range(0, len(motif1)):
			d = 0
			for nt in motif1[i]:
				d += abs(motif1[i][nt] - motif2[i][nt])
			score += 2-d
		distances.append(score)		
	else:
		if len(motif1) > len(motif2):
			big = motif1
			short = motif2
		elif len(motif2) > len(motif1):
			big = motif2
			short = motif1
		max_score = 2*len(big)
		for i in range(0,len(big)-len(short)+1):
			minmotif = []
			for j in range(0,len(big)):
				minmotif.append(background)
			for k in range(0,len(short)):
				minmotif[i+k] = short[k]
			score = 0
			for l in range(len(minmotif)):
				d = 0
				for nt in minmotif[l]:
					#manhattan similarity: allows for highest score
					d += abs(minmotif[l][nt]-big[l][nt])
				score += 2-d
			distances.append(score)
	bestfit = 0
	fitindex = 0
	for i in range(0,len(distances)):
		if i == 0: 
			bestfit = distances[i]
			fitindex = i
		else:
			if distances[i] >= bestfit:
				bestfit = distances[i]
				fitindex = i
	return bestfit, bestfit/max_score
		
def local_motcompare(motif1, motif2):
	'''Finds local similarity score between the meme motif and the jaspar 
	motif by finding the best matching region between the two motifs.'''
	distances = []
	if len(motif1) == 0 or len(motif2)==0:
		score = 0
		distances.append(score)
		max_score = 1
	elif len(motif1) == len(motif2):
		max_score = 2*len(motif1)
		score = 0
		for i in range(len(motif1)):
			d = 0
			for nt in motif1[i]:
				d += abs(motif1[i][nt]-motif2[i][nt])
			score += 2-d
		distances.append(score)
	else:
		if len(motif1) > len(motif2):
			big = motif1
			short = motif2
		elif len(motif2) > len(motif1):
			big = motif2
			short = motif1
		max_score = 2*len(big)
		for i in range(0,len(big)-len(short)+1):
			score = 0
			window = big[i:i+len(short)]
			for j in range(len(window)):
				d = 0
				for nt in window[j]:
					d += abs(window[j][nt]-short[j][nt])
				score += 2-d
			distances.append(score)
	bestfit = 0
	fitindex = 0
	for i in range(len(distances)):
		if i == 0: 
			bestfit = distances[i]
			fitindex = i
		else: 
			if distances[i] > bestfit:
				bestfit = distances[i]
				fitindex = i
	return bestfit, bestfit/max_score
	
def score_motifbit(motif):
	'''Finds the informational content of the motif by summing up the information content of each individual position  '''
	score = 0
	max_score = 2*len(motif)
	for i in range(len(motif)):
		entropy = 0
		pseudos = 0
		for nt in motif[i]:
			prob = motif[i][nt]
			if prob == 0:
				continue
			entropy -= prob* math.log2(prob)
		bits = 2 - entropy - pseudos*.1
		score += bits
	if max_score > 0:
		p = score/max_score
	else:
		p = ''
	return score, p

	
	
	
if __name__ == '__main__':
	'''
	m = [
	{'A':.25, 'C':.25, 'G':.25, 'T':.25},
	{'A':1, 'C':0, 'G':0, 'T':0},
	{'A':.75, 'C':0, 'G':.25, 'T':0}
	]
	m1 = []
	m2 = [
	{'A':1, 'C':0, 'G':0, 'T':0},
	{'A':1, 'C':0, 'G':0, 'T':0},
	{'A':1, 'C':0, 'G':0, 'T':0}
	
	]
	m3 = [
	{'A':.25, 'C':.25, 'G':.25, 'T':.25},
	{'A':.25, 'C':.25, 'G':.25, 'T':.25},
	{'A':.25, 'C':.25, 'G':.25, 'T':.25},
	]
	print(score_motifbit(m))
	print(score_motifbit(m1))
	print(score_motifbit(m2))
	print(score_motifbit(m3))

	
	m1 = []
	m2 = [
		{'A':0, 'C':0, 'G':0, 'T':1},
		{'A':.22, 'C':.28, 'G':0, 'T':.5}
		]
	m3 = [
		{'A':.25, 'C':.25, 'G':.5, 'T':0},
		{'A':1, 'C':0, 'G':0, 'T':0},
		{'A':0, 'C':1, 'G':0, 'T':0},
		
		]
	

	print(read_testmotif('testmotif0265.fa'))
	#print(read_memetxt('meme_out/meme.txt'))
	#print(memepwm('meme_out/meme.txt'))
	#print(read_testmotif('testmotif0265.fa'))
	
	#assert(local_motcompare(m1, m1) == 2) # maximum value or 1?
	#assert(local_motcompare(m1, m2) == 0) # minimum valu
	
	dl11 = local_motcompare(m2, m2)
	dl12 = local_motcompare(m1, m1)
	dl21 = local_motcompare(m2, m3)
	dl22 = local_motcompare(m3, m3)
	print(dl11, dl12, dl21, dl22)
	
	dg11 = global_motcompare(m1, m3, {'A': 0.25, 'C': 0.25, 'G':0.25, 'T':0.25})
	dg12 = global_motcompare(m1, m2, {'A': 0.25, 'C': 0.25, 'G':0.25, 'T':0.25})
	dg21 = global_motcompare(m2, m3, {'A': 0.25, 'C': 0.25, 'G':0.25, 'T':0.25})
	dg22 = global_motcompare(m2, m2, {'A': 0.25, 'C': 0.25, 'G':0.25, 'T':0.25})
	dg33 = global_motcompare(m3, m3, {'A': 0.25, 'C': 0.25, 'G':0.25, 'T':0.25})
	print(dg11)
	print(dg12)
	print(dg21)
	print(dg22)
	print(dg33)
	
	dl11 = compare_motifs(m1, m1)
	dl12 = local_motcompare(m1, m2)
	dl21 = local_motcompare(m2, m1)
	dl22 = local_motcompare(m2, m2)
	print(dl11, dl12, dl21, dl22)

	site = []
	'''
	#mpos = 149
	#mw = 50
	#jpos = 120
	#jw = 50
	#fp, posdis, fl, overlap = (pos_accuracy(mpos,jpos,mw,jw))
	#print('fp',fp)
	#print('pd',posdis)
	#print('fl',fl)
	#print('ol',overlap)
	
	#print(pos_accuracy(mpos,jpos,mw,jw))
	#print('1',read_testmotif('BS667tf1.fa'))
	#print('2',read_testmotif('BS667tf2.fa'))
	#print(read_testmotif('667_25.fa'))
	
		
			
		
			
			
#test code	
#motif1 = read_JASPAR('Jaspar/MA0265.1.jaspar')
#motif2 = read_JASPAR('Jaspar/MA0266.1.jaspar')
#print(motifcompare(motif1,motif2))
#print(globalcompare(motif1,motif2))

	
	
	