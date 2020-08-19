import random
import math

def read_JASPAR(jasparfile):
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
	fp = 0
	fl = 0
	posdis = 0
	if jpos == '':
		fp += 1
		fl += 1
		posdis = 100
	elif jpos != '':
		mpos = int(mpos)
		jpos = int(jpos)
		mend = mpos + mw
		jend = jpos + jw
		if mend >= jpos and mpos < jend :
			posdis = mpos-jpos
		elif jend >= mpos and jpos < mend:
			posdis = mpos-jpos
		else:
			fl += 1	
			posdis = 98
	return fp, posdis, fl
	
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

def global_motcompare(motif1,motif2,background): #background info needs to be a dictionary
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

	'''
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
	
	
	
	#assert(local_motcompare(m1, m1) == 2) # maximum value or 1?
	#assert(local_motcompare(m1, m2) == 0) # minimum valu
	'''
	dl11 = local_motcompare(m2, m2)
	dl12 = local_motcompare(m1, m1)
	dl21 = local_motcompare(m2, m3)
	dl22 = local_motcompare(m3, m3)
	print(dl11, dl12, dl21, dl22)
	'''
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
	'''
	dl11 = compare_motifs(m1, m1)
	dl12 = local_motcompare(m1, m2)
	dl21 = local_motcompare(m2, m1)
	dl22 = local_motcompare(m2, m2)
	print(dl11, dl12, dl21, dl22)

	site = []
	'''
	

		
		
		
			
		
			
			
#test code	
#motif1 = read_JASPAR('Jaspar/MA0265.1.jaspar')
#motif2 = read_JASPAR('Jaspar/MA0266.1.jaspar')
#print(motifcompare(motif1,motif2))
#print(globalcompare(motif1,motif2))

	
	
	