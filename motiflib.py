import random

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
	
def pos_accuracy(mpos,jpos,mw,jw):
	fn = 0
	fp = 0
	fl = 0
	posdis = 0
	if jpos == '' and mpos != 'NA':
		fp += 1
		fl += 1
		posdis = 100
	elif jpos != '' and mpos == 'NA':
		fn += 1
		fl += 1
		posdis = 99
	elif jpos != '' and mpos!= 'NA':
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
	return mpos,jpos, fn, fp, posdis, fl
	
#uses manhattan distance (edit distance) to compare two pwms
#delete this function?
		
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
	score = 0 ### fix placement?
	#is this the best way?
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

if __name__ == '__main__':
	
	bkgd = {'A': 0.25, 'C': 0.25, 'G':0.25, 'T':0.25}
	motifs = [
		[],
		[{'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25}],
		[{'A':0.50, 'C':0.50, 'G':0.00, 'T':0.00}],
		[{'A':1.00, 'C':0.00, 'G':0.00, 'T':0.00}],
		[
			{'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
			{'A':1.00, 'C':0.00, 'G':0.00, 'T':0.00}
		],
	]

	for i in range(len(motifs)):
		for j in range(len(motifs)):
			lij = local_motcompare(motifs[i], motifs[j])
			lji = local_motcompare(motifs[j], motifs[i])
			gij = global_motcompare(motifs[i], motifs[j], bkgd)
			gji = global_motcompare(motifs[i], motifs[j], bkgd)
			assert(lij == lji)
			assert(gij == gji)
			print(i, j, lij, gij)
	
	test1 = pos_accuracy(16,16,16,20)
	print(test1)

	
	