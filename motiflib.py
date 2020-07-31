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
	
def calcmotif(sites):
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
	
def positionalaccuracy(mpos,jpos,mw,jw):
	fn = 0
	fp = 0
	fl = 0
	posdis = 0
	if jpos == '' and mpos != 'NA':
		fp += 1
		posdis = 100
	elif jpos != '' and mpos == 'NA':
		fn += 1
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
def motifcompare(motif1, motif2):
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

def globalcompare(motif1,motif2): #how should input background info?
	background = {'A': 0.25, 'C': 0.25, 'G':0.25, 'T':0.25}
	distances = []
	if len(motif1) > len(motif2):
		max = motif1
		min = motif2
	elif len(motif2) > len(motif1):
		max = motif2
		min = motif1
	for i in range(0,len(max)-len(min)+1):
		minmotif = []
		d = 0
		for j in range(0,len(max)):
			minmotif.append(background)
			#print(minmotif)
		for k in range(0,len(min)):
			minmotif[i+k] = min[k]
		#print(minmotif)
		for l in range(len(minmotif)):
			for nt in minmotif[l]:
				#manhattan similarity: allows for highest score
				d += 1-abs(minmotif[l][nt]-max[l][nt])
		distances.append(d)
	#print(distances)
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
		return bestfit
		
def localcompare(motif1, motif2):
	distances = []
	if len(motif1) > len(motif2):
		max = motif1
		min = motif2
	elif length(motif2) > len(motif1):
		max = motif2
		min = motif1
	for i in range(0,len(max)-len(min)):
		d = 0
		window = max[i:i+len(min)]
		for j in range(len(window)):
			for nt in window[j]:
				d += 1 - abs(window[j][nt]-min[j][nt])
		distances.append(d)
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
	return bestfit



		
		



		
		
		
		
			
		
			
			
#test code	
motif1 = read_JASPAR('Jaspar/MA0265.1.jaspar')
motif2 = read_JASPAR('Jaspar/MA0266.1.jaspar')
#print(motifcompare(motif1,motif2))
print(globalcompare(motif1,motif2))

	
	
	