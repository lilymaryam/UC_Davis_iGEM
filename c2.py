

NT = ['A', 'C', 'G', 'T']

IUPAC = {
	'A': {'A':1, 'C':0, 'G':0, 'T':0},
	'C': {'A':0, 'C':1, 'G':0, 'T':0},
	'G': {'A':0, 'C':0, 'G':1, 'T':0},
	'T': {'A':0, 'C':0, 'G':0, 'T':1},
	'R': {'A':0.5, 'C':0.0, 'G':0.5, 'T':0.0},
	'Y': {'A':0.0, 'C':0.5, 'G':0.0, 'T':0.5},
	'M': {'A':0.5, 'C':0.5, 'G':0.0, 'T':0.0},
	'K': {'A':0.0, 'C':0.0, 'G':0.5, 'T':0.5},
	'W': {'A':0.5, 'C':0.0, 'G':0.0, 'T':0.5},
	'S': {'A':0.0, 'C':0.5, 'G':0.5, 'T':0.0},
	'B': {'A':0.01, 'C':0.33, 'G':0.33, 'T':0.33},
	'D': {'A':0.33, 'C':0.01, 'G':0.33, 'T':0.33},
	'H': {'A':0.33, 'C':0.33, 'G':0.01, 'T':0.33},
	'V': {'A':0.33, 'C':0.33, 'G':0.33, 'T':0.01},
	'N': {'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
}

def manhattan(d1, d2):
	d = 0
	for nt in NT:
		d += abs(d1[nt] - d2[nt])
	return d

def motif2seq(m):
	seq = []
	for i in range(len(m)):
		min = 1e9
		best = None
		for s in IUPAC:
			d = manhattan(m[i], IUPAC[s])
			if d < min:
				min = d
				best = s
		seq.append(best)
	return ''.join(seq)
	
def anti_motif(m):
	a = []
	for i in range(len(m)):
		d = {}
		d['A'] = m[i]['T']
		d['C'] = m[i]['G']
		d['G'] = m[i]['C']
		d['T'] = m[i]['A']
		a.append(d)
	reversed(a)
	return a

def motif_similarity(m1, m2):
	m = []
	for i in range(len(m1)):
		m.append([])
		for j in range(len(m2)): m[i].append([])
	
	# init first row and column
	for i in range(len(m1)): m[i][0] = 2 - manhattan(m1[i], m2[0])
	for j in range(len(m2)): m[0][j] = 2 -manhattan(m1[0], m2[j])

	# compute diagonals and find best score
	max_score = 0	
	for i in range(1, len(m1)):
		for j in range(1, len(m2)):
			m[i][j] = m[i-1][j-1] + 2 - manhattan(m1[i], m2[j])
			if m[i][j] > max_score:
				max_score = m[i][j]
	
	return max_score



m0 = [
	{'A':1, 'C':0, 'G':0, 'T':0},
	{'A':0, 'C':1, 'G':0, 'T':0},
	{'A':0, 'C':0, 'G':1, 'T':0},
	{'A':0, 'C':0, 'G':0, 'T':1},
	{'A':0.5, 'C':0.0, 'G':0.5, 'T':0.0},
	{'A':0.0, 'C':0.5, 'G':0.0, 'T':0.5},
	{'A':0.5, 'C':0.5, 'G':0.0, 'T':0.0},
	{'A':0.0, 'C':0.0, 'G':0.5, 'T':0.5},
	{'A':0.5, 'C':0.0, 'G':0.0, 'T':0.5},
	{'A':0.0, 'C':0.5, 'G':0.5, 'T':0.0},
	{'A':0.01, 'C':0.33, 'G':0.33, 'T':0.33},
	{'A':0.33, 'C':0.01, 'G':0.33, 'T':0.33},
	{'A':0.33, 'C':0.33, 'G':0.01, 'T':0.33},
	{'A':0.33, 'C':0.33, 'G':0.33, 'T':0.01},
	{'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25},
]

print(motif2seq(m0))

m1 = [
	{'A':1.0, 'C': 0.0, 'G': 0.0, 'T': 0.0},
	{'A':0.0, 'C': 1.0, 'G': 0.0, 'T': 0.0},
	{'A':0.0, 'C': 0.0, 'G': 1.0, 'T': 0.0},
	{'A':0.0, 'C': 0.0, 'G': 0.5, 'T': 0.5},
	{'A':0.0, 'C': 0.0, 'G': 0.5, 'T': 0.5},
	{'A':0.0, 'C': 0.0, 'G': 0.5, 'T': 0.5},
]

m2 = [
	{'A':0.0, 'C': 1.0, 'G': 0.0, 'T': 0.0},
	{'A':0.0, 'C': 0.0, 'G': 1.0, 'T': 0.0},
	{'A':0.0, 'C': 0.0, 'G': 0.2, 'T': 0.8},
]

print(motif_similarity(m1, m2))
print(motif_similarity(m1, anti_motif(m2)))
print(motif_similarity(m0, m1))


