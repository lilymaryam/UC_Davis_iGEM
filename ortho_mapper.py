import json
import os
import sys
import math

def ortholog(file, parent):
	match = {}
	with open(file) as fp:
		for line in fp.readlines():
			if line.startswith('#'): continue
			f = line.split()
			query = f[0]
			sbjct = f[1]
			score = int(f[5])
			s1, c1, p1 = query.split(':')
			p2, s2 = sbjct.split('|')
			if p2 not in parent: continue
			
			c2 = parent[p2]
			if c2 not in match: match[c2] = {}
			if p1 not in match[c2]: match[c2][p1] = score
	
	for c2 in match:
		stuff = []
		for p1 in match[c2]:
			stuff.append(match[c2][p1])
		
	# assign matching cluster based on max count then max score
	max_count = 0
	max_score = 0
	cluster = None
	for c2 in match:
		count = len(match[c2])
		score = 0
		for p1 in match[c2]: score += match[c2][p1]
		if count > max_count:
			max_count = count
			max_score = score
			cluster = c2
		elif count == max_count:
			if score > max_score:
				max_count = count
				max_score = score
				cluster = c2
	
	return cluster, max_count, max_score

def score(n, k):
	if k > n: k = n
	if n == 0: return 0
	return (k / n) * math.log2(n)

data = None
with open('clusters.json') as fp:
	data = json.load(fp)

for s1 in data['species']:
	with open(f'ORTHOLOGS/{s1}.txt', 'w') as fp:
		for c1 in data['species'][s1]:
			for s2 in data['species']:
				if s1 == s2: continue
				file = f'BLASTOUT/{s1}-{c1}-{s2}.blastp'
				c2, k, s = ortholog(file, data['parent_cluster'])
				l1 = 0
				if c1 not in data['clusters']: l1 = 0
				else: l1 = len(data['clusters'][c1])
				l2 = 0
				if c2 not in data['clusters']: l2 = 0
				else: l2 = len(data['clusters'][c2])
				
				s = score(l1, k) + score(l2, k)
				fp.write(f'{s1}:{c1} {s2}:{c2} {k} {l1} {l2} {s:.3f}\n')


