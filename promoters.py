#!/usr/bin/env python3

import sys
import os
import pathlib
import json

MIN_GENES_PER_CLUSTER = 4
UPSTREAM = 100

def get_clusters(abbr):
	file = f'BGC/{abbr}__BGC.txt'
	cluster = {}
	with open(file) as fp:
		id = None
		while True:
			line = fp.readline()
			if line.startswith('#'): break
		while True:
			line = fp.readline()
			if line.startswith('#'): break
			if line.startswith('Cluster'):
				f = line.split()
				size = int(f[2])
				gene = int(f[6])
				scaff = f[7]
				desc = f[10]
				id = gene
				cluster[id] = {
					'size':size,
					'desc':desc,
					'scaff':scaff,
					'genes':[]
				}
			else:
				f = line.split()
				cluster[id]['genes'].append(int(f[1]))
	return cluster

def get_seq(abbr, scaff):
	file = f'FASTA/{abbr}.fasta'
	search = f'>{scaff}'
	with open(file) as fp:
		while True:
			line = fp.readline()
			if line.startswith(search):
				seqs = []
				while True:
					line = fp.readline()
					if line.startswith('>'): break
					if line == '': break
					seqs.append(line.rstrip())
				return ''.join(seqs)
	print('whoops')
	sys.exit(1)

genes = {}
with open('gene_coordinates.txt') as fp:
	for line in fp.readlines():
		f = line.split()
		abbr = f[0]
		gid = int(f[1])
		beg = int(f[4])
		end = int(f[5])
		strand = f[3]
		if abbr not in genes: genes[abbr] = {}
		genes[abbr][gid] = (beg, end, strand)

n = 0
with open('public_species.txt') as fp:
	for line in fp.readlines():
		abbr, genus, species = line.split('\t')
		pathlib.Path(f'Promoters/{abbr}').mkdir(parents=True, exist_ok=True)
		
		clusters = get_clusters(abbr)
		for cid in clusters:
			if len(clusters[cid]['genes']) < MIN_GENES_PER_CLUSTER: continue
			with open(f'Promoters/{abbr}/{cid}.fa', 'w') as fp:
				seq = get_seq(abbr, clusters[cid]['scaff'])
				for gid in clusters[cid]['genes']:
					(beg, end, strand) = genes[abbr][gid]
					
					c1, c2 = None, None
					if strand == '+':
						c1 = beg - UPSTREAM
						if c1 < 0: c1 = 0
						c2 = beg
					else:
						c1 = end
						c2 = end + UPSTREAM
						if c2 > len(seq): c2 = len(seq)
					
					fp.write(f'>{abbr}:{cid}:{gid} {UPSTREAM} {strand}\n')
					fp.write(f'{seq[c1:c2]}\n')

		n += 1
		if n > 6: sys.exit() # remove this for full processing
