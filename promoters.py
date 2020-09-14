#!/usr/bin/env python3

import argparse
import pathlib
import sys
import os

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



parser = argparse.ArgumentParser(description='JGI promoter processing')
parser.add_argument('--min_genes', type=int, default=4,
	metavar='<int>', help='minimum genes per cluster [%(default)i]')
parser.add_argument('--upstream', type=int, default=1000,
	metavar='<int>', help='length of promoter region [%(default)i]')
parser.add_argument('--testing', type=int, default=0,
	metavar='<int>', help='short-circuit for testing [%(default)i]')
parser.add_argument('--summary', type=str, default='summary.txt',
	metavar='<path>', help='text file of information [%(default)s]')
arg = parser.parse_args()


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
summary = []
with open('public_species.txt') as fp:
	for line in fp.readlines():
		abbr, genus, species = line.split('\t')
		pathlib.Path(f'Promoters/{abbr}').mkdir(parents=True, exist_ok=True)
		
		clusters = get_clusters(abbr)
		for cid in clusters:
			if len(clusters[cid]['genes']) < arg.min_genes: continue
			with open(f'Promoters/{abbr}/{cid}.fa', 'w') as fp:
				seq = get_seq(abbr, clusters[cid]['scaff'])
				for gid in clusters[cid]['genes']:
					(beg, end, strand) = genes[abbr][gid]
					
					c1, c2 = None, None
					if strand == '+':
						c1 = beg - arg.upstream
						if c1 < 0: c1 = 0
						c2 = beg
					else:
						c1 = end
						c2 = end + arg.upstream
						if c2 > len(seq): c2 = len(seq)
					
					fp.write(f'>{abbr}:{cid}:{gid} {arg.upstream} {strand}\n')
					fp.write(f'{seq[c1:c2]}\n')
			
			size = clusters[cid]['size']
			desc = clusters[cid]['desc']
			scaf = clusters[cid]['scaff']
			
			summary.append((abbr, str(cid), str(size), desc, scaf))

		n += 1
		if arg.testing > 0 and n >= arg.testing: break

with open(arg.summary, 'w') as fp:
	fp.write(f'Species\tClusterID\tGenes\tType\tLocation\n')
	for cluster in summary:
		fp.write(('\t'.join(cluster)))
		fp.write('\n')

