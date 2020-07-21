#!/usr/bin/env python3
"""
(Correct if incorrect)
This program can either take the output of motifapalooza.py or any FASTA file you want,
though we'd mostly be using this on groups of promoter regions where we want to find
the TF binding site. You can customize this with different parameters, described in the
argparse section (lines ~20-38). Things this program outputs: the motif that MEME founds, the p-value
of those motifs (i.e. how likely that these motifs arose by chance), and . . . other stuff?
Either this or an additional program will analyze the MEME-generated motifs and quantify
how "good" MEME was at finding the right motif, i.e. the motif that we either put in there
(in the case of the motifapalooza outputs) or that we know to be there (in the case of the
clusters with known TF binding sites inside).
"""
import argparse
import re
import os
import sys


parser = argparse.ArgumentParser(
	description='Runs MEME on fasta file and reports results.')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<str>', help='fasta file of promoters')
parser.add_argument('--model', required=False, type=str, default='zoops',
	metavar='<int>', help='model type {zoops, oops, anr} [%(default)s]')
parser.add_argument('--nmotifs', required=False, type=int, default=5,
	metavar='<int>', help='maximum number of motifs to find [%(default)i]')
parser.add_argument('--minw', required=False, type=int, default=6,
	metavar='<int>', help='minimum width [%(default)i]')
parser.add_argument('--maxw', required=False, type=int, default=18,
	metavar='<int>', help='maximum width [%(default)i]')
parser.add_argument('--minsites', required=False, type=int, default=0,
	metavar='<int>', help='minimum number of sites per motif [%(default)i]')
parser.add_argument('--maxsites', required=False, type=int, default=4,
	metavar='<int>', help='maximum number of sites per motif [%(default)i]')
parser.add_argument('--background', required=False, type=int, default=0,
	metavar='<int>', help='Markov model for background [%(default)i]')
arg = parser.parse_args()

opts = []
opts.append(f'-dna -revcomp')
opts.append(f'-mod {arg.model}')
opts.append(f'-nmotifs {arg.nmotifs}')
opts.append(f'-minw {arg.minw} -maxw {arg.maxw}')
opts.append(f'-minsites {arg.minsites} -maxsites {arg.maxsites}')
opts.append(f'-markov_order {arg.background}')
opts.append(f'-nostatus')
memeopts = ' '.join(opts)


os.system(f'meme {arg.fasta} {memeopts} 2>/dev/null')

found = {}
with open('meme_out/meme.txt') as fp:

	while True:
		line = fp.readline()
		if line == '': break
		if line.startswith('MOTIF'):
			f = line.split()
			name = f[1]
			sites = []
			
			# fast foward to sites sort by p-value
			while True:
				line = fp.readline()
				if line.startswith('Sequence name'):
					crap = fp.readline()
					break
			
			# read sites
			while True:
				line = fp.readline()
				if line.startswith('---'): break
				f = line.split()
				gene = f[0]
				strand = f[1]
				beg = f[2]
				pval = f[3]
				sites.append((gene, beg, strand, pval))
		
			# save
			found[name] = sites
			
for name in found:
	print(name)
	for sites in found[name]:
		print(f'	{"	".join(sites)}')



"""
alfR GGGWTCGAWCCC
"""