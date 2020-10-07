#!/usr/bin/env python3

import argparse
import random
import motiflib
import math
import os
import make_backgroundseq



extended_help = """
%(prog)s is a program that simulates promoter regions by embedding binding 
sites in DNA sequences modeled from real promoter regions. The sites are 
generated from a position weight matrix representing the frequency distribution
of nucleotides for each position of the motif. The position weight matrix is 
read from a JASPAR formatted file which is a command line argument.
"""

parser = argparse.ArgumentParser(
	description='Embeds instances of a motif into DNA sequences of markov order k',
	epilog=extended_help)
parser.add_argument('--jasparfile', required=True, type=str,
	metavar='<str>', help='jaspar file')
parser.add_argument('--dnafile', required=True, type=str,
	metavar='<str>', help='fasta file of promoters')
parser.add_argument('--markov_order', required=True, type=int,
	metavar='<int>', help='markov order of generated sequence')
parser.add_argument('--numseq', required=False, type=int, default=10,
	metavar='<int>', help='number of sequences to generate [%(default)i]')
parser.add_argument('--seqlen', required=False, type=int, default=100,
	metavar='<int>', help='length of sequences to generate [%(default)i]')
parser.add_argument('--mps', required=False, type=int, default=1,
	metavar='<int>', help='number of motifs per seqeunce [%(default)i]')
parser.add_argument('--freq', required=False, type=float, default=0.9,
	metavar='<float>', help='optional floating point argument [%(default).3f]')	
parser.add_argument('--bothstrands', action='store_true',
	help='on/off switch')
parser.add_argument('--negstrands', action='store_true',
	help='on/off switch')

	
arg = parser.parse_args()

#if arg.bothstrands and arg.negstrands:
	

assert(arg.freq <= 1.0)
if arg.bothstrands and arg.negstrands:
	raise ValueError('Strandedness cannot be both and negative together')


def generate_markovseq(numseq,seqlen,data_file,k):
	dna = make_backgroundseq.organize_dna(arg.dnafile)
	kmers = make_backgroundseq.make_kmerdict(arg.markov_order,dna)
	pool = make_backgroundseq.make_contextpool(kmers)
	seq = make_backgroundseq.generate_seq(arg.markov_order,pool,arg.seqlen)
	seq = seq.lower()
	list_seq = []
	for i in range(len(seq)):
		list_seq.append(seq[i])
	return list_seq


motif = motiflib.read_JASPAR(arg.jasparfile)
assert(len(motif) <= arg.seqlen)
for i in range(0, arg.numseq):
	seq = (generate_markovseq(arg.numseq,arg.seqlen,arg.dnafile,\
	arg.markov_order))
	r = random.random()
	places = []
	if r < arg.freq:
		for j in range(0,arg.mps):
			strand="+"
			site = motiflib.generate_site(motif)
			assert(len(motif)== len(site))
			place = random.randint(0, len(seq)-len(motif))
			if arg.bothstrands or arg.negstrands:
				r = random.random()
				if r < 0.5 and arg.bothstrands:
					strand = '-'
				if arg.negstrands:
					strand = '-'
				if strand == '-':
					site.reverse()
					for l in range(len(site)):
						if site[l] == 'A':
							site[l] = 'T'
						elif site[l] == 'C':
							site[l] = 'G'
						elif site[l] == 'G':
							site[l] = 'C'
						elif site[l] == 'T':
							site[l] = 'A'
			places.append(f'{place} {strand}')
			for k in range(0,len(site)):
				seq[place+k] = site[k]
	print(f'>seq-{i} {places}')
	print(''.join(seq))
