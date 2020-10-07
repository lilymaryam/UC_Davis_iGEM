#!/usr/bin/python3

import argparse
import random 


def organize_dna(dnafile):
	dna = ''
	with open(dnafile) as df:
		for line in df.readlines():
			line = line.strip()
			if not line.startswith('>'):
				dna += line
	return dna


def make_kmerdict(k, dna):
	kmers = {}
	for i in range(len(dna)-k):
		ctx = dna[i:i+k]
		nt = dna[i+k]
		if ctx not in kmers:
			kmers[ctx] = {}
		if nt not in kmers[ctx]:
			kmers[ctx][nt] = 0 
		kmers[ctx][nt] += 1
	if len(kmers) < 4**k:
		raise ValueError('Not enough DNA sequence for markov order')
	else:
		return kmers



def make_contextpool(kmers):
	pool = {}
	for ctx in kmers:
		pool[ctx] = '' 
		for nt in kmers[ctx]:
			pool[ctx] = pool[ctx] + kmers[ctx][nt]*nt
	return pool

def generate_seq(k, pool,seqsize):
	seq = ''
	for i in range(k):
		seq += random.choice('ACGT')
	for i in range(k,seqsize):
		prev = seq[i-k:i]
		if prev in pool:
			seq += random.choice(pool[prev])
	return seq


if __name__ =='__main__':
	parser = argparse.ArgumentParser(
	description='Generates background DNA sequence of Markov Order k.')
	# required arguments
	parser.add_argument('--dnafile', required=True, type=str,
		metavar='<str>', help='fastafile containing promoter DNA')
	parser.add_argument('--markov_order', required=True, type=int,
		metavar='<int>', help='What order markov model')
	parser.add_argument('--seqlen', required=True, type=int,
		metavar='<int>', help='Size of generated sequence')
	arg = parser.parse_args()
	
	
	dna = organize_dna(arg.dnafile)
	kmers = make_kmerdict(arg.markov_order,dna)
	pool = make_contextpool(kmers)
	seq = generate_seq(arg.markov_order,pool,arg.seqlen)
	print(seq)

'''
Test cod

dna = 'TTCAGAATGCTACCGGTA'
UNIT TEST
AA: T
AC: C
AG: A
AT: G
CC: G
CA: G
CG: T
CT: A
GG: T
GA: A
GC: T
GT: A
TT: C
TA: C
TC: A
TG: C
'''	
	
	

