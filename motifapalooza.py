#!/usr/bin/env python3
import argparse
import random

# setup
parser = argparse.ArgumentParser(
	description='Generates random sequence and embeds known fungal TF binding motifs.')
# required arguments
parser.add_argument('--jasparfile', required=True, type=str,
	metavar='<str>', help='jaspar file')
parser.add_argument('--numseq', required=False, type=int, default=10,
	metavar='<int>', help='number of sequences to generate [%(default)i]')
parser.add_argument('--seqlen', required=False, type=int, default=100,
	metavar='<int>', help='length of sequences to generate [%(default)i]')
parser.add_argument('--mps', required=False, type=int, default=1,
	metavar='<int>', help='number of motifs per seqeunce [%(default)i]')
parser.add_argument('--freq', required=False, type=float, default=0.9,
	metavar='<float>', help='optional floating point argument [%(default)f]')	
parser.add_argument('--PA', required=False, type=float, default=0.25,
	metavar='<float>', help='optional floating point argument [%(default)f]')
parser.add_argument('--PC', required=False, type=float, default=0.25,
	metavar='<float>', help='optional floating point argument [%(default)f]')
parser.add_argument('--PG', required=False, type=float, default=0.25,
	metavar='<float>', help='optional floating point argument [%(default)f]')
parser.add_argument('--PT', required=False, type=float, default=0.25,
	metavar='<float>', help='optional floating point argument [%(default)f]')
# switches
parser.add_argument('--bothstrands', action='store_true',
	help='on/off switch')
# finalization
arg = parser.parse_args()

def read_JASPAR(jasparfile):
	motif = []
	linenum = 0
	nt = ['A', 'C', 'G', 'T']
	with open(arg.jasparfile) as jf:
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
				
	
def generate_seq(numseq, seqlen, PA, PC, PG, PT):
	seq = []
	for i in range(0,arg.seqlen):
		r  = random.random()
		if r < PA:             seq += 'a'
		elif r < PA + PC:      seq += 'c'
		elif r < PA + PC + PG: seq += 'g'
		else:                  seq += 't'
	return seq
		
		
motif = read_JASPAR(arg.jasparfile)
for i in range(0, arg.numseq):
	seq = (generate_seq(arg.numseq, arg.seqlen, arg.PA, arg.PC, arg.PG, arg.PT))
	r = random.random()
	places = []
	if r < arg.freq:
		for j in range(0,arg.mps):
			site = generate_site(motif)
			assert(len(motif)== len(site))
			place = random.randint(0, len(seq)-len(motif))
			places.append(f'{place}+')
			for k in range(0,len(site)):
				seq[place+k] = site[k]
	print(f'>seq-{i} {places}')
	print(''.join(seq))

	
		
		
		
	


