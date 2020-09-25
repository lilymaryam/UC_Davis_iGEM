#!/usr/bin/python3

import motiflib
import sys
import os 
import argparse
import statistics
#import markov_memepipe

parser = argparse.ArgumentParser(
	description='Compare meme performance to motif size and info')
# required arguments
parser.add_argument('--memepath', required=True, type=str,
	metavar='<str>', help='path to meme software')
parser.add_argument('--dnafile', required=True, type=str,
	metavar='<str>', help='file containing background dna')
parser.add_argument('--promlength', required=True, type=int,
	metavar='<int>', help='length of promoter')
parser.add_argument('--numseq', required=True, type=int,
	metavar='<int>', help='number of promoters')
parser.add_argument('--markov_order', required=True, type=int,
	metavar='<int>', help='markov_order')
parser.add_argument('--numiterations', required=False, type=int, default=5,
	metavar='<int>', help='iterations for condensed data [%(default)i]')
parser.add_argument('--motiffrequency', required=False, type=float, 
	default=0.9, metavar='<float>', help='background probability of A \
	[%(default).3f]')
#parser.add_argument('--PA', required=False, type=float, default=0.25,
#	metavar='<float>', help='background probability of A [%(default).3f]')
#parser.add_argument('--PC', required=False, type=float, default=0.25,
#	metavar='<float>', help='background probability of C [%(default).3f]')
#parser.add_argument('--PG', required=False, type=float, default=0.25,
#	metavar='<float>', help='background probability of G [%(default).3f]')
#parser.add_argument('--PT', required=False, type=float, default=0.25,
#	metavar='<float>', help='background probability of T [%(default).3f]')
'''
parser.add_argument('--bothstrands', action='store_true',
	help='on/off switch')
parser.add_argument('--negstrands', action='store_true',
	help='on/off switch')
parser.add_argument('--condenseddata', action='store_true',
	help='on/off switch')
'''
'''
parser.add_argument('--rfloat', required=True, type=float,
	metavar='<float>', help='required floating point argument')
# optional arguments with default parameters
parser.add_argument('--dstr', required=False, type=str, default='hello',
	metavar='<str>', help='optional string argument [%(default)s]')
parser.add_argument('--dint', required=False, type=int, default=1,
	metavar='<int>', help='optional integer argument [%(default)i]')
parser.add_argument('--dfloat', required=False, type=float, default=3.14,
	metavar='<float>', help='optional floating point argument [%(default)f]')
# switches
parser.add_argument('--switch', action='store_true',
	help='on/off switch')
# finalization
'''
arg = parser.parse_args()
#assert(len(sys.argv) == 2)
#jasparfile = sys.argv[1]

#motif = motiflib.read_JASPAR(jasparfile)
o = arg.markov_order
n = arg.numseq
p = arg.promlength
freq = arg.motiffrequency
directory = 'Jaspar'
model = ['zoops','oops','anr']


jasparfiles = os.listdir(directory)
jasparfiles.sort()
print('file','length motif','information in bits','percentage bit score','performance','performance percentage','model',sep=', ')
for fl in jasparfiles:
	filepath = f'Jaspar/{fl}'
	motif = motiflib.read_JASPAR(filepath)
	bits,p_bits = motiflib.score_motifbit(motif)
	for m in model:
		modscore = []
		modscorep = []
		for i in range(arg.numiterations):
			tmpfile = f'/tmp/testmotif{fl}_{i}.fa' 
			cmd = f'python3 markov_motifapalooza.py --jasparfile {filepath} \
			--dnafile {arg.dnafile} --markov_order {o} --numseq {n} --seqlen {p} --freq {freq} > {tmpfile} '
			os.system(cmd)
			cmd = f'{arg.memepath} {tmpfile} -dna -markov_order {o} -mod {m}'
			os.system(cmd)
			memetxt = 'meme_out/meme.txt'
			motifs,motif_stats = motiflib.memepwm(memetxt)
			for mot in motifs:
				performance = motiflib.local_motcompare(mot,motif)
				modscore.append(performance[0])
				modscorep.append(performance[1])
		print(tmpfile, len(motif), bits,p_bits,statistics.mean(modscore),statistics.mean(modscorep),m,sep=',')
	#print(filepath, len(motif),bits,p_bits, sep=',')
	
