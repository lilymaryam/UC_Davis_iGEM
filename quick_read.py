#!/usr/bin/env python3

import motiflib
import argparse
import weeder_reader

parser = argparse.ArgumentParser(
	description='Compare meme output to weeder output.')
# required arguments
parser.add_argument('--memeout', required=True, type=str,
	metavar='<str>', help='meme output')
parser.add_argument('--weederout', required=True, type=str,
	metavar='<str>', help='weeder output')
parser.add_argument('--jasparfile', required=True, type=str,
	metavar='<str>', help='jaspar file')
#parser.add_argument('--rfloat', required=True, type=float,
#	metavar='<float>', help='required floating point argument')
# optional arguments with default parameters
#parser.add_argument('--dstr', required=False, type=str, default='hello',
#	metavar='<str>', help='optional string argument [%(default)s]')
#parser.add_argument('--dint', required=False, type=int, default=1,
#	metavar='<int>', help='optional integer argument [%(default)i]')
#parser.add_argument('--dfloat', required=False, type=float, default=3.14,
#	metavar='<float>', help='optional floating point argument [%(default)f]')
# switches
#parser.add_argument('--switch', action='store_true',
#	help='on/off switch')
# finalization
arg = parser.parse_args()

#print(arg.memeout)
#print(arg.weederout)
#print(arg.jasparfile)

background = {'A': .25,'C': .25,'G': .25, 'T':.25}

#with open(arg.memeout) as mo:
motifs, motif_info = motiflib.memepwm(arg.memeout)
#print(motifs)
#print(motif_info)
#with open(arg.jasparfile) as jf:
jpwm = motiflib.read_JASPAR(arg.jasparfile)
#print(jpwm)
weederpwm = weeder_reader.weederpwm(arg.weederout)
print(len(weederpwm))

for i in range(len(weederpwm)):
	lscore = motiflib.local_motcompare(weederpwm[i],jpwm)
	gscore = motiflib.global_motcompare(weederpwm[i],jpwm,background)
	print(f'MAT{i+1}',gscore,lscore)
	