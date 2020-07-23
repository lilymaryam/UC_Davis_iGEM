#!/usr/bin/env python3

import argparse
import motiflib


parser = argparse.ArgumentParser(
	description='Do motif instances look like motifs?')
parser.add_argument('--file', required=True, type=str,
	metavar='<str>', help='motif file in JASPAR format')
parser.add_argument('--seqs', required=False, type=int, default=5,
	metavar='<int>', help='sequences to generate [%(default)i]')
arg = parser.parse_args()

motif = motiflib.read_JASPAR(arg.file)
for i in range(arg.seqs):
	site = motiflib.generate_site(motif)
	print(site)

# if you collected the sites back into a motif, would it look like the motif?
# how would you measure the similarity?
# you might want to run this many times to look at the variance
