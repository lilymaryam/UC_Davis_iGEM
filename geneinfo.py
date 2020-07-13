#!/usr/bin/env python3

import gzip
import sys

genome = {}
with gzip.open('gene_coordinates.txt.gz', 'rt') as fp:
	for line in fp.readlines():
		species, gid, scf, st, gb, ge, cb, ce, fs, s1, s2 = line.split()
		if species not in genome: genome[species] = []
		d = None
		if st == '+': d = int(cb) - int(gb)
		else:         d = int(ge) - int(ce)
		genome[species].append(d)

for species in genome:
	with open(f'{species}.utr5.txt', 'w') as fp:
		for length in genome[species]:
			fp.write(f'{length}\n')


