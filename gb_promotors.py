#!/usr/bin/env python3

import sys
import gzip
import re

gene = {}
beg = None
end = None
strand = None
seqs = []
with gzip.open(sys.argv[1], 'rt') as fp:
	while True:
		line = fp.readline()
		if line.startswith('ORIGIN'): break
		if line.startswith('     gene'):
			m1 = re.search('(\d+)\.\.(\d+)', line)
			if m1:
				beg, end = m1.groups()
				m2 = re.search('complement', line)
				if m2: strand = '+'
				else:  strand = '-'
		if line.startswith('                     /gene'):
			m = re.search('gene="(\S+)"', line)
			if m:
				name = m.groups()[0]
				if name not in gene:
					gene[name] = (beg, end, strand)
	
	while True:
		line = fp.readline()
		if line.startswith('//'): break
		f = line.split()
		s = ''.join(f[1:]).upper()
		seqs.append(s)

seq = ''.join(seqs)
for name in gene:
	beg, end, stand = gene[name]
	print(name, beg, end, strand)
	# now extract the promoter and make fasta file

