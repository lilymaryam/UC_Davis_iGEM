#!/usr/bin/env python3

import sys
import gzip
import re
import argparse

parser = argparse.ArgumentParser(
	description='Extracts promoters from GenBank file.')
parser.add_argument('--genbank', required=True, type=str,
	metavar='<str>', help='genbank file')
parser.add_argument('--upstream', required=False, type=int, default=500,
	metavar='<int>', help='length of promoter region [%(default)i]')
arg = parser.parse_args()

gene = {}
beg = None
end = None
strand = None
seqs = []
with gzip.open(arg.genbank, 'rt') as fp:
	while True:
		line = fp.readline()
		if line.startswith('ORIGIN'): break
		if line.startswith('     gene'):
			m1 = re.search('<?(\d+)\.\.>?(\d+)', line)
			if m1:
				beg, end = m1.groups()
				m2 = re.search('complement', line)
				if m2: strand = '+'
				else:  strand = '-'
		if line.startswith('                     /gene'):
			m = re.search('/gene="(\S+)"', line)
			if m:
				name = m.groups()[0]
				if name not in gene:
					gene[name] = (int(beg), int(end), strand)
		if line.startswith('                     /locus_tag'):
			m = re.search('/locus_tag="(\S+)"', line)
			if m:
				name = m.groups()[0]
				if name not in gene:
					gene[name] = (int(beg), int(end), strand)
	
	while True:
		line = fp.readline()
		if line.startswith('//'): break
		f = line.split()
		s = ''.join(f[1:]).upper()
		seqs.append(s)

seq = ''.join(seqs)
for name in gene:
	beg, end, strand = gene[name]
	
	c1, c2 = None, None
	if strand == '+':
		c1 = beg - arg.upstream
		if c1 < 0: c1 = 0
		c2 = beg
	else:
		c1 = end
		c2 = end + arg.upstream
		if c2 > len(seq): c2 = len(seq)
	
	if c2 - c1 < arg.upstream: continue
	
	print(f'>{name} {arg.upstream} {strand}')
	print(f'{seq[c1:c2]}')
