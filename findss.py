#!/usr/bin/env python3
import sys
import gzip

with gzip.open(sys.argv[1]) as f:
	for line in f.readlines():
		line = line.strip()
		print(line)
		line = line.split()
		print(line)
		line = line.decode(line)
		strand = line[3]
		print(type(strand))
		strand = strand.strip()
		print(strand)
		#if strand[1] == '+': print(yes)
		
		
	
