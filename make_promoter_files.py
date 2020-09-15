#!/usr/bin/python3

import argparse
import findss
import statistics

parser = argparse.ArgumentParser(
	description='Takes a JGI file and converts it into expected form.')
# required arguments
parser.add_argument('--jgi_fasta', required=True, type=str,
	metavar='<str>', help='required string argument')
parser.add_argument('--jgi_genecoords', required=True, type=str,
	metavar='<str>', help='text file containing transcription and translation start sites')
parser.add_argument('--promoter_size', required=True, type=int,
	metavar='<int>', help='number of base pairs in the promoter region')
# optional arguments with default parameters
'''
parser.add_argument('--dstr', required=False, type=str, default='hello',
	metavar='<str>', help='optional string argument [%(default)s]')
parser.add_argument('--dint', required=False, type=int, default=1,
	metavar='<int>', help='optional integer argument [%(default)i]')
parser.add_argument('--dfloat', required=False, type=float, default=3.14,
	metavar='<float>', help='optional floating point argument [%(default)f]')
'''
# switches
'''
parser.add_argument('--switch', action='store_true',
	help='on/off switch')
# finalization
'''
arg = parser.parse_args()

p5,p3,n5,n3 = findss.return_distances(arg.jgi_genecoords)
histograms = findss.create_histograms(p5,p3,n5,n3)
avgp5,avgp3,avgn5,avgn3 = findss.find_histstats(p5,p3,n5,n3)

distance = round(statistics.mean([avgp5,avgn5]))
print(distance)

#distance = 10
with open(arg.jgi_fasta) as jf:
	with open('file.txt','w') as fl:
		while True:
			line = jf.readline()
			if line == '': break
			if line.startswith('>'):
				fl.write(line)
			else:
				fl.write(line[distance:(distance+arg.promoter_size)]+'\n')

				
with open('file.txt') as fl:
	for line in fl.readlines():
		line = line.strip()
		print(line)
		