  
#!/usr/bin/env python3

import sys
import gzip
import re
import argparse

parser = argparse.ArgumentParser(
	description='Extracts Intergenic Sequences from Gff File')
parser.add_argument('--gff', required=True, type=str,
	metavar='<str>', help='Gff File')
#parser.add_argument('--upstream', required=False, type=int, default=500,
#	metavar='<int>', help='length of promoter region [%(default)i]')
parser.add_argument('--sequence', required=False, type=str,
        metavar='<str>', help='Sequence file associated with Gff')
parser.add_argument('--index', required=True, type=str,
        metavar='<str>', help='file containing gene IDs of interest')
parser.add_argument('--upstream', required=False, type=int, default=500,
                    metavar='<int>', help='length of promoter region [%(default)i]')
arg = parser.parse_args()
chrom = {}
gene = {}
beg = 0
end = 0
strand = None
seqs = []
next1 = 0
promoter = ''
x=0
with open(arg.index, 'rt') as fp:
	with open(arg.gff, 'rt') as gf:
		with open(arg.sequence, 'rt') as seq:
			#creates a dictionary indexed by chromosome number with sequences contained within
			for line in seq:
				if line[0] == '>' and line[1] == 'C':
					chr = re.search('(>Chr\w+)', line)
					chrind = chr.group()
					print(chrind)
					chrom[chrind] = ''
				chrom[chrind] += line.strip()

			while True:
				#reads a file of orthologous cluster gene ID's
				lineindex = fp.readline()
				print(lineindex)
				geneid = re.search('(\D+)(\d+)\D+(\d+)', lineindex)
				if geneid == None:
					break
				spec, beg, end = geneid.groups()
				next1 = int(beg)

				while int(next1)<=int(end):
					linegff = gf.readline()
					strand = re.search('.	(\W)	.	(ID=AN%s-T-E)' %next1, linegff)

					if strand:
						direc, blank = strand.groups()
						coordinates = re.search('\s+(\d+)\s+\d+', linegff)
						chromosome = re.search('Chr\w+', linegff)
						chromosome = '>' + chromosome.group()
						start =(coordinates.group(1))
						next1 = next1 + 1
						promstart = int(start) - arg.upstream
						header = chromosome + '  ' + blank + '  ' 'strand =  ' + direc

						while promstart <= int(start):
							promoter = promoter + chrom[chromosome][promstart]
							promstart = promstart + 1
						#output here
						print(header)
						print(promoter)