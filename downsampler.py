import random
import sys
import argparse
#import motif

parser = argparse.ArgumentParser(
	description='Downsample from an existing promoter file')
parser.add_argument('--fasta', required=True, type=str,
	metavar='<str>', help='fasta file of promoters')
parser.add_argument('--step', required=False, type=int, default=2,
	metavar='<int>', help='step down value [%(default)i]')
parser.add_argument('--floor', required=False, type=int, default=2,
	metavar='<int>', help='step down value [%(default)i]')
parser.add_argument('--samples', required=False, type=int, default=3,
    metavar='<int>', help='number of times at each step [%(default)i]')
arg = parser.parse_args()

## Part 1: run meme with whole dataset
# pwm1 = motif.run_meme(...)

## Part 2: run meme with downsampled datasets, compare to whole run
seqs = []
with open(arg.fasta) as fp:
	defline = None
	for line in fp.readlines():
		if line.startswith('>'): defline = line
		else: seqs.append((defline, line))

subtract = 0
while True:
	subtract += arg.step
	size = len(seqs) - subtract
	if size < arg.floor: break
	#print(f'data set size = {size}')
	for i in range(arg.samples):
		#print(f'set {i}')
		set = random.sample(seqs, size)
		
		# write temp file
		with open(f'downsampler.temp.{size}.{i}.fasta', 'w') as fp:
			for s in set:
				fp.write(s[0])
				fp.write(s[1])
		
		# run meme
		# pwm2 = motif.run_meme(...)
		
		# compare meme output to the whole set
		#if motif.distance(pwm1, pwm2) < threshold:...
	


"""
promnumb = 0
outseq = ''
match = 0
count = 0
up = arg.iterate

with open(arg.fasta, 'rt') as fast:
    for l in fast:
           if l[0] == ">":
               count += 1
with open(arg.fasta, 'rt') as fast:

    if arg.iterate:
        while arg.size + up <= count:
            sample = random.sample(range(0, count), arg.size + up)
            for line in fast:
                if line[0] == ">":
                    match = 0
                    header = line
                    for s in sample:
                        if promnumb == s:
                            match = 1
                            outseq += header
                    promnumb += 1
                if match == 1:
                    outseq += line.strip(header) + '\n'
            up += arg.iterate
            print(outseq)

    else:

        sample = random.sample(range(0, count), arg.size)

        for line in fast:

            if line[0] == ">":
                match = 0
                header = line

                for s in sample:

                    if promnumb == s:
                        match = 1
                        outseq += header
                promnumb += 1

            if match == 1:
                outseq += line.strip(header) + '\n'
        print(outseq)
"""
