import random
import sys
import argparse

parser = argparse.ArgumentParser(
	description='Downsample from an existing promoter file')

parser.add_argument('--size', required=True, type=int, default = 1,
                    metavar='<int>', help='Sample Size [%(default)i]')

parser.add_argument('--fasta', required=True, type=str,
	metavar='<str>', help='fasta file of promoters')

parser.add_argument('--iterate', required=False, type=int,
                    metavar='<int>', help='Iterates by this number from --size -> original fasta size')

arg = parser.parse_args()
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

